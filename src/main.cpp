#include <chrono>
#include <iostream>
#include <iomanip>
#include <thread>
#include <random>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/process/context.hpp>
#include <boost/process/operations.hpp>
#include <boost/process/child.hpp>
#include "ligand.hpp"
#include "io_service_pool.hpp"
#include "operation.hpp"
using namespace boost;
using namespace boost::filesystem;

//! Represents a thread safe counter.
template <typename T>
class safe_counter
{
public:
	//! Initializes the counter to 0 and its expected hit value to z.
	void init(const T z);

	//! Increments the counter by 1 in a thread safe manner, and wakes up the calling thread waiting on the internal mutex.
	void increment();

	//! Waits until the counter reaches its expected hit value.
	void wait();
private:
	mutex m;
	condition_variable cv;
	T n; //!< Expected hit value.
	T i; //!< Counter value.
};

template <typename T>
void safe_counter<T>::init(const T z)
{
	n = z;
	i = 0;
}

template <typename T>
void safe_counter<T>::increment()
{
	lock_guard<mutex> guard(m);
	if (++i == n) cv.notify_one();
}

template <typename T>
void safe_counter<T>::wait()
{
	unique_lock<mutex> lock(m);
	if (i < n) cv.wait(lock);
}

template class safe_counter<size_t>;

int main(int argc, char* argv[])
{
	// Initialize the default path to log files. They will be reused when calling idock.
	const path default_log_path = "log.txt";
	const path default_csv_path = "log.csv";

	path initial_generation_csv_path, initial_generation_folder_path, fragment_folder_path, idock_config_path, output_folder_path, csv_path;
	size_t num_threads, seed, num_elitists, num_additions, num_subtractions, num_crossovers, max_failures, max_rotatable_bonds, max_atoms, max_heavy_atoms, max_hb_donors, max_hb_acceptors;
	fl max_mw;

	// Process program options.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_additions = 20;
		const size_t default_num_subtractions = 20;
		const size_t default_num_crossovers = 20;
		const size_t default_num_elitists = 10;
		const size_t default_max_failures = 1000;
		const size_t default_max_rotatable_bonds = 30;
		const size_t default_max_atoms = 100;
		const size_t default_max_heavy_atoms = 80;
		const size_t default_max_hb_donors = 5;
		const size_t default_max_hb_acceptors = 10;
		const fl default_max_mw = 500;

		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("initial_generation_csv", value<path>(&initial_generation_csv_path)->required(), "path to initial generation csv")
			("initial_generation_folder", value<path>(&initial_generation_folder_path)->required(), "path to initial generation folder")
			("fragment_folder", value<path>(&fragment_folder_path)->required(), "path to folder of fragments in PDBQT format")
			("idock_config", value<path>(&idock_config_path)->required(), "path to idock configuration file")
			;

		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output results")
			("csv", value<path>(&csv_path)->default_value(default_csv_path), "summary file in csv format")
			;

		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("elitists", value<size_t>(&num_elitists)->default_value(default_num_elitists), "number of elite ligands to carry over")
			("additions", value<size_t>(&num_additions)->default_value(default_num_additions), "number of child ligands created by addition")
			("subtractions", value<size_t>(&num_subtractions)->default_value(default_num_subtractions), "number of child ligands created by subtraction")
			("crossovers", value<size_t>(&num_crossovers)->default_value(default_num_crossovers), "number of child ligands created by crossover")
			("max_failures", value<size_t>(&max_failures)->default_value(default_max_failures), "maximum number of operational failures to tolerate")
			("max_rotatable_bonds", value<size_t>(&max_rotatable_bonds)->default_value(default_max_rotatable_bonds), "maximum number of rotatable bonds")
			("max_atoms", value<size_t>(&max_atoms)->default_value(default_max_atoms), "maximum number of atoms")
			("max_heavy_atoms", value<size_t>(&max_heavy_atoms)->default_value(default_max_heavy_atoms), "maximum number of heavy atoms")
			("max_hb_donors", value<size_t>(&max_hb_donors)->default_value(default_max_hb_donors), "maximum number of hydrogen bond donors")
			("max_hb_acceptors", value<size_t>(&max_hb_acceptors)->default_value(default_max_hb_acceptors), "maximum number of hydrogen bond acceptors")
			("max_mw", value<fl>(&max_mw)->default_value(default_max_mw), "maximum molecular weight")
			("help", "help information")
			("version", "version information")
			("config", value<path>(), "options can be loaded from a configuration file")
			;

		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// If no command line argument is supplied, simply print the usage and exit.
		if (argc == 1)
		{
			cout << all_options;
			return 0;
		}

		// Parse command line arguments.
		variables_map vm;
		store(parse_command_line(argc, argv, all_options), vm);

		// If no command line argument is supplied or help is requested, print the usage and exit.
		if (argc == 1 || vm.count("help"))
		{
			cout << all_options;
			return 0;
		}

		// If version is requested, print the version and exit.
		if (vm.count("version"))
		{
			cout << "1.0.0" << endl;
			return 0;
		}

		// If a configuration file is presented, parse it.
		if (vm.count("config"))
		{
			boost::filesystem::ifstream config_file(vm["config"].as<path>());
			store(parse_config_file(config_file, all_options), vm);
		}

		// Notify the user of parsing errors, if any.
		vm.notify();

		using namespace boost::filesystem;

		// Validate initial generation csv.
		if (!exists(initial_generation_csv_path))
		{
			cerr << "Initial generation csv " << initial_generation_csv_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(initial_generation_csv_path))
		{
			cerr << "Initial generation csv " << initial_generation_csv_path << " is not a regular file\n";
			return 1;
		}

		// Validate initial generation folder.
		if (!exists(initial_generation_folder_path))
		{
			cerr << "Initial generation folder " << initial_generation_folder_path << " does not exist\n";
			return 1;
		}
		if (!is_directory(initial_generation_folder_path))
		{
			cerr << "Initial generation folder " << initial_generation_folder_path << " is not a directory\n";
			return 1;
		}

		// Validate fragment folder.
		if (!exists(fragment_folder_path))
		{
			cerr << "Fragment folder " << fragment_folder_path << " does not exist\n";
			return 1;
		}
		if (!is_directory(fragment_folder_path))
		{
			cerr << "Fragment folder " << fragment_folder_path << " is not a directory\n";
			return 1;
		}

		// Validate idock configuration file.
		if (!exists(idock_config_path))
		{
			cerr << "idock configuration file " << idock_config_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(idock_config_path))
		{
			cerr << "idock configuration file " << idock_config_path << " is not a regular file\n";
			return 1;
		}

		// Validate output folder.
		remove_all(output_folder_path);
		if (!create_directories(output_folder_path))
		{
			cerr << "Failed to create output folder " << output_folder_path << endl;
			return 1;
		}

		// Validate csv_path.
		if (is_directory(csv_path))
		{
			cerr << "csv path " << csv_path << " is a directory\n";
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			cerr << "Option threads must be 1 or greater\n";
			return 1;
		}
		if (max_mw <= 0)
		{
			cerr << "Option max_mw must be positive\n";
			return 1;
		}
	}
	catch (const std::exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	// The number of ligands (i.e. population size) is equal to the number of elitists plus mutants plus children.
	const size_t num_children = num_additions + num_subtractions + num_crossovers;
	const size_t num_ligands = num_elitists + num_children;
	const fl num_elitists_inv = static_cast<fl>(1) / num_elitists;

	// Initialize a pointer vector to dynamically hold and destroy generated ligands.
	boost::ptr_vector<ligand> ligands;
	ligands.resize(num_ligands);

	// Parse the initial generation csv to get initial elite ligands.
	{
		boost::filesystem::ifstream in(initial_generation_csv_path);
		string line;
		line.reserve(80);
		getline(in, line); // Ligand,Conf,FE1,FE2,FE3,FE4,FE5,FE6,FE7,FE8,FE9
		for (size_t i = 0; i < num_elitists; ++i)
		{
			// Check if there are sufficient initial elite ligands.
			if (!getline(in, line))
			{
				cerr << "Failed to construct initial generation because the initial generation csv " << initial_generation_csv_path << " contains less than " << num_elitists << " ligands.\n";
				return 1;
			}

			// Parse the elite ligand.
			const size_t comma1 = line.find(',', 1);
			ligands.replace(i, new ligand(initial_generation_folder_path / (line.substr(0, comma1) + ".pdbqt")));

			// Parse the free energy and ligand efficiency.
			const size_t comma2 = line.find(',', comma1 + 2);
			const size_t comma3 = line.find(',', comma2 + 6);
			const size_t comma4 = line.find(',', comma3 + 6);
			ligands[i].fe = lexical_cast<fl>(line.substr(comma2 + 1, comma3 - comma2 - 1));
			ligands[i].le = lexical_cast<fl>(line.substr(comma3 + 1, comma4 - comma3 - 1));
		}
	}

	// Scan the fragment folder to obtain a list of fragments.
	cout << "Scanning fragment folder " << fragment_folder_path << endl;
	vector<path> fragments;
	fragments.reserve(1000); // A fragment folder typically consists of <= 1000 fragments.
	for (directory_iterator dir_iter(fragment_folder_path), end_dir_iter; dir_iter != end_dir_iter; ++dir_iter)
	{
		// Skip non-regular files such as folders.
		if (!is_regular_file(dir_iter->status())) continue;
		// Save the fragment path.
		fragments.push_back(dir_iter->path());
	}
	cout << "Found " << fragments.size() << " fragments\n";

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << endl;
	mt19937_64 eng(seed);

	// Initialize a ligand validator.
	const validator v(max_rotatable_bonds, max_atoms, max_heavy_atoms, max_hb_donors, max_hb_acceptors, max_mw);

	// Initialize the number of failures. The program will stop if num_failures reaches max_failures.
	atomic<size_t> num_failures(0);

	// Reserve storage for operation tasks.
	operation op(ligands, num_elitists, fragments, v, max_failures, num_failures);

	// Initialize ligand filenames.
	vector<string> ligand_filenames;
	ligand_filenames.reserve(num_ligands);
	for (size_t i = 1; i <= num_ligands; ++i)
	{
		ligand_filenames.push_back(lexical_cast<string>(i) + ".pdbqt");
	}

	// Find the full path to idock executable.
	const path idock_path = path(boost::process::find_executable_in_path("idock")).make_preferred();
	cout << "Using idock executable at " << idock_path << endl;

	// Initialize arguments to idock.
	vector<string> idock_args(12);
	idock_args[0]  = "--ligand_folder";
	idock_args[2]  = "--output_folder";
	idock_args[4]  = "--log";
	idock_args[6]  = "--csv";
	idock_args[8]  = "--seed";
	idock_args[9]  = lexical_cast<string>(seed);
	idock_args[10] = "--config";
	idock_args[11] = idock_config_path.string();

	// Initialize process context.
	const boost::process::context ctx;

	// Initialize an io service pool and create worker threads for later use.
	cout << "Creating an io service pool of " << num_threads << " worker thread" << ((num_threads == 1) ? "" : "s") << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize csv file for dumping statistics.
	boost::filesystem::ofstream csv(csv_path);
	csv << "generation,ligand,parent 1,connector 1,parent 2,connector 2,free energy in kcal/mol,ligand efficiency in kcal/mol,no. of rotatable bonds,no. of atoms,no. of heavy atoms,no. of hydrogen bond donors,no. of hydrogen bond acceptors,molecular weight in g/mol\n";

	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(3);
	for (size_t generation = 1; true; ++generation)
	{
		cout << "Running generation " << generation << endl;

		// Initialize the paths to current generation folder and its two subfolders.
		const path generation_folder(output_folder_path / lexical_cast<string>(generation));
		const path ligand_folder(generation_folder / "ligand");
		const path output_folder(generation_folder / "output");

		// Create a new folder and two subfolders for current generation.
		create_directory(generation_folder);
		create_directory(ligand_folder);
		create_directory(output_folder);

		// Create addition, subtraction and crossover tasks.
		cnt.init(num_children);
		for (size_t i = 0; i < num_additions; ++i)
		{
			io.post([&,i]()
			{
				op.addition_task(num_elitists + i, ligand_folder / ligand_filenames[i], eng());
				cnt.increment();
			});
		}
		for (size_t i = num_additions; i < num_additions + num_subtractions; ++i)
		{
			io.post([&,i]()
			{
				op.subtraction_task(num_elitists + i, ligand_folder / ligand_filenames[i], eng());
				cnt.increment();
			});
		}
		for (size_t i = num_additions + num_subtractions; i < num_children; ++i)
		{
			io.post([&,i]()
			{
				op.crossover_task(num_elitists + i, ligand_folder / ligand_filenames[i], eng());
				cnt.increment();
			});
		}
		cnt.wait();

		// Check if the maximum number of failures has been reached.
		if (num_failures >= max_failures)
		{
			cout << "The number of failures has reached " << max_failures << endl;
			return 0;
		}

		// Invoke idock.
		idock_args[1] = ligand_folder.string();
		idock_args[3] = output_folder.string();
		idock_args[5] = (generation_folder / default_log_path).string();
		idock_args[7] = (generation_folder / default_csv_path).string();
		const int exit_code = create_child(idock_path.string(), idock_args, ctx).wait();
		if (exit_code)
		{
			cout << "idock exited with code " << exit_code << endl;
			return 1;
		}

		// Parse docked ligands to obtain predicted free energy and docked coordinates, and save the updated ligands into the ligand subfolder.
		for (size_t i = 0; i < num_children; ++i)
		{
			ligands[num_elitists + i].update(output_folder / ligand_filenames[i]);
		}

		// Sort ligands in ascending order of efficacy.
		ligands.sort();

		// Write summaries to csv and calculate average statistics.
		for (size_t i = 0; i < num_ligands; ++i)
		{
			const ligand& l = ligands[i];
			csv << generation
				<< ',' << l.p
				<< ',' << l.parent1
				<< ',' << l.connector1
				<< ',' << l.parent2
				<< ',' << l.connector2
				<< ',' << l.fe
				<< ',' << l.le
				<< ',' << l.num_rotatable_bonds
				<< ',' << l.num_atoms
				<< ',' << l.num_heavy_atoms
				<< ',' << l.num_hb_donors
				<< ',' << l.num_hb_acceptors
				<< ',' << l.mw
				<< endl;
		}

		// Calculate average statistics of elite ligands.
		fl avg_mw = 0, avg_fe = 0, avg_le = 0, avg_rotatable_bonds = 0, avg_atoms = 0, avg_heavy_atoms = 0, avg_hb_donors = 0, avg_hb_acceptors = 0;
		for (size_t i = 0; i < num_elitists; ++i)
		{
			const ligand& l = ligands[i];
			avg_mw += l.mw;
			avg_fe += l.fe;
			avg_le += l.le;
			avg_rotatable_bonds += l.num_rotatable_bonds;
			avg_atoms += l.num_atoms;
			avg_heavy_atoms += l.num_heavy_atoms;
			avg_hb_donors += l.num_hb_donors;
			avg_hb_acceptors += l.num_hb_acceptors;
		}
		avg_mw *= num_elitists_inv;
		avg_fe *= num_elitists_inv;
		avg_le *= num_elitists_inv;
		avg_rotatable_bonds *= num_elitists_inv;
		avg_atoms *= num_elitists_inv;
		avg_heavy_atoms *= num_elitists_inv;
		avg_hb_donors *= num_elitists_inv;
		avg_hb_acceptors *= num_elitists_inv;
		cout << "Failures |  Avg FE |  Avg LE |  Avg HA | Avg MWT | Avg NRB | Avg HBD | Avg HBA\n"
		    << setw(8) << num_failures << "   "
			<< setw(7) << avg_fe << "   "
			<< setw(7) << avg_le << "   "
			<< setw(7) << avg_heavy_atoms << "   "
			<< setw(7) << avg_mw << "   "
			<< setw(7) << avg_rotatable_bonds << "   "
			<< setw(7) << avg_hb_donors << "   "
			<< setw(7) << avg_hb_acceptors << endl;
	}
}
