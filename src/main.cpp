#include <chrono>
#include <iostream>
#include <iomanip>
#include <random>
#include <boost/program_options.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/process.hpp>
#include "io_service_pool.hpp"
#include "safe_counter.hpp"
#include "ligand.hpp"
using namespace std::chrono;
using namespace boost::filesystem;
using namespace boost::process;
using namespace boost::process::initializers;

int main(int argc, char* argv[])
{
	// Declare program option variables.
	path idock_example_path, out_path;
	size_t seed, num_threads, num_elitists, num_children, num_generations, nrb_lb, nrb_ub, hbd_lb, hbd_ub, hba_lb, hba_ub;
	double mms_lb, mms_ub;

	// Process program options.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_out_path = ".";
		const size_t default_seed = duration_cast<seconds>(system_clock::now().time_since_epoch()).count();
		const size_t default_num_threads = boost::thread::hardware_concurrency();
		const size_t default_num_children = 20;
		const size_t default_num_elitists = 10;
		const size_t default_num_generations = 8;
		const size_t default_nrb_lb =  0;
		const size_t default_nrb_ub = 10;
		const size_t default_hbd_lb =  0;
		const size_t default_hbd_ub =  5;
		const size_t default_hba_lb =  0;
		const size_t default_hba_ub = 10;
		const double default_mms_lb = 300;
		const double default_mms_ub = 500;

		// Set up the description of options.
		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("idock_example", value<path>(&idock_example_path)->required(), "path to an idock example folder")
			;
		options_description output_options("output (optional)");
		output_options.add_options()
			("out", value<path>(&out_path)->default_value(default_out_path), "folder of generated ligands")
			;
		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "# of worker threads to use")
			("elitists", value<size_t>(&num_elitists)->default_value(default_num_elitists), "# of elite ligands to carry over")
			("children", value<size_t>(&num_children)->default_value(default_num_children), "# of child ligands created from elite ligands")
			("generations", value<size_t>(&num_generations)->default_value(default_num_generations), "# of GA generations to run")
			("mms_lb", value<double>(&mms_lb)->default_value(default_mms_lb), "minimum molecular mass in Dalton")
			("mms_ub", value<double>(&mms_ub)->default_value(default_mms_ub), "maximum molecular mass in Dalton")
			("nrb_lb", value<size_t>(&nrb_lb)->default_value(default_nrb_lb), "minimum # of rotatable bonds")
			("nrb_ub", value<size_t>(&nrb_ub)->default_value(default_nrb_ub), "maximum # of rotatable bonds")
			("hbd_lb", value<size_t>(&hbd_lb)->default_value(default_hbd_lb), "minimum # of hydrogen bond donors")
			("hbd_ub", value<size_t>(&hbd_ub)->default_value(default_hbd_ub), "maximum # of hydrogen bond donors")
			("hba_lb", value<size_t>(&hba_lb)->default_value(default_hba_lb), "minimum # of hydrogen bond acceptors")
			("hba_ub", value<size_t>(&hba_ub)->default_value(default_hba_ub), "maximum # of hydrogen bond acceptors")
			("help", "this help information")
			("version", "version information")
			("config", value<path>(), "configuration file to load options from")
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

		// Validate idock example folder.
		if (!exists(idock_example_path))
		{
			cerr << "Option idock_example " << idock_example_path << " does not exist" << endl;
			return 1;
		}
		if (!is_directory(idock_example_path))
		{
			cerr << "Option idock_example " << idock_example_path << " is not a directory" << endl;
			return 1;
		}

		// Validate output folder.
		if (!exists(out_path) && !create_directories(out_path))
		{
			cerr << "Failed to create output folder " << out_path << endl;
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			cerr << "Option threads must be 1 or greater" << endl;
			return 1;
		}
	}
	catch (const std::exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	// Calculate the number of ligands (i.e. population size), which is the sum of number of elitists and children.
	const size_t num_ligands = num_elitists + num_children;

	// Initialize a pointer vector to dynamically hold and destroy generated ligands.
	boost::ptr_vector<ligand> ligands;
	ligands.resize(num_ligands);

	// Extract initial elite ligands from the idock example folder.
	cout << "Extracting " << num_elitists << " elite ligands from " << idock_example_path << endl;
	const path default_log_path = "log.csv";
	{
		const path idock_example_log_path = idock_example_path / default_log_path;
		boost::filesystem::ifstream ifs(idock_example_log_path);
		string line;

		// Skip the log file's header line, which is Ligand,nConfs,idock score (kcal/mol),RF-Score (pKd)
		getline(ifs, line);

		// Extract tuples of <Ligand, idock score, RF-Score> from the idock log file.
		vector<tuple<string, double, double>> records;
		for (const boost::char_separator<char> sep(","); getline(ifs, line);)
		{
			boost::tokenizer<boost::char_separator<char>> tokens(line, sep);
			auto tok_iter = tokens.begin();
			const auto stem = *tok_iter++;
			if (!stoul(*tok_iter++)) continue;
			records.emplace_back(stem, stod(*tok_iter++), stod(*tok_iter++));
			assert(tok_iter == tokens.end());
		}

		// Check if there are sufficient log records.
		const size_t num_records = records.size();
		if (num_records < num_elitists)
		{
			cerr << "Failed to construct an initial generation because the idock example log " << idock_example_log_path << " contains less than " << num_elitists << " ligands." << endl;
			return 1;
		}
		assert(num_elitists <= num_records);

		// Sort the idock log records by their idock score.
		sort(records.begin(), records.end(), [](const tuple<string, double, double>& record0, const tuple<string, double, double>& record1)
		{
			return get<1>(record0) < get<1>(record1);
		});

		// Find crossoverable elite ligands.
		for (size_t i = 0, j = 0; i < num_elitists && j < num_records; ++j)
		{
			// Create an elite ligand.
			const auto record = records[j];
			unique_ptr<ligand> elitist(new ligand(idock_example_path / (get<0>(record) + ".pdbqt")));

			// Skip the elite ligand if it is not crossoverable.
			if (!elitist->crossoverable()) continue;

			// Save the elite ligand in the dynamic vector.
			ligands.replace(i, elitist.release());
			ligands[i].id_score = get<1>(record);
			ligands[i].rf_score = get<2>(record);
			++i;
		}
		assert(ligands[num_elitists].p.empty());

		// Check if there are sufficient initial elite ligands.
		if (ligands[num_elitists - 1].p.empty())
		{
			cerr << "Failed to construct an initial generation because the idock example log " << idock_example_log_path << " contains less than " << num_elitists << " crossoverable ligands whose number of rotatable bonds is at least 1." << endl;
			return 1;
		}
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Seeding a random number generator with " << seed << endl;
	mt19937_64 rng(seed);

	// Initialize ligand filenames.
	vector<string> filenames(num_children);
	for (size_t i = 0; i < num_children; ++i)
	{
		filenames[i] = to_string(i) + ".pdbqt";
	}

	// Find the full path to idock executable.
	const path idock_path = path(search_path("idock")).make_preferred();
	cout << "Using idock executable at " << idock_path << endl;

	// Initialize arguments to idock.
	vector<string> idock_args(11);
	idock_args[0] = idock_path.string();
	idock_args[1] = "--config";
	idock_args[2] = "idock.conf";
	idock_args[3] = "--ligand";
	idock_args[5] = "--out";
	idock_args[7] = "--seed";
	idock_args[8] = to_string(seed);
	idock_args[9] = "--seed";
	idock_args[10]= num_threads;

	// Initialize an io service pool and create worker threads for later use.
	cout << "Creating an io service pool of " << num_threads << " worker thread" << (num_threads == 1 ? "" : "s") << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize the log file for writing statistics.
	boost::filesystem::ofstream log(out_path / default_log_path);
	log.setf(ios::fixed, ios::floatfield);
	log << "generation,ligand,parent 1,connector 1,parent 2,connector 2,idock score (kcal/mol),RF-Score (pKd),molecular mass (Da),rotatable bonds,hydrogen bond donors,hydrogen bond acceptors\n" << setprecision(2);

	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(2);
	for (size_t generation = 1; generation <= num_generations; ++generation)
	{
		cout << "Running GA generation " << generation << endl;

		// Initialize the paths to current generation folder and its two subfolders.
		const path generation_folder(out_path / to_string(generation));
		const path idock_folder(generation_folder / "idock"); // Used to save idock output.

		// Create a new folder and two subfolders for current generation.
		create_directory(generation_folder);
		create_directory(idock_folder);

		// Run crossover tasks.
		cout << "Executing " << num_children << " crossover operations in parallel" << endl;
		cnt.init(num_children);
		for (size_t i = 0; i < num_children; ++i)
		{
			const size_t s = rng();
			io.post([&, i, s]()
			{
				const size_t index = num_elitists + i;

				// Initialize a Mersenne Twister random number generator.
				mt19937_64 rng(s);
				uniform_int_distribution<size_t> uniform_elitist(0, num_elitists - 1);

				// Obtain references to two randomly selected parent ligands.
				ligand& l1 = ligands[uniform_elitist(rng)];
				ligand& l2 = ligands[uniform_elitist(rng)];

				// Obtain a random mutable atom from the two parent ligands respectively.
				const size_t g1 = uniform_int_distribution<size_t>(1, l1.num_rotatable_bonds)(rng);
				const size_t g2 = uniform_int_distribution<size_t>(1, l2.num_rotatable_bonds)(rng);

				ligands.replace(index, new ligand(generation_folder / filenames[i], l1, l2, g1, g2));
				ligands[index].save();
				cnt.increment();
			});
		}
		cnt.wait();

		// Invoke idock to dock generated ligands and save the docked conformations in the idock subfolder.
		cout << "Calling idock in the working directory of " << idock_example_path << endl;
		idock_args[4] = absolute(generation_folder).string();
		idock_args[6] = absolute(idock_folder).string();
		const auto exit_code = wait_for_exit(execute(start_in_dir(idock_example_path.string()), set_args(idock_args), inherit_env(), throw_on_error()));
		if (exit_code)
		{
			cerr << "idock exited with code " << exit_code << endl;
			return 1;
		}

		// Parse docked ligands to obtain predicted free energy and docked coordinates, and save the updated ligands.
		cout << "Refining ligands from docking results" << endl;
		for (size_t i = 0; i < num_children; ++i)
		{
			ligands[num_elitists + i].update(idock_folder / filenames[i]);
		}
		remove_all(idock_folder);

		// Sort ligands in ascending order of free energy.
		ligands.sort();

		// Write summaries to log.
		for (const auto& l : ligands)
		{
			log << generation
				<< ',' << l.p
				<< ',' << l.parent1
				<< ',' << l.connector1
				<< ',' << l.parent2
				<< ',' << l.connector2
				<< ',' << l.id_score
				<< ',' << l.rf_score
				<< ',' << l.mol_mass
				<< ',' << l.num_rotatable_bonds
				<< ',' << l.num_hb_donors
				<< ',' << l.num_hb_acceptors
				<< endl;
		}
	}
	io.wait();
}
