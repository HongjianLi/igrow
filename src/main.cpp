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
using namespace boost;
using namespace boost::filesystem;
using namespace boost::process;
using namespace boost::process::initializers;

int main(int argc, char* argv[])
{
	// Initialize the default path to log files. They will be reused when calling idock.
	const path default_log_path = "log.csv";

	path initial_generation_csv_path, initial_generation_folder_path, fragment_folder_path, idock_config_path, output_folder_path, log_path;
	size_t num_threads, seed, num_elitists, num_crossovers, max_failures, max_rotatable_bonds, max_hb_donors, max_hb_acceptors;
	double max_mw;

	// Process program options.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_num_crossovers = 20;
		const size_t default_num_elitists = 10;
		const size_t default_max_failures = 1000;
		const size_t default_max_rotatable_bonds = 30;
		const size_t default_max_hb_donors = 5;
		const size_t default_max_hb_acceptors = 10;
		const double default_max_mw = 500;

		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("initial_generation_csv", value<path>(&initial_generation_csv_path)->required(), "path to initial generation csv")
			("initial_generation_folder", value<path>(&initial_generation_folder_path)->required(), "path to initial generation folder")
			("idock_config", value<path>(&idock_config_path)->required(), "path to idock configuration file")
			;
		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output results")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file in csv format")
			;
		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("elitists", value<size_t>(&num_elitists)->default_value(default_num_elitists), "number of elite ligands to carry over")
			("crossovers", value<size_t>(&num_crossovers)->default_value(default_num_crossovers), "number of child ligands created by crossover")
			("max_failures", value<size_t>(&max_failures)->default_value(default_max_failures), "maximum number of operational failures to tolerate")
			("max_rotatable_bonds", value<size_t>(&max_rotatable_bonds)->default_value(default_max_rotatable_bonds), "maximum number of rotatable bonds")
			("max_hb_donors", value<size_t>(&max_hb_donors)->default_value(default_max_hb_donors), "maximum number of hydrogen bond donors")
			("max_hb_acceptors", value<size_t>(&max_hb_acceptors)->default_value(default_max_hb_acceptors), "maximum number of hydrogen bond acceptors")
			("max_mw", value<double>(&max_mw)->default_value(default_max_mw), "maximum molecular weight")
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

		// Validate initial generation csv.
		if (!exists(initial_generation_csv_path))
		{
			cerr << "Initial generation csv " << initial_generation_csv_path << " does not exist" << endl;
			return 1;
		}
		if (!is_regular_file(initial_generation_csv_path))
		{
			cerr << "Initial generation csv " << initial_generation_csv_path << " is not a regular file" << endl;
			return 1;
		}

		// Validate initial generation folder.
		if (!exists(initial_generation_folder_path))
		{
			cerr << "Initial generation folder " << initial_generation_folder_path << " does not exist" << endl;
			return 1;
		}
		if (!is_directory(initial_generation_folder_path))
		{
			cerr << "Initial generation folder " << initial_generation_folder_path << " is not a directory" << endl;
			return 1;
		}

		// Validate fragment folder.
		if (!exists(fragment_folder_path))
		{
			cerr << "Fragment folder " << fragment_folder_path << " does not exist" << endl;
			return 1;
		}
		if (!is_directory(fragment_folder_path))
		{
			cerr << "Fragment folder " << fragment_folder_path << " is not a directory" << endl;
			return 1;
		}

		// Validate idock configuration file.
		if (!exists(idock_config_path))
		{
			cerr << "idock configuration file " << idock_config_path << " does not exist" << endl;
			return 1;
		}
		if (!is_regular_file(idock_config_path))
		{
			cerr << "idock configuration file " << idock_config_path << " is not a regular file" << endl;
			return 1;
		}

		// Validate output folder.
		remove_all(output_folder_path);
		if (!create_directories(output_folder_path))
		{
			cerr << "Failed to create output folder " << output_folder_path << endl;
			return 1;
		}

		// Validate log_path.
		if (is_directory(log_path))
		{
			cerr << "log path " << log_path << " is a directory" << endl;
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			cerr << "Option threads must be 1 or greater" << endl;
			return 1;
		}
		if (max_mw <= 0)
		{
			cerr << "Option max_mw must be positive" << endl;
			return 1;
		}
	}
	catch (const std::exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	// The number of ligands (i.e. population size) is equal to the number of elitists plus mutants plus children.
	const size_t num_children = num_crossovers;
	const size_t num_ligands = num_elitists + num_children;
	const double num_elitists_inv = static_cast<double>(1) / num_elitists;

	// Initialize a pointer vector to dynamically hold and destroy generated ligands.
	ptr_vector<ligand> ligands;
	ligands.resize(num_ligands);

	// Parse the initial generation csv to get initial elite ligands.
	{
		boost::filesystem::ifstream ifs(initial_generation_csv_path);
		string line;
		line.reserve(80);
		getline(ifs, line); // Ligand,pKd1,pKd2,pKd3,pKd4,pKd5,pKd6,pKd7,pKd8,pKd9
		for (size_t i = 0; i < num_elitists; ++i)
		{
			// Check if there are sufficient initial elite ligands.
			if (!getline(ifs, line))
			{
				cerr << "Failed to construct initial generation because the initial generation csv " << initial_generation_csv_path << " contains less than " << num_elitists << " ligands." << endl;
				return 1;
			}

			// Parse the elite ligand.
			const size_t comma1 = line.find(',', 1);
			ligands.replace(i, new ligand(initial_generation_folder_path / (line.substr(0, comma1) + ".pdbqt")));

			// Parse the free energy and ligand efficiency.
			const size_t comma2 = line.find(',', comma1 + 2);
			ligands[i].fe = stod(line.substr(comma1 + 1, comma2 - comma1 - 1));
		}
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << endl;
	mt19937_64 eng(seed);

	// Initialize a ligand validator.
	const validator v(max_rotatable_bonds, max_hb_donors, max_hb_acceptors, max_mw);

	// Initialize the number of failures. The program will stop if num_failures reaches max_failures.
	atomic<size_t> num_failures(0);

	// Initialize ligand filenames.
	vector<string> ligand_filenames;
	ligand_filenames.reserve(num_ligands);
	for (size_t i = 1; i <= num_ligands; ++i)
	{
		ligand_filenames.push_back(to_string(i) + ".pdbqt");
	}

	// Find the full path to idock executable.
	const path idock_path = path(search_path("idock")).make_preferred();
	cout << "Using idock executable at " << idock_path << endl;

	// Initialize arguments to idock.
	vector<string> idock_args(10);
	idock_args[0] = "--input_folder";
	idock_args[2] = "--output_folder";
	idock_args[4] = "--log";
	idock_args[6] = "--seed";
	idock_args[7] = to_string(seed);
	idock_args[8] = "--config";
	idock_args[9] = idock_config_path.string();

	// Initialize an io service pool and create worker threads for later use.
	cout << "Creating an io service pool of " << num_threads << " worker thread" << (num_threads == 1 ? "" : "s") << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize log file for dumping statistics.
	boost::filesystem::ofstream log(log_path);
	log << "generation,ligand,parent 1,connector 1,parent 2,connector 2,free energy (kcal/mol),rotatable bonds,hydrogen bond donors,hydrogen bond acceptors,molecular weight (g/mol)\n";

	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(3);
	for (size_t generation = 1; true; ++generation)
	{
		cout << "Running generation " << generation << endl;

		// Initialize the paths to current generation folder and its two subfolders.
		const path generation_folder(output_folder_path / to_string(generation));
		const path  input_folder(generation_folder /  "input");
		const path output_folder(generation_folder / "output");

		// Create a new folder and two subfolders for current generation.
		create_directory(generation_folder);
		create_directory( input_folder);
		create_directory(output_folder);

		// Create addition, subtraction and crossover tasks.
		cnt.init(num_children);
		for (size_t i = 0; i < num_children; ++i)
		{
			const size_t s = eng();
			io.post([&, i, s]()
			{
				const size_t index = num_elitists + i;

				// Initialize a Mersenne Twister random number generator.
				mt19937_64 eng(s);
				uniform_int_distribution<size_t> uniform_elitist(0, num_elitists - 1);

				// Create a child ligand by crossover.
				do
				{
					// Obtain constant references to the two parent ligands.
					ligand& l1 = ligands[uniform_elitist(eng)];
					ligand& l2 = ligands[uniform_elitist(eng)];
					while (!(l1.crossover_feasible() && l2.crossover_feasible()))
					{
						l1 = ligands[uniform_elitist(eng)];
						l2 = ligands[uniform_elitist(eng)];
					}

					// Obtain a random mutable atom from the two parent ligands respectively.
					const size_t g1 = uniform_int_distribution<size_t>(1, l1.num_rotatable_bonds)(eng);
					const size_t g2 = uniform_int_distribution<size_t>(1, l2.num_rotatable_bonds)(eng);

					ligands.replace(index, new ligand(input_folder / ligand_filenames[i], l1, l2, g1, g2));
					if (v(ligands[index]))
					{
						// Save the newly created child ligand.
						ligands[index].save();
						break;
					}
				} while (++num_failures < max_failures);
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
		idock_args[1] =  input_folder.string();
		idock_args[3] = output_folder.string();
		idock_args[5] = (generation_folder / default_log_path).string();
		const auto exit_code = wait_for_exit(execute(run_exe(idock_path), set_args(idock_args), throw_on_error()));
		if (exit_code)
		{
			cerr << "idock exited with code " << exit_code << endl;
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
		for (const auto& l : ligands)
		{
			log << generation
				<< ',' << l.p
				<< ',' << l.parent1
				<< ',' << l.connector1
				<< ',' << l.parent2
				<< ',' << l.connector2
				<< ',' << l.fe
				<< ',' << l.num_rotatable_bonds
				<< ',' << l.num_hb_donors
				<< ',' << l.num_hb_acceptors
				<< ',' << l.mw
				<< endl;
		}

		// Calculate average statistics of elite ligands.
		double avg_mw = 0, avg_fe = 0, avg_le = 0, avg_rotatable_bonds = 0, avg_hb_donors = 0, avg_hb_acceptors = 0;
		for (size_t i = 0; i < num_elitists; ++i)
		{
			const ligand& l = ligands[i];
			avg_mw += l.mw;
			avg_fe += l.fe;
			avg_rotatable_bonds += l.num_rotatable_bonds;
			avg_hb_donors += l.num_hb_donors;
			avg_hb_acceptors += l.num_hb_acceptors;
		}
		avg_mw *= num_elitists_inv;
		avg_fe *= num_elitists_inv;
		avg_rotatable_bonds *= num_elitists_inv;
		avg_hb_donors *= num_elitists_inv;
		avg_hb_acceptors *= num_elitists_inv;
		cout << "Failures |  Avg FE | Avg MWT | Avg NRB | Avg HBD | Avg HBA\n"
		    << setw(8) << num_failures << "   "
			<< setw(7) << avg_fe << "   "
			<< setw(7) << avg_mw << "   "
			<< setw(7) << avg_rotatable_bonds << "   "
			<< setw(7) << avg_hb_donors << "   "
			<< setw(7) << avg_hb_acceptors << endl;
	}
}
