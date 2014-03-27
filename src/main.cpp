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

	path idock_example_folder_path, output_folder_path, log_path;
	size_t num_threads, seed, num_elitists, num_crossovers, num_generations;

	// Process program options.
	try
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const size_t default_seed = chrono::system_clock::now().time_since_epoch().count();
		const size_t default_num_threads = thread::hardware_concurrency();
		const size_t default_num_crossovers = 20;
		const size_t default_num_elitists = 10;
		const size_t default_num_generations = 8;

		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("idock_example_folder", value<path>(&idock_example_folder_path)->required(), "path to an idock example folder")
			;
		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output results")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file in csv format")
			;
		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("elitists", value<size_t>(&num_elitists)->default_value(default_num_elitists), "number of elite ligands to carry over")
			("crossovers", value<size_t>(&num_crossovers)->default_value(default_num_crossovers), "number of child ligands created by crossover")
			("generations", value<size_t>(&num_generations)->default_value(default_num_generations), "number of generations")
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

		// Validate idock example folder.
		if (!exists(idock_example_folder_path))
		{
			cerr << "idock example folder " << idock_example_folder_path << " does not exist" << endl;
			return 1;
		}
		if (!is_directory(idock_example_folder_path))
		{
			cerr << "idock example folder " << idock_example_folder_path << " is not a directory" << endl;
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
	}
	catch (const std::exception& e)
	{
		cerr << e.what() << endl;
		return 1;
	}

	// The number of ligands (i.e. population size) is equal to the number of elitists plus children.
	const size_t num_children = num_crossovers;
	const size_t num_ligands = num_elitists + num_children;

	// Initialize a pointer vector to dynamically hold and destroy generated ligands.
	ptr_vector<ligand> ligands;
	ligands.resize(num_ligands);

	// Parse the idock example folder to get initial elite ligands.
	{
		const path idock_example_log_path = idock_example_folder_path / default_log_path;
		const path idock_example_output_path = idock_example_folder_path / "output";
		boost::filesystem::ifstream ifs(idock_example_log_path);
		string line;
		getline(ifs, line); // Ligand,Energy1,Energy2,Energy3,Energy4,Energy5,Energy6,Energy7,Energy8,Energy9
		size_t i = 0;
		while (getline(ifs, line) && i < num_elitists)
		{
			// Parse the elite ligand.
			const size_t comma1 = line.find(',', 1);
			unique_ptr<ligand> elitist(new ligand(idock_example_output_path / (line.substr(0, comma1) + ".pdbqt")));

			// Skip elite ligands that are not feasible of crossover.
			if (!elitist->crossover_feasible()) continue;
			ligands.replace(i, elitist.release());

			// Parse the free energy.
			const size_t comma2 = line.find(',', comma1 + 2);
			ligands[i].fe = stod(line.substr(comma1 + 1, comma2 - comma1 - 1));

			// Move to the next position.
			++i;
		}

		// Check if there are sufficient initial elite ligands.
		if (i < num_elitists)
		{
			cerr << "Failed to construct initial generation because the idock example log " << idock_example_log_path << " contains less than " << num_elitists << " crossoverable ligands." << endl;
			return 1;
		}
	}

	// Initialize a Mersenne Twister random number generator.
	cout << "Using random seed " << seed << endl;
	mt19937_64 rng(seed);

	// Initialize ligand filenames.
	vector<string> ligand_filenames(num_children);
	for (size_t i = 0; i < num_children; ++i)
	{
		ligand_filenames[i] = to_string(i) + ".pdbqt";
	}

	// Find the full path to idock executable.
	const path idock_path = path(search_path("idock")).make_preferred();
	cout << "Using idock executable at " << idock_path << endl;

	// Initialize arguments to idock.
	vector<string> idock_args(11);
	idock_args[0] = idock_path.string();
	idock_args[1] = "--config";
	idock_args[2] = "idock.conf";
	idock_args[3] = "--input_folder";
	idock_args[5] = "--output_folder";
	idock_args[7] = "--log";
	idock_args[9] = "--seed";
	idock_args[10]= to_string(seed);

	// Initialize an io service pool and create worker threads for later use.
	cout << "Creating an io service pool of " << num_threads << " worker thread" << (num_threads == 1 ? "" : "s") << endl;
	io_service_pool io(num_threads);
	safe_counter<size_t> cnt;

	// Initialize log file for writing statistics.
	boost::filesystem::ofstream log(log_path);
	log << "generation,ligand,parent 1,connector 1,parent 2,connector 2,free energy (kcal/mol),rotatable bonds,hydrogen bond donors,hydrogen bond acceptors,molecular mass (daltons)\n";

	cout.setf(ios::fixed, ios::floatfield);
	cout << setprecision(3);
	for (size_t generation = 1; generation <= num_generations; ++generation)
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

		// Create crossover tasks.
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

				// Obtain constant references to the two parent ligands.
				ligand& l1 = ligands[uniform_elitist(rng)];
				ligand& l2 = ligands[uniform_elitist(rng)];

				// Obtain a random mutable atom from the two parent ligands respectively.
				const size_t g1 = uniform_int_distribution<size_t>(1, l1.num_rotatable_bonds)(rng);
				const size_t g2 = uniform_int_distribution<size_t>(1, l2.num_rotatable_bonds)(rng);

				ligands.replace(index, new ligand(input_folder / ligand_filenames[i], l1, l2, g1, g2));
				ligands[index].save();
				cnt.increment();
			});
		}
		cnt.wait();

		// Invoke idock.
		idock_args[4] = absolute( input_folder).string();
		idock_args[6] = absolute(output_folder).string();
		idock_args[8] = absolute(generation_folder / default_log_path).string();
		const auto exit_code = wait_for_exit(execute(start_in_dir(idock_example_folder_path.string()), set_args(idock_args), inherit_env(), throw_on_error()));
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
				<< ',' << l.fe
				<< ',' << l.num_rotatable_bonds
				<< ',' << l.num_hb_donors
				<< ',' << l.num_hb_acceptors
				<< ',' << l.ma
				<< endl;
		}
	}
}
