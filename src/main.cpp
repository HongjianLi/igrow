/*

   Copyright (c) 2011, The Chinese University of Hong Kong

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

	   http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

 */

/**
 * \mainpage igrow
 *
 * \section introduction Introduction
 * igrow is a multithreaded virtual screening tool for flexible ligand docking.
 *
 * \section features Features
 * igrow is inspired by AutoGrow. It uses either idock or AutoDock Vina as backend docking engine.
 * igrow supports more types of chemical synthesis such as halogen replacement and branch replacement in addition to hydrogen replacement.
 * igrow digests ligands and fragments in pdbqt format, saving the effort of frequently calling the prepare_ligand4 python script.
 * igrow invents its own thread pool in order to reuse threads and maintain a high CPU utilization throughout the entire synhsizing procedure. The thread pool parallelizes the creation of mutants and children in each generation.
 * igrow utilizes flyweight pattern for caching fragments and dynamic pointer vector for caching and sorting ligands.
 * igrow traces the sources of generated ligands and dumps the statistics in csv format so that users can easily get to know how the ligands are synthesized from the initial ligand and fragments.
 *
 * \section availability Availability
 * igrow is free and open source available at https://GitHub.com/HongjianLi/igrow under Apache License 2.0. Both x86 and x64 binaries for Linux and Windows are provided.
 *
 * \author Hongjian Li, The Chinese University of Hong Kong.
 * \date 22 December 2011
 *
 * Copyright (C) 2011 The Chinese University of Hong Kong.
 */

#include <boost/thread/thread.hpp>
#include <boost/program_options.hpp>
#include <boost/random.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/process/context.hpp>
#include <boost/process/operations.hpp>
#include <boost/process/child.hpp>
#include "seed.hpp"
#include "fstream.hpp"
#include "tee.hpp"
#include "ligand.hpp"

// Choose the appropriate Mersenne Twister engine for random number generation on 32-bit or 64-bit platform.
#if defined(__x86_64) || defined(__x86_64__) || defined(__amd64) || defined(__amd64__) || defined(_M_X64) || defined(_M_AMD64)
typedef boost::random::mt19937_64 mt19937eng;
#else
typedef boost::random::mt19937 mt19937eng;
#endif

int
main(int argc, char* argv[])
{
	std::cout << "igrow 1.0\n";

	using namespace igrow;
	path fragment_folder_path, initial_ligand_path, docking_program_path, docking_config_path, output_folder_path, log_path, csv_path;
	size_t num_threads, seed, num_generations, num_elitists, num_mutants, num_children, max_atoms, max_hb_donors, max_hb_acceptors;
	fl max_mw, max_logp;
	bool idock;

	// Process program options.
	{
		using namespace boost::program_options;

		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const path default_log_path = "log.txt";
		const path default_csv_path = "log.csv";
		const unsigned int concurrency = boost::thread::hardware_concurrency();
		const unsigned int default_num_threads = concurrency ? concurrency : 1;
		const size_t default_seed = random_seed();
		const size_t default_num_generations = 8;
		const size_t default_num_mutants = 20;
		const size_t default_num_children = 20;
		const size_t default_num_elitists = 10;
		const size_t default_max_atoms = 80;
		const size_t default_max_hb_donors = 5;
		const size_t default_max_hb_acceptors = 10;
		const fl default_max_mw = 500;
		const fl default_max_logp = 5;

		options_description input_options("input (required)");
		input_options.add_options()
			("fragment_folder", value<path > (&fragment_folder_path)->required(), "path to folder of fragments in PDB format")
			("initial_ligand", value<path > (&initial_ligand_path)->required(), "path to initial ligand in PDB format")
			("docking_program", value<path > (&docking_program_path)->required(), "path to igrow or vina executable")
			("docking_config", value<path > (&docking_config_path)->required(), "path to igrow or vina configuration file")
			;

		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path > (&output_folder_path)->default_value(default_output_folder_path), "folder of output results")
			("log", value<path > (&log_path)->default_value(default_log_path), "log file in plain text")
			("csv", value<path > (&csv_path)->default_value(default_csv_path), "summary file in csv format")
			;

		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("threads", value<size_t > (&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t > (&seed)->default_value(default_seed), "explicit non-negative random seed")
			("generations", value<size_t > (&num_generations)->default_value(default_num_generations), "number of GA generations")
			("elitists", value<size_t > (&num_elitists)->default_value(default_num_elitists), "number of elite ligands to carry over")
			("mutants", value<size_t > (&num_mutants)->default_value(default_num_mutants), "number of child ligands generated by mutation")
			("children", value<size_t > (&num_children)->default_value(default_num_children), "number of child ligands generated by crossover")
			("max_atoms", value<size_t > (&max_atoms)->default_value(default_max_atoms), "maximum number of atoms")
			("max_hb_donors", value<size_t > (&max_hb_donors)->default_value(default_max_hb_donors), "maximum number of hydrogen bond donors")
			("max_hb_acceptors", value<size_t > (&max_hb_acceptors)->default_value(default_max_hb_acceptors), "maximum number of hydrogen bond acceptors")
			("max_mw", value<fl > (&max_mw)->default_value(default_max_mw), "maximum molecular weight")
			("max_logp", value<fl > (&max_logp)->default_value(default_max_logp), "maximum logP")
			("config", value<path > (), "options can be loaded from a configuration file")
			;

		options_description all_options;
		all_options.add(input_options).add(output_options).add(miscellaneous_options);

		// If no command line argument is supplied, simply print the usage and exit.
		if (argc == 1)
		{
			std::cout << all_options;
			return 0;
		}

		// Parse command line arguments.
		try
		{
			variables_map vm;
			store(parse_command_line(argc, argv, all_options), vm);
			variable_value config_value = vm["config"];
			if (!config_value.empty()) // If a configuration file is presented, parse it.
			{
				ifstream config_file(config_value.as<path > ());
				store(parse_config_file(config_file, all_options), vm);
			}
			vm.notify(); // Notify the user if there are any parsing errors.
		}
		catch (const std::exception& e)
		{
			std::cerr << e.what() << '\n';
			return 1;
		}

		// Validate fragment folder.
		if (!exists(fragment_folder_path))
		{
			std::cerr << "Fragment folder " << fragment_folder_path << " does not exist\n";
			return 1;
		}
		if (!is_directory(fragment_folder_path))
		{
			std::cerr << "Fragment folder " << fragment_folder_path << " is not a directory\n";
			return 1;
		}

		// Validate initial ligand.
		if (!exists(initial_ligand_path))
		{
			std::cerr << "Initial ligand " << initial_ligand_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(initial_ligand_path))
		{
			std::cerr << "Initial ligand " << initial_ligand_path << " is not a regular file\n";
			return 1;
		}

		// Validate docking program.
		if (!exists(docking_program_path))
		{
			std::cerr << "Docking program " << docking_program_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(docking_program_path))
		{
			std::cerr << "Docking program " << docking_program_path << " is not a regular file\n";
			return 1;
		}

		// Determine if the docking program is igrow or vina.
		const string docking_program = docking_program_path.stem().string();
		if (docking_program == "idock")
		{
			idock = true;
		}
		else if (docking_program == "vina")
		{
			idock = false;
		}
		else
		{
			std::cerr << "Docking program must be either igrow or vina\n";
			return 1;
		}

		// Validate docking configuration file.
		if (!exists(docking_config_path))
		{
			std::cerr << "Docking configuration file " << docking_config_path << " does not exist\n";
			return 1;
		}
		if (!is_regular_file(docking_config_path))
		{
			std::cerr << "Docking configuration file " << docking_config_path << " is not a regular file\n";
			return 1;
		}

		// Validate output_folder.
		remove_all(output_folder_path);
		if (!create_directories(output_folder_path))
		{
			std::cerr << "Failed to create output folder " << output_folder_path << '\n';
			return 1;
		}

		// Validate miscellaneous options.
		if (num_threads < 1)
		{
			std::cerr << "Option threads must be 1 or greater\n";
			return 1;
		}
		if (num_generations < 1)
		{
			std::cerr << "Option generations must be 1 or greater\n";
			return 1;
		}
		if (num_mutants < 1)
		{
			std::cerr << "Option mutants must be 1 or greater\n";
			return 1;
		}
		if (num_children < 0)
		{
			std::cerr << "Option children must be 0 or greater\n";
			return 1;
		}
		if (num_elitists < 0)
		{
			std::cerr << "Option carryovers must be 0 or greater\n";
			return 1;
		}
		if (max_atoms <= 0)
		{
			std::cerr << "Option max_atoms must be 1 or greater\n";
			return 1;
		}
	}

	try
	{
		// Initialize the log.
		std::cout << "Logging to " << log_path.string() << '\n';
		igrow::tee log(log_path);

		// Initialize a Mersenne Twister random number generator.
		log << "Using random seed " << seed << '\n';
		mt19937eng eng(seed);

		// Scan the fragment folder to obtain a list of fragments.
		log << "Scanning fragment folder " << fragment_folder_path.string() << '\n';
		vector<path> fragment_paths;
		fragment_paths.reserve(1000); // A fragment library typically consists of <= 1000 fragments.
		{
			using namespace boost::filesystem;
			const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
			for (directory_iterator dir_iter(fragment_folder_path); dir_iter != end_dir_iter; ++dir_iter)
			{
				// Skip non-regular files such as folders.
				if (!is_regular_file(dir_iter->status())) continue;

				// Save the fragment path.
				fragment_paths.push_back(dir_iter->path());
			}
		}
		const size_t num_fragments = fragment_paths.size();
		log << num_fragments << " fragments found\n";

		// Initialize random number generators.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t> > uniform_fragment_gen(eng, uniform_int_distribution<size_t>(0, num_fragments - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t> > uniform_elitist_gen(eng, uniform_int_distribution<size_t>(0, num_elitists - 1));

		// Initialize process context.
		using namespace boost::process;
		context ctx;

		// Initialize arguments to igrow or vina.
		using namespace boost;
		vector<string> docking_args(8);
		docking_args[0] = "--config";
		docking_args[1] = docking_config_path.string();
		docking_args[6] = "--seed";
		docking_args[7] = lexical_cast<string>(seed);
		if (idock)
		{
			docking_args[2] = "--ligand_folder";
			docking_args[4] = "--output_folder";			
		}
		else
		{
			docking_args[2] = "--ligand";
			docking_args[4] = "--out";
		}

		// The number of ligands is equal to the number of elitists plus mutants plus children.
		const size_t num_ligands = num_elitists + num_mutants + num_children;
		ptr_vector<ligand> ligands(num_ligands);

		// Initialize csv file for dumping statistics.
		ofstream csv(csv_path);
		csv << "generation,ligand,free energy in kcal/mol,no. of heavy atoms,no. of hydrogen bond donors,no. of hydrogen bond acceptors,molecular weight,logp\n";

		for (size_t current_generation = 1; current_generation <= num_generations; ++current_generation)
		{
			log << "Running generation " << current_generation << '\n';

			// Initialize the paths to current generation folder and its two subfolders.
			const path current_generation_folder_path(output_folder_path / path(lexical_cast<string > (current_generation)));
			const path current_pdbqt_folder_path (current_generation_folder_path / path("pdbqt"));
			const path current_output_folder_path(current_generation_folder_path / path("output"));

			// Create a new folder and four subfolders for current generation.
			create_directory(current_generation_folder_path);
			create_directory(current_pdbqt_folder_path);
			create_directory(current_output_folder_path);

			// Generate ligands and save them into the pdbqt subfolder.
			if (current_generation == 1)
			{
				// Parse and dump the initial ligand.
				ligands.push_back(new ligand(initial_ligand_path));
				ligands.front().save(current_pdbqt_folder_path / "1.pdbqt");

				// Create mutants in parallel.
				for (size_t i = 1; i != num_ligands; ++i)
				{
					/*const*/ ligand lig(initial_ligand_path);					
				    lig.mutate(ligand_flyweight(fragment_paths[uniform_fragment_gen()]));
					// Check ligand validity.
					lig.save(current_pdbqt_folder_path / (lexical_cast<string > (i + 1) + ".pdbqt"));
				}
			}
			else
			{
				// TODO: abstract into tasks for parallel execution.
				for (size_t i = 0; i < num_mutants; ++i)
				{
					ligand lig = ligands[uniform_elitist_gen()];
				    lig.mutate(ligand_flyweight(fragment_paths[uniform_fragment_gen()]));
					ligands.replace(num_elitists + i, &lig);
				}

				// TODO: abstract into tasks for parallel execution.
				for (size_t i = 0; i < num_children; ++i)
				{
					//const path previous_generation_folder_path(output_folder_path / path(lexical_cast<string > (current_generation - 1)));
					//	    Interaction interact;
					//	    Ligand child = interact.mate(m1, m2);
					//	    if (child.valid());
					//ligands.replace(num_elitists + num_mutants + i, &lig);
				}
			}

			// Call either igrow or vina to dock ligands to predict their free energy and then parse the docking log.
			// TODO: Parse resultant pdbqt file instead.
			ctx.work_dir = current_generation_folder_path.string();
			if (idock)
			{
				// Invoke idock.
				log << "Calling idock to dock " << num_ligands << " ligands\n";
				docking_args[3] = "pdbqt";
				docking_args[5] = "output";
				create_child(docking_program_path.string(), docking_args, ctx).wait();

				// Parse idock log.
				//ifstream log(current_log_folder_path / path("log"));
				//log.seekg(103);
				//string line;
				//while (getline(log, line))
				//{
				//	if (line[0] == ' ') break; //   index |       ligand |   progress | conf | top 5 conf free energy in kcal/mol
				//}
				//while (getline(log, line))
				//{
				//	const size_t ligand_id = right_cast<size_t>(line, 11, 22);
				//	const fl free_energy = right_cast<fl>(line, 46, 51);
				//	ligands[ligand_id].free_energy = free_energy;
				//}
				//log.close();
			}
			else
			{
				// Invoke vina.
				log << "Calling vina to dock " << num_ligands << " ligands\n";
				for (size_t i = 1; i <= num_ligands; ++i)
				{
					docking_args[3] = (path("pdbqt")  / path(lexical_cast<string > (i) + ".pdbqt")).string();
					docking_args[5] = (path("output") / path(lexical_cast<string > (i) + ".pdbqt")).string();
					create_child(docking_program_path.string(), docking_args, ctx).wait();
				}

				// Parse vina log.
				//for (size_t i = 1; i <= num_ligands; ++i)
				//{
				//	ifstream log(current_log_folder_path / path(lexical_cast<string > (i) + ".log"));
				//	log.seekg(1177);
				//	string line;
				//	while (getline(log, line))
				//	{
				//		if (line[3] == '1')
				//		{
				//			const size_t start = line.find_first_not_of(' ', 8);
				//			const fl free_energy = lexical_cast<fl > (line.substr(start, 17 - start));
				//			ligands[i - 1].free_energy = free_energy;
				//		}
				//	}
				//	log.close();
				//}
			}

			// Sort ligands in ascending order of predicted free energy
			ligands.sort();

			// Write csv
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
}
