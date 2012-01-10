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
 * \date 10 January 2012
 *
 * Copyright (C) 2011-2012 The Chinese University of Hong Kong.
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

int main(int argc, char* argv[])
{
	std::cout << "igrow 1.0\n";

	using namespace igrow;
	path fragment_folder_path, initial_ligand_path, docking_program_path, docking_config_path, output_folder_path, log_path, csv_path;
	size_t num_threads, seed, num_generations, num_elitists, num_mutants, num_crossovers, max_failures, max_rotatable_bonds, max_atoms, max_heavy_atoms, max_hb_donors, max_hb_acceptors;
	fl max_mw, max_logp, min_logp;
	bool idock; ///< True if the docking program is idock, false if vina.

	// Initialize the default path to log files. They will be reused when calling idock.
	const path default_log_path = "log.txt";
	const path default_csv_path = "log.csv";

	// Process program options.
	{
		// Initialize the default values of optional arguments.
		const path default_output_folder_path = "output";
		const string default_docking_program = "idock";
		const unsigned int concurrency = boost::thread::hardware_concurrency();
		const size_t default_num_threads = concurrency ? concurrency : 1;
		const size_t default_seed = random_seed();
		const size_t default_num_generations = 8;
		const size_t default_num_mutants = 20;
		const size_t default_num_crossovers = 20;
		const size_t default_num_elitists = 10;
		const size_t default_max_failures = 1000;
		const size_t default_max_rotatable_bonds = 30;
		const size_t default_max_atoms = 100;
		const size_t default_max_heavy_atoms = 80;
		const size_t default_max_hb_donors = 5;
		const size_t default_max_hb_acceptors = 10;
		const fl default_max_mw = 500;
		const fl default_max_logp = 5;
		const fl default_min_logp = -5;

		using namespace boost::program_options;
		options_description input_options("input (required)");
		input_options.add_options()
			("fragment_folder", value<path>(&fragment_folder_path)->required(), "path to folder of fragments in pdbqt format")
			("initial_ligand", value<path>(&initial_ligand_path)->required(), "path to initial ligand in pdbqt format")
			("docking_config", value<path>(&docking_config_path)->required(), "path to docking configuration file")
			;

		options_description output_options("output (optional)");
		output_options.add_options()
			("output_folder", value<path>(&output_folder_path)->default_value(default_output_folder_path), "folder of output results")
			("log", value<path>(&log_path)->default_value(default_log_path), "log file in plain text")
			("csv", value<path>(&csv_path)->default_value(default_csv_path), "summary file in csv format")
			;

		string docking_program;
		options_description miscellaneous_options("options (optional)");
		miscellaneous_options.add_options()
			("docking_program", value<string>(&docking_program)->default_value(default_docking_program), "either idock or vina")
			("threads", value<size_t>(&num_threads)->default_value(default_num_threads), "number of worker threads to use")
			("seed", value<size_t>(&seed)->default_value(default_seed), "explicit non-negative random seed")
			("generations", value<size_t>(&num_generations)->default_value(default_num_generations), "number of GA generations")
			("elitists", value<size_t>(&num_elitists)->default_value(default_num_elitists), "number of elite ligands to carry over")
			("mutants", value<size_t>(&num_mutants)->default_value(default_num_mutants), "number of child ligands created by mutation")
			("crossovers", value<size_t>(&num_crossovers)->default_value(default_num_crossovers), "number of child ligands created by crossover")
			("max_failures", value<size_t>(&max_failures)->default_value(default_max_failures), "maximum number of operational failures to tolerate")
			("max_rotatable_bonds", value<size_t>(&max_rotatable_bonds)->default_value(default_max_rotatable_bonds), "maximum number of rotatable bonds")
			("max_atoms", value<size_t>(&max_atoms)->default_value(default_max_atoms), "maximum number of atoms")
			("max_heavy_atoms", value<size_t>(&max_heavy_atoms)->default_value(default_max_heavy_atoms), "maximum number of heavy atoms")
			("max_hb_donors", value<size_t>(&max_hb_donors)->default_value(default_max_hb_donors), "maximum number of hydrogen bond donors")
			("max_hb_acceptors", value<size_t>(&max_hb_acceptors)->default_value(default_max_hb_acceptors), "maximum number of hydrogen bond acceptors")
			("max_mw", value<fl>(&max_mw)->default_value(default_max_mw), "maximum molecular weight")
			("max_logp", value<fl>(&max_logp)->default_value(default_max_logp), "maximum logP")
			("min_logp", value<fl>(&min_logp)->default_value(default_min_logp), "minimum logP")
			("config", value<path>(), "options can be loaded from a configuration file")
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
				ifstream config_file(config_value.as<path>());
				store(parse_config_file(config_file, all_options), vm);
			}
			vm.notify(); // Notify the user if there are any parsing errors.

			// Determine if the docking program is idock or vina.
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
				std::cerr << "Docking program must be either idock or vina\n";
				return 1;
			}

			// Find the full path to the docking executable.
			docking_program_path = boost::process::find_executable_in_path(docking_program);
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

		// Validate output folder.
		remove_all(output_folder_path);
		if (!create_directories(output_folder_path))
		{
			std::cerr << "Failed to create output folder " << output_folder_path << '\n';
			return 1;
		}

		// Validate miscellaneous options.
		if (!num_threads)
		{
			std::cerr << "Option threads must be 1 or greater\n";
			return 1;
		}
		if (!num_generations)
		{
			std::cerr << "Option generations must be 1 or greater\n";
			return 1;
		}
		if (!num_mutants)
		{
			std::cerr << "Option mutants must be 1 or greater\n";
			return 1;
		}
		if (max_mw <= 0)
		{
			std::cerr << "Option max_mw must be positive\n";
			return 1;
		}
		if (min_logp > max_logp)
		{
			std::cerr << "Option max_mw must be larger than or equal to option min_mw\n";
			return 1;
		}
	}

	try
	{
		// Initialize the log.
		std::cout << "Logging to " << log_path.string() << '\n';
		igrow::tee log(log_path);

		// Scan the fragment folder to obtain a list of fragments.
		log << "Scanning fragment folder " << fragment_folder_path.string() << '\n';
		vector<path> fragments;
		fragments.reserve(1000); // A fragment folder typically consists of <= 1000 fragments.
		{
			using namespace boost::filesystem;
			const directory_iterator end_dir_iter; // A default constructed directory_iterator acts as the end iterator.
			for (directory_iterator dir_iter(fragment_folder_path); dir_iter != end_dir_iter; ++dir_iter)
			{
				// Skip non-regular files such as folders.
				if (!is_regular_file(dir_iter->status())) continue;

				// Save the fragment path.
				fragments.push_back(dir_iter->path());
			}
		}
		const size_t num_fragments = fragments.size();
		log << num_fragments << " fragments found\n";

		// Initialize a Mersenne Twister random number generator.
		log << "Using random seed " << seed << '\n';
		mt19937eng eng(seed);

		// Initialize random number generators for obtaining a random fragment and a random elitist.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_fragment_gen(eng, uniform_int_distribution<size_t>(0, num_fragments - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_elitist_gen (eng, uniform_int_distribution<size_t>(0, num_elitists  - 1));

		// Initialize a ligand validator.
		const validator v(max_rotatable_bonds, max_atoms, max_heavy_atoms, max_hb_donors, max_hb_acceptors, max_mw, max_logp, min_logp);

		// Initialize the number of failures. The program will stop if num_failures reaches max_failures.
		size_t num_failures = 0;

		// The number of ligands (i.e. population size) is equal to the number of elitists plus mutants plus children.
		const size_t num_ligands = num_elitists + num_mutants + num_crossovers;

		// Initialize a pointer vector to dynamically hold and destroy generated ligands.
		using namespace boost;
		ptr_vector<ligand> ligands(num_ligands);

		// Initialize constant strings.
		const string pdbqt_extension_string = ".pdbqt";
		const string ligand_folder_string = "ligand";
		const string output_folder_string = "output";
		const string maximum_failures_reached_string = "The number of failures has reached " + lexical_cast<string>(max_failures);
		const string docking_failed_string = "Docking failed with exit code ";
		const char comma = ',';

		// Initialize process context.
		const boost::process::context ctx;

		// Initialize arguments to docking program.
		vector<string> docking_args;
		if (idock)
		{
			// Initialize argument to idock.
			docking_args.resize(12);
			docking_args[3]  = lexical_cast<string>(seed); // idock supports 64-bit seed.
			docking_args[4]  = "--ligand_folder";
			docking_args[6]  = "--output_folder";
			docking_args[8]  = "--log";
			docking_args[10] = "--csv";
		}
		else
		{
			// Initialize argument to vina.
			docking_args.resize(8);
			docking_args[3] = lexical_cast<string>(int(seed)); // AutoDock Vina does not support 64-bit seed.
			docking_args[4] = "--ligand";
			docking_args[6] = "--out";
		}
		docking_args[0] = "--config";
		docking_args[1] = docking_config_path.string();
		docking_args[2] = "--seed";

		// Initialize csv file for dumping statistics.
		ofstream csv(csv_path);
		csv << "generation,ligand,parent 1,connector 1,parent 2,connector 2,efficacy,free energy in kcal/mol,no. of rotatable bonds,no. of atoms,no. of heavy atoms,no. of hydrogen bond donors,no. of hydrogen bond acceptors,molecular weight,logP\n";

		for (size_t generation = 1; generation <= num_generations; ++generation)
		{
			log << "Running generation " << generation << '\n';

			// Initialize the paths to current generation folder and its two subfolders.
			const path generation_folder(output_folder_path / lexical_cast<string>(generation));
			const path ligand_folder(generation_folder / ligand_folder_string);
			const path output_folder(generation_folder / output_folder_string);

			// Create a new folder and two subfolders for current generation.
			create_directory(generation_folder);
			create_directory(ligand_folder);
			create_directory(output_folder);

			// Create ligands and save them into the ligand folder.
			if (generation == 1)
			{
				// Parse the initial ligand.
				ligands.push_back(new ligand(initial_ligand_path));
				ligand& initial_ligand = ligands.back();

				// Set the parent of 1/1.pdbqt to the initial ligand.
				initial_ligand.parent1 = initial_ligand_path;

				// Save the initial ligand as 1.pdbqt into the ligand folder of generation 1.
				initial_ligand.save(ligand_folder / "1.pdbqt");

				// Create mutants in parallel.
				for (size_t i = 2; i <= num_ligands; ++i)
				{
					// Create a child ligand by mutation.
					do
					{
						ligands.push_back(new ligand(initial_ligand, ligand_flyweight(fragments[uniform_fragment_gen()]), eng, operation_mutation));
						if (v(ligands.back())) break;
						ligands.pop_back();
						if (num_failures++ == max_failures)
						{
							log << maximum_failures_reached_string << '\n';
							return 1;
						}
					} while (true);

					// Save the newly created child ligand.
					ligands.back().save(ligand_folder / (lexical_cast<string>(i) + pdbqt_extension_string));
				}
			}
			else
			{
				// TODO: abstract into tasks for parallel execution.
				// For the elitists, update its parents.
				for (size_t i = 0; i < num_elitists; ++i)
				{
					ligand& l = ligands[i];
					l.parent1 = l.p;
					l.connector1 = 0;
					l.parent2.clear();
					l.connector2 = 0;
					l.p = ligand_folder / (lexical_cast<string>(i + 1) + pdbqt_extension_string);
				}

				// TODO: abstract into tasks for parallel execution.
				// Create child ligands by mutating an elitist with a fragment.
				for (size_t i = num_elitists; i < num_elitists + num_mutants; ++i)
				{					
					do
					{
						ligands.replace(i, new ligand(ligands[uniform_elitist_gen()], ligand_flyweight(fragments[uniform_fragment_gen()]), eng, operation_mutation));
						if (v(ligands[i])) break;
						if (num_failures++ == max_failures)
						{
							log << maximum_failures_reached_string << '\n';
							return 1;
						}
					} while (true);

					// Save the newly created child ligand.
					ligands.back().save(ligand_folder / (lexical_cast<string>(i + 1) + pdbqt_extension_string));
				}

				// TODO: abstract into tasks for parallel execution.
				// Create child ligands by crossovering two elitists.
				for (size_t i = num_elitists + num_mutants; i < num_ligands; ++i)
				{					
					do
					{
						ligands.replace(i, new ligand(ligands[uniform_elitist_gen()], ligands[uniform_elitist_gen()], eng, operation_crossover));
						if (v(ligands[i])) break;
						if (num_failures++ == max_failures)
						{
							log << maximum_failures_reached_string << '\n';
							return 1;
						}
					} while (true);

					// Save the newly created child ligand.
					ligands.back().save(ligand_folder / (lexical_cast<string>(i + 1) + pdbqt_extension_string));
				}
			}

			// Call the docking program to dock ligands.
			if (idock)
			{
				// Invoke idock.
				log << "Calling idock to dock " << (generation == 1 ? num_ligands : num_mutants + num_crossovers) << " ligands\n";
				docking_args[5]  = ligand_folder.string();
				docking_args[7]  = output_folder.string();
				docking_args[9]  = (generation_folder / default_log_path).string();
				docking_args[11] = (generation_folder / default_csv_path).string();
				const int exit_code = create_child(docking_program_path.string(), docking_args, ctx).wait();
				if (exit_code)
				{
					log << docking_failed_string << exit_code << '\n';
					return 1;
				}
			}
			else
			{
				// Invoke vina.
				log << "Calling vina to dock " << (generation == 1 ? num_ligands : num_mutants + num_crossovers) << " ligands\n";
				for (size_t i = (generation == 1 ? 1 : num_elitists + 1); i <= num_ligands; ++i)
				{
					const string filename = lexical_cast<string>(i) + pdbqt_extension_string;
					docking_args[5] = (ligand_folder / filename).string();
					docking_args[7] = (output_folder / filename).string();
					const int exit_code = create_child(docking_program_path.string(), docking_args, ctx).wait();
					if (exit_code)
					{
						log << docking_failed_string << exit_code << '\n';
						return 1;
					}
				}
			}

			// Parse output ligands to obtain predicted free energy and docked coordinates.
			for (size_t i = (generation == 1 ? 0 : num_elitists); i < num_ligands; ++i)
			{
				ligand& l = ligands[i];
				string line;
				ifstream in(output_folder / (lexical_cast<string>(i + 1) + pdbqt_extension_string));
				getline(in, line); // MODEL        1 or MODEL 1
				getline(in, line); // REMARK     FREE ENERGY PREDICATED BY IDOCK:    -4.07 KCAL/MOL or REMARK VINA RESULT:      -9.8      0.000      0.000
				l.evaluate_efficacy(idock ? right_cast<fl>(line, 45, 52) : right_cast<fl>(line, 21, 29));				
				for (size_t i = 0; true;)
				{
					getline(in, line);
					if (line[0] == 'A') // ATOM
					{
						BOOST_ASSERT(l.atoms[i].number == right_cast<size_t>(line, 7, 11));
						l.atoms[i++].coordinate = vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54));						
					}
					else if (line[0] == 'T') // TORSDOF
					{
						BOOST_ASSERT(i == l.atoms.size());
						break;					
					}
				}
				in.close();
			}

			// Sort ligands in ascending order of efficacy.
			ligands.sort();

			// Write summaries to csv.
			for (size_t i = 0; i < num_ligands; ++i)
			{
				const ligand& l = ligands[i];
				csv << generation
					<< comma << l.p
					<< comma << l.parent1
					<< comma << l.connector1
					<< comma << l.parent2
					<< comma << l.connector2
					<< comma << l.efficacy
					<< comma << l.free_energy
					<< comma << l.num_rotatable_bonds
					<< comma << l.num_atoms
					<< comma << l.num_heavy_atoms
					<< comma << l.num_hb_donors
					<< comma << l.num_hb_acceptors
					<< comma << l.mw
					<< comma << l.logp
					<< '\n';
			}
		}
	}
	catch (const std::exception& e)
	{
		std::cerr << e.what() << '\n';
		return 1;
	}
}
