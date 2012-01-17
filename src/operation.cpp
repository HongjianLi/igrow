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

#include "operation.hpp"

namespace igrow
{
	void operation::mutation_task(const size_t index, const path& ligand_path, const path& output_path, const size_t seed, const size_t num_elitists)
	{
		// Initialize a Mersenne Twister random number generator.
		mt19937eng eng(seed);

		// Initialize random number generators for obtaining a random fragment and a random elitist.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_elitist_gen(eng, uniform_int_distribution<size_t>(0, num_elitists - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_fragment_gen(eng, uniform_int_distribution<size_t>(0, num_fragments - 1));

		// Create a child ligand by mutation.
		do
		{
			// Obtain constant references to the two parent ligands.
			const ligand& l1 = ligands[uniform_elitist_gen()];
			const ligand& l2 = ligand_flyweight(fragments[uniform_fragment_gen()]);

			// Obtain a random mutable atom from the two parent ligands respectively.
			const size_t g1 = variate_generator<mt19937eng, uniform_int_distribution<size_t>>(eng, uniform_int_distribution<size_t>(0, l1.mutable_atoms.size() - 1))();
			const size_t g2 = variate_generator<mt19937eng, uniform_int_distribution<size_t>>(eng, uniform_int_distribution<size_t>(0, l2.mutable_atoms.size() - 1))();

			ligands.replace(index, new ligand(l1, l2, g1, g2));
			if (v(ligands[index])) break;
			if (num_failures++ >= max_failures) throw maximum_failures_reached_error(max_failures);
		} while (true);

		// Save the newly created child ligand.
		ligand& l = ligands[index];
		l.save(ligand_path);
		l.p =  output_path;
	}

	void operation::crossover_task(const size_t index, const path& ligand_path, const path& output_path, const size_t seed, const size_t num_elitists)
	{
		// Initialize a Mersenne Twister random number generator.
		mt19937eng eng(seed);

		// Initialize random number generators for obtaining a random fragment and a random elitist.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_elitist_gen(eng, uniform_int_distribution<size_t>(0, num_elitists - 1));

		// Create a child ligand by mutation.
		do
		{
			// Obtain constant references to the two parent ligands.
			const ligand& l1 = ligands[uniform_elitist_gen()];
			const ligand& l2 = ligands[uniform_elitist_gen()];

			// Obtain a random mutable atom from the two parent ligands respectively.
			const size_t f1 = variate_generator<mt19937eng, uniform_int_distribution<size_t>>(eng, uniform_int_distribution<size_t>(0, l1.num_rotatable_bonds - 1))();
			const size_t f2 = variate_generator<mt19937eng, uniform_int_distribution<size_t>>(eng, uniform_int_distribution<size_t>(0, l2.num_rotatable_bonds - 1))();
			const size_t g1 = variate_generator<mt19937eng, uniform_int_distribution<size_t>>(eng, uniform_int_distribution<size_t>(0, 1))();
			const size_t g2 = variate_generator<mt19937eng, uniform_int_distribution<size_t>>(eng, uniform_int_distribution<size_t>(0, 1))();

			ligands.replace(index, new ligand(l1, l2, f1, f2, g1, g2));
			if (v(ligands[index])) break;
			if (num_failures++ >= max_failures) throw maximum_failures_reached_error(max_failures);
		} while (true);

		// Save the newly created child ligand.
		ligand& l = ligands[index];
		l.save(ligand_path);
		l.p =  output_path;
	}
}
