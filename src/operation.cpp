#include <random>
#include "operation.hpp"

void operation::addition_task(const size_t index, const path& p, const size_t seed)
{
	// Initialize a Mersenne Twister random number generator.
	mt19937_64 eng(seed);
	uniform_int_distribution<size_t> uniform_elitist(0, num_elitists - 1);
	uniform_int_distribution<size_t> uniform_fragment(0, num_fragments - 1);

	// Create a child ligand by addition.
	do
	{
		// Obtain references to the two parent ligands.
		ligand& l1 = ligands[uniform_elitist(eng)];
		ligand l2 = ligand_flyweight(fragments[uniform_fragment(eng)]);
		while (!(l1.addition_feasible() && l2.addition_feasible()))
		{
			l1 = ligands[uniform_elitist(eng)];
			l2 = ligand_flyweight(fragments[uniform_fragment(eng)]);
		}

		// Obtain a random mutable atom from the two parent ligands respectively.
		const size_t g1 = uniform_int_distribution<size_t>(0, l1.mutable_atoms.size() - 1)(eng);
		const size_t g2 = uniform_int_distribution<size_t>(0, l2.mutable_atoms.size() - 1)(eng);

		ligands.replace(index, new ligand(p, l1, l2, g1, g2));
		if (v(ligands[index]))
		{
			// Save the newly created child ligand.
			ligands[index].save();
			break;
		}
	} while (++num_failures < max_failures);
}

void operation::subtraction_task(const size_t index, const path& p, const size_t seed)
{
	// Initialize a Mersenne Twister random number generator.
	mt19937_64 eng(seed);
	uniform_int_distribution<size_t> uniform_elitist(0, num_elitists - 1);

	// Create a child ligand by subtraction.
	do
	{
		// Obtain reference to the parent ligand.
		ligand& l1 = ligands[uniform_elitist(eng)];
		while (!l1.subtraction_feasible())
		{
			l1 = ligands[uniform_elitist(eng)];
		}

		// Obtain a random mutable atom from the two parent ligands respectively.
		const size_t g1 = uniform_int_distribution<size_t>(1, l1.num_rotatable_bonds)(eng);

		ligands.replace(index, new ligand(p, l1, g1));
		if (v(ligands[index]))
		{
			// Save the newly created child ligand.
			ligands[index].save();
			break;
		}
	} while (++num_failures < max_failures);
}

void operation::crossover_task(const size_t index, const path& p, const size_t seed)
{
	// Initialize a Mersenne Twister random number generator.
	mt19937_64 eng(seed);
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

		ligands.replace(index, new ligand(p, l1, l2, g1, g2, true));
		if (v(ligands[index]))
		{
			// Save the newly created child ligand.
			ligands[index].save();
			break;
		}
	} while (++num_failures < max_failures);
}
