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

#pragma once
#ifndef IGROW_LIGAND_HPP
#define IGROW_LIGAND_HPP

#include <boost/filesystem/path.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include "atom.hpp"

namespace igrow
{

// Choose the appropriate Mersenne Twister engine for random number generation on 32-bit or 64-bit platform.
#if defined(__x86_64) || defined(__x86_64__) || defined(__amd64) || defined(__amd64__) || defined(_M_X64) || defined(_M_AMD64)
	typedef boost::random::mt19937_64 mt19937eng;
#else
	typedef boost::random::mt19937 mt19937eng;
#endif

	/// Represents a ROOT or a BRANCH in PDBQT structure.
	class frame
	{
	public:
		size_t parent; ///< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
		size_t rotorX; ///< Index pointing to the parent frame atom which forms a rotatable bond with the first atom of current frame, a.k.a. rotor Y.
		size_t begin; ///< The inclusive beginning index to the atoms of the current frame.
		size_t end; ///< The exclusive ending index to the atoms of the current frame.
		vector<size_t> branches; ///< Child branches.

		/// Constructs a frame, and initializes its parent frame and beginning atom index.
		explicit frame(const size_t parent, const size_t rotorX, const size_t begin) : parent(parent), rotorX(rotorX), begin(begin)
		{
			BOOST_ASSERT(rotorX <= begin); // The equal sign holds only for ROOT frame.
			branches.reserve(4); // A frame typically consists of <= 4 branch frames.
		}

		/// Copy constructor.
		frame(const frame& f) : parent(f.parent), rotorX(f.rotorX), begin(f.begin), end(f.end), branches(f.branches) {}

		/// Move constructor.
		frame(frame&& f) : parent(f.parent), rotorX(f.rotorX), begin(f.begin), end(f.end), branches(static_cast<vector<size_t>&&>(f.branches)) {}
	};

	enum operation
	{
		operation_mutation,
		operation_crossover
	};

	using boost::filesystem::path;

	/// Represents a ligand.
	class ligand
	{
	public:
		path p; ///< The path to the current ligand.
		path parent1; ///< The first parent ligand to synthesize the current ligand.
		path parent2; ///< The second parent ligand, if any, to synthesize the current ligand.
		size_t connector1; ///< The serial number of the connecting atom of parent 1.
		size_t connector2; ///< The serial number of the connecting atom of parent 2.
		vector<frame> frames; ///< Frames.
		vector<atom> atoms; ///< Atoms.
		vector<size_t> mutable_atoms; ///< Hydrogens or halogens.
		size_t num_rotatable_bonds; ///< Number of rotatable bonds.
		size_t num_atoms; ///< Number of atoms.
		size_t num_heavy_atoms; ///< Number of heavy atoms.
		size_t num_hb_donors; ///< Number of hydrogen bond donors.
		size_t num_hb_acceptors; ///< Number of hydrogen bond acceptors.
		fl mw; ///< Molecular weight.
		fl logp; ///< Predicted LogP obtained by external XLOGP3.
		fl free_energy; ///< Predicted free energy obtained by external docking.
		fl efficacy; ///< Ligand efficacy

		/// Constructs a ligand by parsing a given ligand file in pdbqt.
		explicit ligand(const path& p);

		/// Constructs a ligand by either mutation or crossover.
		explicit ligand(const ligand& l1, const ligand& l2, const mt19937eng& eng, const operation op);

		/// Updates the path and saves the current ligand to a file in pdbqt format.
		void save(const path& p);

		/// Returns true if the current ligand is able to perform mutation.
		bool mutation_feasible() const
		{
			return mutable_atoms.size() > 0;
		}

		/// Returns true if the current ligand is able to perform crossover.
		bool crossover_feasible() const
		{
			return num_rotatable_bonds > 0;
		}

		/// Recalculates ligand efficacy, defined as free_energy / num_heavy_atoms. This definition contradicts our conventional definition, but it works fine for sorting ligands.
		void evaluate_efficacy(const fl free_energy)
		{
			this->free_energy = free_energy;
			efficacy = free_energy / num_heavy_atoms;
		}

		/// Compares the efficacy of the current ligand and the other ligand for sorting ptr_vector<ligand>.
		const bool operator<(const ligand& l) const
		{
			return efficacy < l.efficacy;
		}

		/// Gets the frame to which a given atom belongs to.
		const size_t get_frame(const size_t atom_idx) const
		{
			BOOST_ASSERT(num_rotatable_bonds == frames.size() - 1);
			for (size_t i = 0; i < num_rotatable_bonds; ++i)
			{
				if (atom_idx < frames[i].end) return i;
			}
			return num_rotatable_bonds;
		}

	private:
		/// Mutates ligand 1 with ligand 2.
		void mutate(const ligand& l1, const ligand& l2, const mt19937eng& eng);

		/// Crossovers ligand 1 with ligand 2.
		void crossover(const ligand& l1, const ligand& l2, const mt19937eng& eng);
	};

	/// For extracting the path out of a ligand.
	class ligand_path_extractor
	{
	public:
		const path& operator()(const ligand& l) const
		{
			return l.p;
		}
	};

	/// Define flyweight type for ligand.
	using namespace boost::flyweights;
	typedef	flyweight<key_value<path, ligand, ligand_path_extractor>, no_tracking> ligand_flyweight;

	/// Represents a ligand validator.
	class validator
	{
	public:
		validator(const size_t max_rotatable_bonds, const size_t max_atoms, const size_t max_heavy_atoms, const size_t max_hb_donors, const size_t max_hb_acceptors, const fl max_mw, const fl max_logp, const fl min_logp) : max_rotatable_bonds(max_rotatable_bonds), max_atoms(max_atoms), max_heavy_atoms(max_heavy_atoms), max_hb_donors(max_hb_donors), max_hb_acceptors(max_hb_acceptors), max_mw(max_mw), max_logp(max_logp), min_logp(min_logp) {}

		const bool operator()(const ligand& l) const
		{
			if (l.num_rotatable_bonds > max_rotatable_bonds) return false;
			if (l.num_atoms > max_atoms) return false;
			if (l.num_heavy_atoms > max_heavy_atoms) return false;
			if (l.num_hb_donors > max_hb_donors) return false;
			if (l.num_hb_acceptors > max_hb_acceptors) return false;
			if (l.mw > max_mw) return false;
			if (l.logp > max_logp) return false;
			if (l.logp < min_logp) return false;
			return true;
		}

	private:
		const size_t max_rotatable_bonds;
		const size_t max_atoms;
		const size_t max_heavy_atoms;
		const size_t max_hb_donors;
		const size_t max_hb_acceptors;
		const fl max_mw;
		const fl max_logp;
		const fl min_logp;
	};
}

#endif
