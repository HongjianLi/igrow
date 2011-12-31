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
		vector<size_t> branches; ///< Child branches.
		vector<atom> atoms; ///< Heavy atoms.

		/// Constructs a frame, and relates it to its parent frame.
		explicit frame(const size_t parent) : parent(parent)
		{
			branches.reserve(5); // A frame typically consists of <= 5 branch frames.
			atoms.reserve(20); // A frame typically consists of <= 20 atoms.
		}

		/// Copy constructor.
		frame(const frame& f) : parent(f.parent), rotorX(f.rotorX), branches(f.branches), atoms(f.atoms) {}

		/// Move constructor.
		frame(frame&& f) : parent(f.parent), rotorX(f.rotorX), branches(static_cast<vector<size_t>&&>(f.branches)), atoms(static_cast<vector<atom>&&>(f.atoms)) {}
	};

	using boost::filesystem::path;

	// Represents a ligand.
	class ligand
	{
	public:
		path p; ///< The path to the current ligand.
		path parent1; ///< The first parent ligand to synthesize the current ligand.
		path parent2; ///< The second parent ligand, if any, to synthesize the current ligand.
		size_t connector1; ///< The serial number of the connecting atom of parent 1.
		size_t connector2; ///< The serial number of the connecting atom of parent 2.
		vector<frame> frames; ///< Ligand frames.
		vector<atom_index> mutable_atoms; ///< Hydrogens or halogens.
		size_t num_heavy_atoms; ///< Number of heavy atoms.
		size_t num_hb_donors; ///< Number of hydrogen bond donors.
		size_t num_hb_acceptors; ///< Number of hydrogen bond acceptors.
		fl mw; ///< Molecular weight.
		fl logp; ///< Predicted LogP obtained by external XLOGP3.
		fl free_energy; ///< Predicted free energy obtained by external docking.
		fl efficacy; ///< Ligand efficacy
		bool mutation_feasible;  // True if the current ligand is able to perform mutation.
		bool crossover_feasible; // True if the current ligand is able to perform crossover.
		
		/// Constructs a ligand by parsing a given ligand file in pdbqt.
		explicit ligand(const path& p);

		/// Updates the path and saves the current ligand to a file in pdbqt format.
		void save(const path& p);

		/// Mutates the current ligand.
		ligand* mutate(const ligand& other, const mt19937eng& eng) const;

		/// Recalculates ligand efficacy, defined as free_energy / num_heavy_atoms. This definition contradicts our conventional definition, but it works fine for sorting ligands.
		void evaluate_efficacy(const fl free_energy);

		/// For sorting ptr_vector<ligand>.
		const bool operator<(const ligand& l) const
		{
			return efficacy < l.efficacy;
		}
	
	private:
		/// Constructs an empty ligand.
		ligand() {}
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
		validator(const size_t max_heavy_atoms, const size_t max_hb_donors, const size_t max_hb_acceptors, const fl max_mw, const fl max_logp, const fl min_logp) : max_heavy_atoms(max_heavy_atoms), max_hb_donors(max_hb_donors), max_hb_acceptors(max_hb_acceptors), max_mw(max_mw), max_logp(max_logp), min_logp(min_logp) {}

		const bool operator()(const ligand& l) const
		{
			if (l.num_heavy_atoms > max_heavy_atoms) return false;
			if (l.num_hb_donors > max_hb_donors) return false;
			if (l.num_hb_acceptors > max_hb_acceptors) return false;
			if (l.mw > max_mw) return false;
			if (l.logp > max_logp) return false;
			if (l.logp < min_logp) return false;
			return true;
		}

	private:
		const size_t max_heavy_atoms;
		const size_t max_hb_donors;
		const size_t max_hb_acceptors;
		const fl max_mw;
		const fl max_logp;
		const fl min_logp;
	};
}

#endif
