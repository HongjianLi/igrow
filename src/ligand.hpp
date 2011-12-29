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
	/// Represents a ROOT or a BRANCH in PDBQT structure.
	class frame
	{
	public:
		size_t parent; ///< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
		size_t rotorX; ///< Index pointing to the parent frame atom which forms a rotatable bond with the first atom of current frame, a.k.a. rotor Y.
		vector<size_t> branches; ///< Child branches.
		vector<atom> atoms; ///< Heavy atoms.
		vector<mutation_point> mutation_points; ///< Hydrogens or halogens.

		/// Constructs a frame, and relates it to its parent frame.
		explicit frame(const size_t parent) : parent(parent)
		{
			branches.reserve(5); // A frame typically consists of <= 5 branch frames.
			atoms.reserve(20); // A frame typically consists of <= 20 atoms.
			mutation_points.reserve(5); // A frame typically consists of <= 5 mutation points.
		}
		
		/// Copy constructor.
		frame(const frame& f) : parent(f.parent), rotorX(f.rotorX), branches(f.branches), atoms(f.atoms), mutation_points(f.mutation_points) {}
		
		/// Move constructor.
		frame(frame&& f) : parent(f.parent), rotorX(f.rotorX), branches(static_cast<vector<size_t>&&>(f.branches)), atoms(static_cast<vector<atom>&&>(f.atoms)), mutation_points(static_cast<vector<mutation_point>&&>(f.mutation_points)) {}
	};

	using boost::filesystem::path;

	// Represents a ligand.

	class ligand
	{
	public:
		const path p; ///< The path to the fragment.
		vector<frame> frames; ///< Ligand frames.
		size_t num_heavy_atoms; ///< Number of heavy atoms.
		size_t num_hb_donors; ///< Number of hydrogen bond donors.
		size_t num_hb_acceptors; ///< Number of hydrogen bond acceptors.
		fl mw; ///< Molecular weight.		
		fl logp; ///< Predicted LogP obtained by external XLOGP3.
		fl free_energy; ///< Predicted free energy obtained by external docking.
		fl efficacy; ///< Ligand efficacy

		ligand() {}
		explicit ligand(const path& p);
		
		/// Saves the current ligand to a file in pdbqt format.
		void save(const path& p) const;

		/// Mutates the current ligand.
		ligand* mutate(const ligand& lig) const;

		/// Recalculates ligand efficacy, defined as free_energy / num_heavy_atoms.
		void evaluate_efficacy();

		/// For sorting ptr_vector<ligand>.
		const bool operator<(const ligand& lig) const
		{
			return efficacy < lig.efficacy;
		}
	};

	/// For extracting the path out of a ligand.
	class ligand_path_extractor
	{
	public:
		const path& operator()(const ligand& lig) const
		{
			return lig.p;
		}
	};

	/// Define flyweight type for ligand.
	using namespace boost::flyweights;
	typedef	flyweight<key_value<path, ligand, ligand_path_extractor>, no_tracking> ligand_flyweight;
}

#endif
