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

#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include "atom.hpp"

namespace igrow
{
	/// Represents a parsing error.
	class parsing_error : public std::domain_error
	{
	public:
		/// Constructs a parsing error.
		parsing_error(const path& file, const size_t line, const string& reason) : std::domain_error("Error parsing \"" + file.filename().string() + "\" on line " + boost::lexical_cast<string>(line) + ": " + reason) {}
	};

	/// Represents a ROOT or a BRANCH in PDBQT structure.
	class frame
	{
	public:
		size_t parent; ///< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
		size_t rotorX; ///< Serial number of the parent frame atom which forms a rotatable bond with rotorY.
		size_t rotorY; ///< Serial number of the current frame atom which forms a rotatable bond with rotorX.
		size_t begin; ///< The inclusive beginning index to the atoms of the current frame.
		size_t end; ///< The exclusive ending index to the atoms of the current frame.
		vector<size_t> branches; ///< Indexes to child branches.

		/// Constructs a frame, and initializes its parent frame, rotor connectors, and beginning atom index.
		explicit frame(const size_t parent, const size_t rotorX, const size_t rotorY, const size_t begin) : parent(parent), rotorX(rotorX), rotorY(rotorY), begin(begin) {}

		/// Copy constructor.
		frame(const frame& f) : parent(f.parent), rotorX(f.rotorX), rotorY(f.rotorY), begin(f.begin), end(f.end), branches(f.branches) {}

		/// Move constructor.
		frame(frame&& f) : parent(f.parent), rotorX(f.rotorX), rotorY(f.rotorY), begin(f.begin), end(f.end), branches(static_cast<vector<size_t>&&>(f.branches)) {}
		
		/// Copy assignment.
		frame& operator=(const frame& f)
		{
			if (this != &f)
			{
				parent = f.parent;
				rotorX = f.rotorX;
				rotorY = f.rotorY;
				begin = f.begin;
				end = f.end;
				branches = f.branches;
			}
			return *this;
		}	
		
		/// Move assignment.
		frame& operator=(frame&& f)
		{
			if (this != &f)
			{
				parent = f.parent;
				rotorX = f.rotorX;
				rotorY = f.rotorY;
				begin = f.begin;
				end = f.end;
				branches = static_cast<vector<size_t>&&>(f.branches);
			}
			return *this;
		}
	};

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
		size_t max_atom_number; ///< Maximum atom serial number.
		size_t num_rotatable_bonds; ///< Number of rotatable bonds.
		size_t num_atoms; ///< Number of atoms.
		size_t num_heavy_atoms; ///< Number of heavy atoms.
		size_t num_hb_donors; ///< Number of hydrogen bond donors.
		size_t num_hb_acceptors; ///< Number of hydrogen bond acceptors.
		fl mw; ///< Molecular weight.
		fl logp; ///< Predicted LogP obtained by external XLOGP3.
		fl free_energy; ///< Predicted free energy obtained by external docking.
		fl efficacy; ///< Ligand efficacy
		
		explicit ligand() {}

		/// Constructs a ligand by parsing a given ligand file in pdbqt.
		/// @exception parsing_error Thrown when error parsing the ligand file.
		explicit ligand(const path& p);

		/// Constructs a ligand by mutation.
		explicit ligand(const ligand& l1, const ligand& l2, const size_t g1, const size_t g2);

		/// Constructs a ligand by crossover.
		explicit ligand(const ligand& l1, const ligand& l2, const size_t f1idx, const size_t f2idx, const size_t g1, const size_t g2);
		
		/// Saves the current ligand to a file in pdbqt format.
		void save(const path& p) const;

		/// Gets the frame and index to which a atom belongs to given its serial number.
		std::pair<size_t, size_t> get_frame(const size_t srn) const;

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

		/// Evaluates ligand efficacy from free energy.
		void evaluate_efficacy(const fl free_energy)
		{
			this->free_energy = free_energy;
			efficacy = free_energy * pow(static_cast<fl>(num_heavy_atoms), static_cast<fl>(-0.4));
		}

		/// Compares the efficacy of the current ligand and the other ligand for sorting ptr_vector<ligand>.
		const bool operator<(const ligand& l) const
		{
			return efficacy < l.efficacy;
		}
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
