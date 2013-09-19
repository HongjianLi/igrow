/*

   Copyright (c) 2012, The Chinese University of Hong Kong

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
using boost::filesystem::path;

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

	/// Copy assignment operator.
	frame& operator=(const frame& f)
	{
		this->parent = f.parent;
		this->rotorX = f.rotorX;
		this->rotorY = f.rotorY;
		this->begin = f.begin;
		this->end = f.end;
		this->branches = f.branches;
		return *this;
	}

	/// Move assignment operator.
	frame& operator=(frame&& f)
	{
		this->parent = f.parent;
		this->rotorX = f.rotorX;
		this->rotorY = f.rotorY;
		this->begin = f.begin;
		this->end = f.end;
		this->branches = static_cast<vector<size_t>&&>(f.branches);
		return *this;
	}
};

/// Represents a ligand.
class ligand
{
public:
	path p; ///< Path to the current ligand.
	path parent1; ///< Parent ligand 1.
	path parent2; ///< Parent ligand 2.
	string connector1; ///< The connecting bond of parent 1.
	string connector2; ///< The connecting bond of parent 2.
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
	fl fe; ///< Predicted free energy obtained by external docking.
	fl le; ///< Predicted ligand efficiency obtained by external docking.
	explicit ligand() {}

	/// Constructs a ligand by parsing a given ligand file in PDBQT.
	/// @exception parsing_error Thrown when error parsing the ligand file.
	explicit ligand(const path& p);

	/// Constructs a ligand by addition.
	explicit ligand(const path& p, const ligand& l1, const ligand& l2, const size_t g1, const size_t g2);

	/// Constructs a ligand by subtraction.
	explicit ligand(const path& p, const ligand& l1, const size_t g1);

	/// Constructs a ligand by crossover.
	explicit ligand(const path& p, const ligand& l1, const ligand& l2, const size_t g1, const size_t g2, const bool dummy);

	/// Saves the current ligand to a file in PDBQT format.
	void save() const;

	/// Parse the docked ligand to obtain predicted free energy and docked coordinates.
	void update(const path& p);

	/// Gets the frame and index to which a atom belongs to given its serial number.
	pair<size_t, size_t> get_frame(const size_t srn) const;

	/// Returns true if the current ligand is able to perform addition.
	bool addition_feasible() const
	{
		return mutable_atoms.size() > 0;
	}

	/// Returns true if the current ligand is able to perform subtraction.
	bool subtraction_feasible() const
	{
		return num_rotatable_bonds > 0;
	}

	/// Returns true if the current ligand is able to perform crossover.
	bool crossover_feasible() const
	{
		return num_rotatable_bonds > 0;
	}

	/// Compares the efficacy of the current ligand and the other ligand for sorting ptr_vector<ligand>.
	bool operator<(const ligand& l) const
	{
		return fe < l.fe;
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
	validator(const size_t max_rotatable_bonds, const size_t max_atoms, const size_t max_heavy_atoms, const size_t max_hb_donors, const size_t max_hb_acceptors, const fl max_mw) : max_rotatable_bonds(max_rotatable_bonds), max_atoms(max_atoms), max_heavy_atoms(max_heavy_atoms), max_hb_donors(max_hb_donors), max_hb_acceptors(max_hb_acceptors), max_mw(max_mw) {}

	bool operator()(const ligand& l) const
	{
		if (l.num_rotatable_bonds > max_rotatable_bonds) return false;
		if (l.num_atoms > max_atoms) return false;
		if (l.num_heavy_atoms > max_heavy_atoms) return false;
		if (l.num_hb_donors > max_hb_donors) return false;
		if (l.num_hb_acceptors > max_hb_acceptors) return false;
		if (l.mw > max_mw) return false;
		return true;
	}

private:
	const size_t max_rotatable_bonds;
	const size_t max_atoms;
	const size_t max_heavy_atoms;
	const size_t max_hb_donors;
	const size_t max_hb_acceptors;
	const fl max_mw;
};

#endif
