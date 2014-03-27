#pragma once
#ifndef IGROW_LIGAND_HPP
#define IGROW_LIGAND_HPP

#include <boost/filesystem/path.hpp>
#include "atom.hpp"
using namespace boost::filesystem;

//! Represents a ROOT or a BRANCH in PDBQT structure.
class frame
{
public:
	size_t parent; //!< Frame array index pointing to the parent of current frame. For ROOT frame, this field is not used.
	size_t rotorX; //!< Serial number of the parent frame atom which forms a rotatable bond with rotorY.
	size_t rotorY; //!< Serial number of the current frame atom which forms a rotatable bond with rotorX.
	size_t begin; //!< The inclusive beginning index to the atoms of the current frame.
	size_t end; //!< The exclusive ending index to the atoms of the current frame.
	vector<size_t> branches; //!< Indexes to child branches.

	//! Constructs a frame, and initializes its parent frame, rotor connectors, and beginning atom index.
	explicit frame(const size_t parent, const size_t rotorX, const size_t rotorY, const size_t begin) : parent(parent), rotorX(rotorX), rotorY(rotorY), begin(begin) {}
};

//! Represents a ligand.
class ligand
{
public:
	path p; //!< Path to the current ligand.
	path parent1; //!< Parent ligand 1.
	path parent2; //!< Parent ligand 2.
	string connector1; //!< The connecting bond of parent 1.
	string connector2; //!< The connecting bond of parent 2.
	vector<frame> frames; //!< Frames.
	vector<atom> atoms; //!< Atoms.
	size_t max_atom_number; //!< Maximum atom serial number.
	size_t num_rotatable_bonds; //!< Number of rotatable bonds.
	size_t num_atoms; //!< Number of atoms.
	size_t num_hb_donors; //!< Number of hydrogen bond donors.
	size_t num_hb_acceptors; //!< Number of hydrogen bond acceptors.
	double mm; //!< Molecular mass.
	double fe; //!< Predicted free energy obtained by external docking.

	//! Constructs an empty ligand.
	explicit ligand() {}

	//! Constructs a ligand by parsing a given ligand file in PDBQT.
	explicit ligand(const path& p);

	//! Constructs a ligand by crossover.
	explicit ligand(const path& p, const ligand& l1, const ligand& l2, const size_t g1, const size_t g2);

	//! Returns true if the current ligand is able to perform crossover.
	bool crossover_feasible() const;

	//! Compares the efficacy of the current ligand and the other ligand for sorting ptr_vector<ligand>.
	bool operator<(const ligand& l) const;

	//! Saves the current ligand to a file in PDBQT format.
	void save() const;

	//! Parses the docked ligand to obtain predicted free energy and docked coordinates.
	void update(const path& p);

	//! Gets the frame and index to which a atom belongs to given its serial number.
	pair<size_t, size_t> get_frame(const size_t srn) const;
};

#endif
