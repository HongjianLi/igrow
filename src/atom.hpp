#pragma once
#ifndef IGROW_ATOM_HPP
#define IGROW_ATOM_HPP

#include <array>
#include <string>
using namespace std;

// Represents an atom.
class atom
{
public:
	static const size_t n = 31; //!< Number of AutoDock4 atom types.
	static const array<string, n> ad_strings; //!< AutoDock4 atom type strings, e.g. H, HD, C, A.
	static const array<double, n> ad_covalent_radii; //!< Covalent radii of AutoDock4 atom types.
	static const array<double, n> ad_atomic_masses; //!< Atomic masses of AutoDock4 atom types.
	string name; //!< Atom name;
	string columns_13_to_30; //!< Columns from 1-based [13, 30] of an ATOM/HETATM line in PDBQT format.
	string columns_55_to_79; //!< Columns from 1-based [55, 79] of an ATOM/HETATM line in PDBQT format.
	size_t srn; //!< Serial number.
	array<double, 3> coord; //!< 3D coordinate.
	size_t ad; //!< AutoDock4 atom type.

	//! Returns the AutoDock4 atom type of the given string.
	static size_t parse_ad_string(const string& ad_string);

	//! Constructs an atoms.
	explicit atom(const string& name, const string& columns_13_to_30, const string columns_55_to_79, const size_t srn, const array<double, 3>& coord, const size_t ad);

	//! Returns covalent radius of the current AutoDock4 atom type.
	double covalent_radius() const;

	//! Returns atomic mass of the current AutoDock4 atom type.
	double atomic_mass() const;

	//! Returns true if the current atom is a hydrogen.
	bool is_hydrogen() const;

	//! Returns true is the current atom is a hydrogen bond donor, i.e. polar hydrogen.
	bool is_hb_donor() const;

	//! Returns true is the current atom is a hydrogen bond acceptor.
	bool is_hb_acceptor() const;

	//! Returns true if the current atom is covalently bonded to a given atom.
	bool is_neighbor(const atom& a) const;
};

#endif
