#include <algorithm>
#include <cassert>
#include "array.hpp"
#include "atom.hpp"

//! AutoDock4 atom type strings, e.g. H, HD, C, A.
const array<string, atom::n> atom::ad_strings =
{
	"H" , //  0
	"HD", //  1
	"C" , //  2
	"A" , //  3
	"N" , //  4
	"NA", //  5
	"OA", //  6
	"S" , //  7
	"SA", //  8
	"Se", //  9
	"P" , // 10
	"F" , // 11
	"Cl", // 12
	"Br", // 13
	"I" , // 14
	"Zn", // 15
	"Fe", // 16
	"Mg", // 17
	"Ca", // 18
	"Mn", // 19
	"Cu", // 20
	"Na", // 21
	"K" , // 22
	"Hg", // 23
	"Ni", // 24
	"Co", // 25
	"Cd", // 26
	"As", // 27
	"Sr", // 28
	"U" , // 29
	"Cs", // 30
};

//! Covalent radii of AutoDock4 atom types. http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements
const array<double, atom::n> atom::ad_covalent_radii =
{
	0.37, //  0 = H
	0.37, //  1 = HD
	0.77, //  2 = C
	0.77, //  3 = A
	0.75, //  4 = N
	0.75, //  5 = NA
	0.73, //  6 = OA
	1.02, //  7 = S
	1.02, //  8 = SA
	1.16, //  9 = Se
	1.06, // 10 = P
	0.71, // 11 = F
	1.99, // 12 = Cl
	1.14, // 13 = Br
	1.33, // 14 = I
	1.31, // 15 = Zn
	1.25, // 16 = Fe
	1.30, // 17 = Mg
	1.74, // 18 = Ca
	1.39, // 19 = Mn
	1.38, // 20 = Cu
	1.54, // 21 = Na
	1.96, // 22 = K
	1.49, // 23 = Hg
	1.21, // 24 = Ni
	1.26, // 25 = Co
	1.48, // 26 = Cd
	1.19, // 27 = As
	1.92, // 28 = Sr
	1.96, // 29 = U
	2.25, // 30 = Cs
};

//! AutoDock4 atomic masses. http://en.wikipedia.org/wiki/Relative_atomic_mass
const array<double, atom::n> atom::ad_atomic_masses =
{
	  1.008,//  0 = HD
	  1.008,//  1 = H
	 12.01, //  2 = C
	 12.01, //  3 = A
	 14.01, //  4 = N
	 14.01, //  5 = NA
	 16.00, //  6 = OA
	 32.07, //  7 = SA
	 32.07, //  8 = S
	 78.96, //  9 = Se
	 30.97, // 10 = P
	 19.00, // 11 = F
	 35.45, // 12 = Cl
	 79.90, // 13 = Br
	126.90, // 14 = I
	 65.38, // 15 = Zn
	 55.85, // 16 = Fe
	 24.31, // 17 = Mg
	 40.08, // 18 = Ca
	 54.94, // 19 = Mn
	 63.55, // 20 = Cu
	 22.99, // 21 = Na
	 39.10, // 22 = K
	200.59, // 23 = Hg
	 58.69, // 24 = Ni
	 58.93, // 25 = Co
	112.41, // 26 = Cd
	 74.92, // 27 = As
	 87.62, // 28 = Sr
	238.03, // 29 = U
	132.91, // 30 = Cs
};

//! Returns the AutoDock4 atom type of the given string.
size_t atom::parse_ad_string(const string& ad_string)
{
	return find(ad_strings.cbegin(), ad_strings.cend(), ad_string) - ad_strings.cbegin();
}

//! Constructs an atoms.
atom::atom(const string& name, const string& columns_13_to_30, const string columns_55_to_79, const size_t srn, const array<double, 3>& coord, const size_t ad) : name(name), columns_13_to_30(columns_13_to_30), columns_55_to_79(columns_55_to_79), srn(srn), coord(coord), ad(ad)
{
}

//! Returns covalent radius of the current AutoDock4 atom type.
double atom::covalent_radius() const
{
	return ad_covalent_radii[ad];
}

//! Returns atomic mass of the current AutoDock4 atom type.
double atom::atomic_mass() const
{
	return ad_atomic_masses[ad];
}

//! Returns true if the current atom is a hydrogen.
bool atom::is_hydrogen() const
{
	return ad <= 1;
}

//! Returns true is the current atom is a hydrogen bond donor, i.e. polar hydrogen.
bool atom::is_hb_donor() const
{
	return !ad;
}

//! Returns true is the current atom is a hydrogen bond acceptor.
bool atom::is_hb_acceptor() const
{
	return 5 <= ad && ad <= 7;
}

//! Returns true if the current atom is covalently bonded to a given atom.
bool atom::is_neighbor(const atom& a) const
{
	assert(this != &a);
	const double r = 1.1 * (covalent_radius() + a.covalent_radius());
	return distance_sqr(coord, a.coord) < r * r;
}
