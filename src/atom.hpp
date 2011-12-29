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
#ifndef IGROW_ATOM_HPP
#define IGROW_ATOM_HPP

#include <boost/lexical_cast.hpp>
#include "common.hpp"
#include "vec3.hpp"

namespace igrow
{
	// AutoDock4 atom types.
	const size_t AD_TYPE_HD   =  0;	///< Polar hydrogen, i.e. bonded to a hetero atom.
	const size_t AD_TYPE_H    =  1;	///< Non-polar hydrogen, i.e. bonded to carbon.
	const size_t AD_TYPE_C    =  2; ///< Carbon, not in a ring.
	const size_t AD_TYPE_A    =  3; ///< Carbon, in a ring.
	const size_t AD_TYPE_N    =  4; ///< Nitrogen, not a hydrogen bond acceptor.
	const size_t AD_TYPE_NA   =  5; ///< Nitrogen, a hydrogen bond acceptor.
	const size_t AD_TYPE_OA   =  6; ///< Oxygen, a hydrogen bond acceptor.
	const size_t AD_TYPE_SA   =  7; ///< Sulfur, a hydrogen bond acceptor.
	const size_t AD_TYPE_S    =  8; ///< Sulfur, not a hydrogen bond acceptor.
	const size_t AD_TYPE_Se   =  9; ///< Selenium.
	const size_t AD_TYPE_P    = 10; ///< Phosphorus.
	const size_t AD_TYPE_F    = 11; ///< Fluorine.
	const size_t AD_TYPE_Cl   = 12; ///< Chlorine.
	const size_t AD_TYPE_Br   = 13; ///< Bromine.
	const size_t AD_TYPE_I    = 14; ///< Iodine.
	const size_t AD_TYPE_Zn   = 15; ///< Zine.
	const size_t AD_TYPE_Fe   = 16; ///< Iron.
	const size_t AD_TYPE_Mg   = 17; ///< Magnesium.
	const size_t AD_TYPE_Ca   = 18; ///< Calcium.
	const size_t AD_TYPE_Mn   = 19; ///< Manganese.
	const size_t AD_TYPE_Cu   = 20; ///< Copper.
	const size_t AD_TYPE_Na   = 21; ///< Sodium.
	const size_t AD_TYPE_K    = 22; ///< Potassium.
	const size_t AD_TYPE_Hg   = 23; ///< Mercury.
	const size_t AD_TYPE_Ni   = 24; ///< Nickel.
	const size_t AD_TYPE_Co   = 25; ///< Cobalt.
	const size_t AD_TYPE_Cd   = 26; ///< Cadmium.
	const size_t AD_TYPE_As   = 27; ///< Arsenic.
	const size_t AD_TYPE_SIZE = 28; ///< Number of supported AutoDock4 atom types.

	const string ad_names[] = ///< AutoDock4 atom type names.
	{
		"HD", //  0 = AD_TYPE_HD
		"H" , //  1 = AD_TYPE_H
		"C" , //  2 = AD_TYPE_C
		"A" , //  3 = AD_TYPE_A
		"N" , //  4 = AD_TYPE_N
		"NA", //  5 = AD_TYPE_NA
		"OA", //  6 = AD_TYPE_OA
		"SA", //  7 = AD_TYPE_SA
		"S" , //  8 = AD_TYPE_S
		"Se", //  9 = AD_TYPE_Se
		"P" , // 10 = AD_TYPE_P
		"F" , // 11 = AD_TYPE_F
		"Cl", // 12 = AD_TYPE_Cl
		"Br", // 13 = AD_TYPE_Br
		"I" , // 14 = AD_TYPE_I
		"Zn", // 15 = AD_TYPE_Zn
		"Fe", // 16 = AD_TYPE_Fe
		"Mg", // 17 = AD_TYPE_Mg
		"Ca", // 18 = AD_TYPE_Ca
		"Mn", // 19 = AD_TYPE_Mn
		"Cu", // 20 = AD_TYPE_Cu
		"Na", // 21 = AD_TYPE_Na
		"K" , // 22 = AD_TYPE_K
		"Hg", // 23 = AD_TYPE_Hg
		"Ni", // 24 = AD_TYPE_Ni
		"Co", // 25 = AD_TYPE_Co
		"Cd", // 26 = AD_TYPE_Cd
		"As", // 27 = AD_TYPE_As
	};

	// http://en.wikipedia.org/wiki/Atomic_radii_of_the_elements_(data_page)
	// http://en.wikipedia.org/wiki/Covalent_radius
	// The above two references have inconsistent values for covalent radius.
	// The following definitions use the first reference, while OpenBabel uses the second.
	const fl ad_covalent_radii[] = ///< AutoDock4 covalent radii, factorized by 1.1 for extra allowance.
	{
		0.407, //  0 = AD_TYPE_HD, 0.407 = 1.1 * 0.37
		0.407, //  1 = AD_TYPE_H , 0.407 = 1.1 * 0.37
		0.847, //  2 = AD_TYPE_C , 0.847 = 1.1 * 0.77
		0.847, //  3 = AD_TYPE_A , 0.847 = 1.1 * 0.77
		0.825, //  4 = AD_TYPE_N , 0.825 = 1.1 * 0.75
		0.825, //  5 = AD_TYPE_NA, 0.825 = 1.1 * 0.75
		0.803, //  6 = AD_TYPE_OA, 0.803 = 1.1 * 0.73
		1.122, //  7 = AD_TYPE_SA, 1.122 = 1.1 * 1.02
		1.122, //  8 = AD_TYPE_S , 1.122 = 1.1 * 1.02
		1.276, //  9 = AD_TYPE_Se, 1.276 = 1.1 * 1.16
		1.166, // 10 = AD_TYPE_P , 1.166 = 1.1 * 1.06
		0.781, // 11 = AD_TYPE_F , 0.781 = 1.1 * 0.71
		1.089, // 12 = AD_TYPE_Cl, 1.089 = 1.1 * 0.99
		1.254, // 13 = AD_TYPE_Br, 1.254 = 1.1 * 1.14
		1.463, // 14 = AD_TYPE_I , 1.463 = 1.1 * 1.33
		1.441, // 15 = AD_TYPE_Zn, 1.441 = 1.1 * 1.31
		1.375, // 16 = AD_TYPE_Fe, 1.375 = 1.1 * 1.25
		1.430, // 17 = AD_TYPE_Mg, 1.430 = 1.1 * 1.30
		1.914, // 18 = AD_TYPE_Ca, 1.914 = 1.1 * 1.74
		1.529, // 19 = AD_TYPE_Mn, 1.529 = 1.1 * 1.39
		1.518, // 20 = AD_TYPE_Cu, 1.518 = 1.1 * 1.38
		1.694, // 21 = AD_TYPE_Na, 1.694 = 1.1 * 1.54
		2.156, // 22 = AD_TYPE_K , 2.156 = 1.1 * 1.96
		1.639, // 23 = AD_TYPE_Hg, 1.639 = 1.1 * 1.49
		1.331, // 24 = AD_TYPE_Ni, 1.331 = 1.1 * 1.21
		1.386, // 25 = AD_TYPE_Co, 1.386 = 1.1 * 1.26
		1.628, // 26 = AD_TYPE_Cd, 1.628 = 1.1 * 1.48
		1.309  // 27 = AD_TYPE_As, 1.309 = 1.1 * 1.19
	};

	const fl ad_atomic_weights[] = ///< AutoDock4 atomic weights.
	{
		  1.008,//  0 = AD_TYPE_HD
		  1.008,//  1 = AD_TYPE_H
		 12.01, //  2 = AD_TYPE_C
		 12.01, //  3 = AD_TYPE_A
		 14.01, //  4 = AD_TYPE_N
		 14.01, //  5 = AD_TYPE_NA
		 16.00, //  6 = AD_TYPE_OA
		 32.07, //  7 = AD_TYPE_SA
		 32.07, //  8 = AD_TYPE_S
		 78.96, //  9 = AD_TYPE_Se
		 30.97, // 10 = AD_TYPE_P
		 19.00, // 11 = AD_TYPE_F
		 35.45, // 12 = AD_TYPE_Cl
		 79.90, // 13 = AD_TYPE_Br
		126.90, // 14 = AD_TYPE_I
		 65.39, // 15 = AD_TYPE_Zn
		 55.84, // 16 = AD_TYPE_Fe
		 24.31, // 17 = AD_TYPE_Mg
		 40.08, // 18 = AD_TYPE_Ca
		 54.94, // 19 = AD_TYPE_Mn
		 63.55, // 20 = AD_TYPE_Cu
		 22.99, // 21 = AD_TYPE_Na
		 39.10, // 22 = AD_TYPE_K
		200.59, // 23 = AD_TYPE_Hg
		 58.69, // 24 = AD_TYPE_Ni
		 58.93, // 25 = AD_TYPE_Co
		112.41, // 26 = AD_TYPE_Cd
		 74.92  // 27 = AD_TYPE_As
	};

	/// Parses right-justified 1-based [i, j] of str into generic type T lexically.
	/// This conversion does not apply to left-justified values.
	template<typename T>
	inline T right_cast(const string& str, const size_t i, const size_t j)
	{
		const size_t start = str.find_first_not_of(' ', i - 1);
		return boost::lexical_cast<T > (str.substr(start, j - start));
	}

	// Represents an atom.
	class atom
	{
	public:
		string columns_13_to_30; ///< Columns from 1-based [13, 30] of an ATOM/HETATM line in pdbqt format.
		string columns_55_to_79; ///< Columns from 1-based [55, 79] of an ATOM/HETATM line in pdbqt format.
		size_t number; ///< Serial number.		
		vec3 coordinate; ///< 3D coordinate.
		size_t ad; ///< AutoDock4 atom type.

		/// Parses AutoDock4 atom type name, and returns AD_TYPE_SIZE if it does not match any supported AutoDock4 atom types.
		static size_t parse_ad_name(const string& ad_name)
		{
			for (size_t i = 0; i < AD_TYPE_SIZE; ++i)
				if (ad_names[i] == ad_name) return i;
			return AD_TYPE_SIZE;
		}

		/// Constructs an atom from an ATOM/HETATM line in pdbqt format.
		explicit atom(const string& line) : columns_13_to_30(line.substr(12, 18)), columns_55_to_79(line.substr(54)), number(right_cast<size_t>(line, 7, 11)), coordinate(vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54))), ad(parse_ad_name(line.substr(77, isspace(line[78]) ? 1 : 2))) {}

		/// Copy constructor.
		atom(const atom& a) : columns_13_to_30(a.columns_13_to_30), columns_55_to_79(a.columns_55_to_79), number(a.number), coordinate(a.coordinate), ad(a.ad) {}

		/// Move constructor.
		atom(atom&& a) : columns_13_to_30(static_cast<string&&>(a.columns_13_to_30)), columns_55_to_79(static_cast<string&&>(a.columns_55_to_79)), number(a.number), coordinate(a.coordinate), ad(a.ad) {}

		/// Returns covalent radius from an AutoDock4 atom type.
		const fl covalent_radius() const
		{
			return ad_covalent_radii[ad];
		}

		/// Returns atomic weight from an AutoDock4 atom type.
		const fl atomic_weight() const
		{
			return ad_atomic_weights[ad];
		}

		/// Returns true if the current atom is a hydrogen.
		const bool is_hydrogen() const
		{
			return (ad <= AD_TYPE_H);
		}

		/// Returns true if the current atom is a halogen.
		const bool is_halogen() const
		{
			return ((ad <= AD_TYPE_F) && (ad <= AD_TYPE_I));
		}

		/// Returns true if the current atom is a mutation point.
		const bool is_mutation_point() const
		{
			return (is_hydrogen() || is_halogen());
		}

		/// Returns true is the current atom is a hydrogen bond donor.
		const bool is_hb_donor() const
		{
			return ((!ad) || ((AD_TYPE_Zn <= ad) && (ad <= AD_TYPE_As)));
		}

		/// Returns true is the current atom is a hydrogen bond acceptor.
		const bool is_hb_acceptor() const
		{
			return ((AD_TYPE_NA <= ad) && (ad <= AD_TYPE_SA));
		}

		const bool is_neighbor(const atom& a) const
		{
			return (distance_sqr(coordinate, a.coordinate) < sqr(covalent_radius() + a.covalent_radius()));
		}
/*
		// planarize surroundings when this atom is found in this orbital
		bool isSP2() const
		{
			return ((element == "C") && (IndexArray.size() == 3));
		}

		// the Euclidean distance to another atom
		double DistanceTo(const atom& other) const
		{
			return (coordinates - other.coordinates).length();
		}
*/
	};

	/// Represents the index to a hydrogen or a halogen together with the index to its neighbor heavy atom.
	class mutation_point
	{
	public:
		size_t point; ///< The index to a mutation point, e.g. hydrogen or halogen.
		size_t neighbor; ///< The index to the neighbor of the current mutation point.

		explicit mutation_point(const size_t point, const size_t neighbor) : point(point), neighbor(neighbor) {}
	};
}

#endif
