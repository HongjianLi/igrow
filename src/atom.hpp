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
	const size_t AD_TYPE_SIZE = 15; ///< Number of supported AutoDock4 atom types.

	const string ad_type_strings[] = ///< AutoDock4 atom type names.
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
		"I"   // 14 = AD_TYPE_I
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
		1.463  // 14 = AD_TYPE_I , 1.463 = 1.1 * 1.33
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
		126.90  // 14 = AD_TYPE_I
	};

	/// Parses right-justified 1-based [i, j] of str into generic type T lexically.
	/// This conversion does not apply to left-justified values.
	template<typename T>
	inline T right_cast(const string& str, const size_t i, const size_t j)
	{
		const size_t start = str.find_first_not_of(' ', i - 1);
		return boost::lexical_cast<T > (str.substr(start, j - start));
	}

	/// Represents the index to an atom.
	class atom_index
	{
	public:
		size_t frame; ///< The index to the frame to which the atom belongs.
		size_t index; ///< The index to the atom within its frame.

		explicit atom_index(const size_t frame, const size_t index) : frame(frame), index(index) {}
	};
	
	// Represents an atom.
	class atom
	{
	public:
		string columns_13_to_30; ///< Columns from 1-based [13, 30] of an ATOM/HETATM line in pdbqt format.
		string columns_55_to_79; ///< Columns from 1-based [55, 79] of an ATOM/HETATM line in pdbqt format.
		size_t number; ///< Serial number.
		vec3 coordinate; ///< 3D coordinate.
		size_t ad; ///< AutoDock4 atom type.
		vector<atom_index> neighbors;

		/// Parses AutoDock4 atom type name, and returns AD_TYPE_SIZE if it does not match any supported AutoDock4 atom types.
		static size_t parse_ad_type_string(const string& ad_type_string)
		{
			for (size_t i = 0; i < AD_TYPE_SIZE; ++i)
				if (ad_type_strings[i] == ad_type_string) return i;
			return AD_TYPE_SIZE;
		}

		/// Constructs an atom from an ATOM/HETATM line in pdbqt format.
		explicit atom(const string& line) : columns_13_to_30(line.substr(12, 18)), columns_55_to_79(line.substr(54)), number(right_cast<size_t>(line, 7, 11)), coordinate(vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54))), ad(parse_ad_type_string(line.substr(77, isspace(line[78]) ? 1 : 2)))
		{
			neighbors.reserve(4); // An atom typically consists of <= 4 neighbors.
		}

		/// Copy constructor.
		atom(const atom& a) : columns_13_to_30(a.columns_13_to_30), columns_55_to_79(a.columns_55_to_79), number(a.number), coordinate(a.coordinate), ad(a.ad), neighbors(a.neighbors) {}

		/// Move constructor.
		atom(atom&& a) : columns_13_to_30(static_cast<string&&>(a.columns_13_to_30)), columns_55_to_79(static_cast<string&&>(a.columns_55_to_79)), number(a.number), coordinate(a.coordinate), ad(a.ad), neighbors(static_cast<vector<atom_index>&&>(a.neighbors)) {}

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

		/// Returns true if the current atom is a mutable atom.
		const bool is_mutable() const
		{
			return (is_hydrogen() || is_halogen());
		}

		/// Returns true is the current atom is a polar hydrogen.
		const bool is_polar_hydrogen() const
		{
			return (!ad);
		}

		/// Returns true is the current atom is a hydrogen bond acceptor.
		const bool is_hb_acceptor() const
		{
			return ((AD_TYPE_NA <= ad) && (ad <= AD_TYPE_SA));
		}

		/// Returns true if the current atom is covalently bonded to a given atom.
		const bool is_neighbor(const atom& a) const
		{
			return (distance_sqr(coordinate, a.coordinate) < sqr(covalent_radius() + a.covalent_radius()));
		}

		/// Returns true if the current atom is sp2 hybridized.
		const bool is_sp2() const
		{
			return (((ad == AD_TYPE_C) || (ad == AD_TYPE_A)) && (neighbors.size() == 3));
		}
	};
}

#endif
