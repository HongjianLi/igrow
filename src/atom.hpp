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
	const fl ad_covalent_radii[] = ///< AutoDock4 covalent radii.
	{
		0.37, //  0 = AD_TYPE_HD
		0.37, //  1 = AD_TYPE_H
		0.77, //  2 = AD_TYPE_C
		0.77, //  3 = AD_TYPE_A
		0.75, //  4 = AD_TYPE_N
		0.75, //  5 = AD_TYPE_NA
		0.73, //  6 = AD_TYPE_OA
		1.02, //  7 = AD_TYPE_SA
		1.02, //  8 = AD_TYPE_S
		1.16, //  9 = AD_TYPE_Se
		1.06, // 10 = AD_TYPE_P
		0.71, // 11 = AD_TYPE_F
		0.99, // 12 = AD_TYPE_Cl
		1.14, // 13 = AD_TYPE_Br
		1.33  // 14 = AD_TYPE_I
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

	/// Parses AutoDock4 atom type name, and returns AD_TYPE_SIZE if it does not match any supported AutoDock4 atom types.
	inline size_t parse_ad_type_string(const string& ad_type_string)
	{
		for (size_t i = 0; i < AD_TYPE_SIZE; ++i)
		{
			if (ad_type_strings[i] == ad_type_string) return i;
		}
		return AD_TYPE_SIZE;
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

		/// Constructs an atoms.
		explicit atom(const string& columns_13_to_30, const string columns_55_to_79, const size_t number, const vec3& coordinate, const size_t ad) : columns_13_to_30(columns_13_to_30), columns_55_to_79(columns_55_to_79), number(number), coordinate(coordinate), ad(ad) {}

		/// Copy constructor.
		atom(const atom& a) : columns_13_to_30(a.columns_13_to_30), columns_55_to_79(a.columns_55_to_79), number(a.number), coordinate(a.coordinate), ad(a.ad) {}

		/// Move constructor.
		atom(atom&& a) : columns_13_to_30(static_cast<string&&>(a.columns_13_to_30)), columns_55_to_79(static_cast<string&&>(a.columns_55_to_79)), number(a.number), coordinate(a.coordinate), ad(a.ad) {}
		
		/// Copy assignment.
		atom& operator=(const atom& a)
		{
			if (this != &a)
			{
				columns_13_to_30 = a.columns_13_to_30;
				columns_55_to_79 = a.columns_55_to_79;
				number = a.number;
				coordinate = a.coordinate;
				ad = a.ad;
			}
			return *this;
		}
		
		/// Move assignment.
		atom& operator=(atom&& a)
		{
			if (this != &a)
			{
				columns_13_to_30 = static_cast<string&&>(a.columns_13_to_30);
				columns_55_to_79 = static_cast<string&&>(a.columns_55_to_79);
				number = a.number;
				coordinate = a.coordinate;
				ad = a.ad;
			}
			return *this;
		}
		
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
			return ((AD_TYPE_F <= ad) && (ad <= AD_TYPE_I));
		}

		/// Returns true if the current atom is a mutable atom.
		const bool is_mutable() const
		{
			return (is_hydrogen() || is_halogen());
		}

		/// Returns true is the current atom is a hydrogen bond donor, i.e. polar hydrogen.
		const bool is_hb_donor() const
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
			BOOST_ASSERT(this != &a);
			return (distance_sqr(coordinate, a.coordinate) < sqr(1.1 * (covalent_radius() + a.covalent_radius())));
		}
	};
}

#endif
