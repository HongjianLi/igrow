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

#ifndef IGROW_ATOM_HPP
#define IGROW_ATOM_HPP

#include <set>
#include "common.hpp"
#include "vec.hpp"

namespace igrow
{
	using std::set;

	// Represents an atom.

	class atom
	{
	public:
		// 3D coordinates
		Vec3d coordinates;

		// the type of atom as in a periodic table
		string element;

		// the name of an atom indicated by the atom type and a number, following the standard of IUPAC
		string name;

		// the index of connected atoms, given by PDB index
		set<int> IndexArray;

		// the index read from the PDB file
		string PDBIndex;

		// the residue this atom belongs tO
		string Residue;

		// a unique ID representing the ligand which could help distinguish in some operations
		int ID;

		// planarize surroundings when this atom is found in this orbital
		bool isSP2() const;

		// the Euclidean distance to another atom
		double DistanceTo(const atom& other) const;

		// parser to read a PDB line stated ATOM or HETATM
		void ReadPDBLine(const string& line);

		// write the information of this atom in the PDB format
		string WritePDBLine(int index);
	};

}

#endif
