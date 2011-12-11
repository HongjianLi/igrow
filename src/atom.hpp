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

#include "vec.hpp"
#include "common.hpp"
#include <set>

class atom {
public:
	// 3D coordinates of an atom indicated by the PDB file
	Vec3d coordinates;
	// the type of atom as in a periodic table
	std::string element;
	// the name of an atom indicated by the atom type and a number, following the standard of IUPAC
	std::string name;
	// the index of connected atoms, given by PDB index
	std::set<int> IndexArray;
	// the index read from the PDB file
	std::string PDBIndex;
	// the residue this atom belongs tp
	std::string Residue;
	// a unique ID representing the ligand which could help distinguish in some operations
	int ID;

	atom();
	virtual ~atom();

	// planarize surroundings when this atom is found in this orbital
	bool isSP2();
	// the Euclidean distance to another atom
	double DistanceTo(atom &other);
	// parser to read a PDB line stated ATOM or HETATM
	void ReadPDBLine(std::string Line);
	// write the information of this atom in the PDB format
	std::string WritePDBLine(int index);
};

#endif
