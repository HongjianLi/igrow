/*

   Copystd::right (c) 2011, The Chinese University of Hong Kong

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

#include <sstream>
#include "atom.hpp"

namespace igrow
{
	bool atom::isSP2() const
	{
		return ((element == "C") && (IndexArray.size() == 3));
	}

	// Eucledian distance between 2 atoms

	double atom::DistanceTo(const atom& other) const
	{
		return (coordinates - other.coordinates).length();
	}

	void atom::ReadPDBLine(const string& line)
	{
		// Parse atom serial number.
		const size_t start = line.find_first_not_of(' ', 6);
		PDBIndex = line.substr(start, 11 - start);

		// Parse atom name.
		name = line.substr(13, 3);

		// Parse residue name.
		Residue = line.substr(16, 4);
		if (Residue == "    ") Residue = " MOL";

		// Parse X,Y,Z coordinates.
		coordinates.n[0] = right_cast<int>(line, 31, 38);
		coordinates.n[1] = right_cast<int>(line, 39, 46);
		coordinates.n[2] = right_cast<int>(line, 47, 54);

		// Parse element symbol.
		if (line[76] != ' ') element = toupper(line[76]);
		element += toupper(line[77]);
	}
}
