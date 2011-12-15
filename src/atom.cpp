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
#include <iomanip>
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
		string tempString, space(" ");
		int start, end;
		// name of atom is at column 11 for at most 5 characters
		tempString = line.substr(11, 5);
		start = tempString.find_first_not_of(space);
		end = tempString.find_last_not_of(space);
		name = tempString.substr(start, end - start + 1);
		// append space to the back if not enough
		if (name.size() == 1)
			name += space;
		if (name.size() == 2)
			name += space;
		// coordinates are recorded as Cartesian space in (x,y,z) in decimal format
		coordinates.n[0] = atof(line.substr(30, 8).c_str());
		coordinates.n[1] = atof(line.substr(38, 8).c_str());
		coordinates.n[2] = atof(line.substr(46, 8).c_str());
		// element type is also available in this PDB file
		if (line.length() >= 79)
		{
			if (isalpha(line[76]))
				element += toupper(line[76]);
			if (isalpha(line[77]))
				element += toupper(line[77]);
			if (isalpha(line[78]))
				element += toupper(line[78]);
		}
		// if element is not found, recover it in the name
		if (element.empty())
		{
			start = 0;
			// retrieve all alphabets in front of number
			while (start != name.length())
			{
				if (isalpha(name[start]))
					element += name[start];
				++start;
			}
		}
		for (int i = 0; i != name.length(); ++i)
			name[i] = toupper(name[i]);
		// read index of the atom
		tempString = line.substr(6, 6);
		start = tempString.find_first_not_of(space);
		end = tempString.find_last_not_of(space);
		PDBIndex = tempString.substr(start, end - start + 1);
		// read residue name, if none is found, assign it
		Residue = line.substr(16, 4);
		if (Residue.find_first_not_of(space) == string::npos)
			Residue = " MOL";
	}

	string atom::WritePDBLine(int index)
	{
		string line("ATOM ");
		// use stream to collect numerous information
		std::ostringstream buffer;
		buffer.str(string());
		// fill space and append index to std::right, taking up 6 slots
		buffer << std::setfill(' ') << std::setw(6);
		buffer << std::right << index;
		line += buffer.str();
		buffer.str(string());
		// append name to std::right, taking up 5 slots
		buffer << std::setfill(' ') << std::setw(5);
		buffer << std::right << name;
		line += buffer.str();
		line += Residue;
		buffer.str(string());
		// set float precision to 2 decimal places, x-coordinate takes up 18 slots
		buffer << std::setfill(' ') << std::setw(18);
		// reduce precision when the coordinate is too small
		buffer.precision(3);
		// prevent overflow on the entry
		/*if (abs(coordinates.n[0]) < 0.01)
			buffer.precision(2);
		else
			buffer.precision(3);
		if (abs(coordinates.n[0]) < 0.001) coordinates.n[0] = 0;*/
		buffer << std::fixed << std::right << coordinates.n[0];
		line += buffer.str();
		buffer.str(string());
		// y-coordinate takes up 8 slots
		buffer << std::setfill(' ') << std::setw(8);
		/*if (abs(coordinates.n[1]) < 0.01)
			buffer.precision(2);
		else
			buffer.precision(3);
		if (abs(coordinates.n[1]) < 0.001) coordinates.n[1] = 0;*/
		buffer << std::fixed << std::right << coordinates.n[1];
		line += buffer.str();
		buffer.str(string());
		// z-coordinate takes up 8 slots
		buffer << std::setfill(' ') << std::setw(8);
		/*if (abs(coordinates.n[2]) < 0.01)
			buffer.precision(2);
		else
			buffer.precision(3);
		if (abs(coordinates.n[2]) < 0.001) coordinates.n[2] = 0;*/
		buffer << std::fixed << std::right << coordinates.n[2];
		line += buffer.str();
		buffer.str(string());
		// append element type that takes up 24 slots
		buffer << std::setfill(' ') << std::setw(24);
		buffer << std::right << element;
		if (element.length() == 1)
			buffer << std::right << "  ";
		else
			buffer << std::right << " ";
		line += buffer.str();
		return line;
	}

}
