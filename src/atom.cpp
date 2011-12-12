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

#include "atom.hpp"
#include <sstream>
#include <iomanip>
#include <stdio.h>

using namespace std;

bool atom::isSP2() {
	bool flag(false);
	string carbon("C");
	// the element is carbon and has 3 connections
	if (!element.compare(carbon) && (IndexArray.size() == 3))
		flag = true;
	// nitrogen could also produce SP2 orbital, needs improvement
	return flag;
}

// Eucledian distance between 2 atoms
double atom::DistanceTo(atom &other) {
	return (coordinates-other.coordinates).length();
}

void atom::ReadPDBLine(string Line) {
	string tempString, space(" ");
	int start, end;
	// name of atom is at column 11 for at most 5 characters
	tempString = Line.substr(11,5);
	start = tempString.find_first_not_of(space);
	end = tempString.find_last_not_of(space);
	name = tempString.substr(start, end-start+1);
	// append space to the back if not enough
	if (name.size() == 1)
		name += space;
	if (name.size() == 2)
		name += space;
	// coordinates are recorded as Cartesian space in (x,y,z) in decimal format
	coordinates.n[0] = atof(Line.substr(30,8).c_str());
	coordinates.n[1] = atof(Line.substr(38,8).c_str());
	coordinates.n[2] = atof(Line.substr(46,8).c_str());
	// element type is also available in this PDB file
	if (Line.length() >= 79) {
		if (isalpha(Line[76]))
			element += toupper(Line[76]);
		if (isalpha(Line[77]))
			element += toupper(Line[77]);
		if (isalpha(Line[78]))
			element += toupper(Line[78]);
	}
	// if element is not found, recover it in the name
	if (element.empty()) {
		start = 0;
		// retrieve all alphabets in front of number
		while (start != name.length()) {
			if (isalpha(name[start]))
				element += name[start];
			++start;
		}
	}
	for (int i = 0; i != name.length(); ++i)
		name[i] = toupper(name[i]);
	// read index of the atom
	tempString = Line.substr(6,6);
	start = tempString.find_first_not_of(space);
	end = tempString.find_last_not_of(space);
	PDBIndex = tempString.substr(start,end-start+1);
	// read residue name, if none is found, assign it
	Residue = Line.substr(16,4);
	if (Residue.find_first_not_of(space) == string::npos)
		Residue = " MOL";
}

string atom::WritePDBLine(int index) {
	string Line("ATOM ");
	// use stream to collect numerous information
	ostringstream buffer;
	buffer.str(string());
	// fill space and append index to right, taking up 6 slots
	buffer << setfill(' ') << setw(6);
	buffer << right << index;
	Line += buffer.str();
	buffer.str(string());
	// append name to right, taking up 5 slots
	buffer << setfill(' ') << setw(5);
	buffer << right << name;
	Line += buffer.str();
	Line += Residue;
	buffer.str(string());
	// set float precision to 2 decimal places, x-coordinate takes up 18 slots
	buffer << setfill(' ') << setw(18);
	// reduce precision when the coordinate is too small
	buffer.precision(3);
	// prevent overflow on the entry
	/*if (abs(coordinates.n[0]) < 0.01)
		buffer.precision(2);
	else
		buffer.precision(3);
	if (abs(coordinates.n[0]) < 0.001) coordinates.n[0] = 0;*/
	buffer << fixed << right << coordinates.n[0];
	Line += buffer.str();
	buffer.str(string());
	// y-coordinate takes up 8 slots
	buffer << setfill(' ') << setw(8);
	/*if (abs(coordinates.n[1]) < 0.01)
		buffer.precision(2);
	else
		buffer.precision(3);
	if (abs(coordinates.n[1]) < 0.001) coordinates.n[1] = 0;*/
	buffer << fixed << right << coordinates.n[1];
	Line += buffer.str();
	buffer.str(string());
	// z-coordinate takes up 8 slots
	buffer << setfill(' ') << setw(8);
	/*if (abs(coordinates.n[2]) < 0.01)
		buffer.precision(2);
	else
		buffer.precision(3);
	if (abs(coordinates.n[2]) < 0.001) coordinates.n[2] = 0;*/
	buffer << fixed << right << coordinates.n[2];
	Line += buffer.str();
	buffer.str(string());
	// append element type that takes up 24 slots
	buffer << setfill(' ') << setw(24);
	buffer << right << element;
	if (element.length() == 1)
		buffer << right << "  ";
	else
		buffer << right << " ";
	Line += buffer.str();
	return Line;
}
