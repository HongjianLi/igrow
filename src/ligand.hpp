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

#ifndef IGROW_LIGAND_HPP
#define IGROW_LIGAND_HPP

#include <map>
#include <list>
#include <set>
#include <boost/filesystem/path.hpp>
#include "atom.hpp"

namespace igrow
{
	using std::set;
	using std::map;
	using boost::filesystem::path;

	// Represents a ligand.

	class ligand
	{
	public:
		static const fl pi;

		// a collection of atoms indicated by their indice
		map<int, atom> atoms;
		int ID;
		// load information of a PDB file in the form of a molecule
		void load(const path& file);
		// produce a PDB file based on this molecule
		void save(const path& file);
		// break the molecule into two while retaining the atoms indicated by the reference
		ligand split(const ligand& ref);
		// add a fragment to the molecule by replacing a hydrogen in the original structure
		void mutate(ligand fragment);
		// check whether all the atoms are no less than a certain threshold distance
		bool hasBadBonds();
		// not implemented
		bool hasStericClashes();
		// the largest index of atoms in this molecule
		int MaxIndex();
		// the smallest index of atoms in this molecule
		int MinIndex();
		// remove atom of index and its connection, return an iterator for further process
		map<int, atom>::reverse_iterator DeleteAtom(int index);
		// a distance between two molecules given by the sum of minimum distance of each atom
		double MolecularDistance(ligand& other);
		// measuring distance between conformation, supported by RMSD
		double IntraMolecularDistance(ligand& other);
		// obtain the index of one hydrogen of this molecule
		int IndexOfRandomHydrogen();
		// move atom of index to origin along with its connected atoms
		void Translate(int index, Vec3d origin);
		// translation based on centre
		void Translate(Vec3d origin);
		// rotate along the line given by v1-v2 using indexed atom as pivot in terms of radian
		void RotateLine(Vec3d v1, Vec3d v2, int index, double radian);
		// rotate along the line normal using indexed atom as pivot in terms of radian
		void RotateLine(Vec3d normal, int index, double radian);
		// rotate along the line passing through CG of the molecule
		void RotateLine(Vec3d normal, double radian);
		// rotate using Euler angles with indexed atom as centre
		void Rotate(int index, Vec3d EulerAngles);
		// rotation based on centre
		void Rotate(Vec3d EulerAngles);
		// add missing hydrogen based on some chemical rules (distance and orientation) when the PDB file is not complete
		void AddHydrogen();
		// method detect whether there is an aromatic ring in the structure, return number of rings and the first ring in argument
		int DetectAromatic(std::list<int>& ring, int index = -1);
		// detect whether there is a ring structure in the molecule
		int DetectRing(std::list<int>& ring, int index = -1);
		// follow more strictly with the chemical rules, the resulting ligand would be more likely feasible
		int synthesis(const path& file);
		// returns true if the ligand is valid.
		bool valid();
		// returns the molecular weight
		double mwt();

	protected:
		// set has no duplication
		// record the indice of atoms visited
		set<int> scanned;
		// record the indice of atoms which are the same in terms of type and position
		set<int> overlap;
		// a set of hydrogen to replace a removed fragment
		map<int, atom> toAddAtoms;

		// the dihedral angle among the four atoms
		double DihedralAngle(atom a1, atom a2, atom a3, atom a4);
		// the angle between v1-v2-v3
		double Angle(Vec3d v1, Vec3d v2, Vec3d v3);
		// return the coordinates of the molecule centre, taking the weight of each atom equally
		// todo: improve it using molecular weight
		Vec3d CentreOfGravity();
		// method to tranverse an atom, decide which part to retain
		void scan_recursive(int index);
		// private method to handle Euler angles based rotation
		inline void Rotate(Vec3d centre, Vec3d EulerAngles);
		// method recursively traverse element within possible ring
		void traverse_ring(std::list<int>& possible_ring, std::map<int, atom>& possible_list, int index);
		// method to split a list of connected index into multiple ring
		int split_ring(std::list<int>& possible_ring, std::list<int>& vector_of_rings, int index = -1);
		// merge the list based on connection information
		void merge_list(std::set<int>& connect, std::vector<std::list<int> >& rings, std::list<int>& input, int limit);
		// splitting fragment in 2 smaller parts
		int replace_bond(std::pair<int, int> bond);
		// scan the fragment from the start index and stop on the end index and return the fragment to target
		void split_fragment(ligand& target, int start, int end);
		// internal method to check number of heteroatom in a ring
		inline bool checkHeteroAtom(std::list<int>& candidate);
		// determine whether both molecules consist ring structure that can be joined
		int JoinRing(ligand& ref);
	};

}

#endif
