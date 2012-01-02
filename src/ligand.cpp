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

#include <iomanip>
#include <boost/random.hpp>
#include "fstream.hpp"
#include "ligand.hpp"

namespace igrow
{
	/// Represents a parsing error.
	class parsing_error : public std::domain_error
	{
	public:
		/// Constructs a parsing error.
		parsing_error(const path& file, const size_t line, const string& reason) : std::domain_error("Error parsing \"" + file.filename().string() + "\" on line " + boost::lexical_cast<string>(line) + ": " + reason) {}
	};

	/// Returns true if a string starts with another string.
	bool starts_with(const string& str, const string& start)
	{
		const size_t start_size = start.size();
		if (str.size() < start_size) return false;
		for (size_t i = 0; i < start_size; ++i)
		{
			if (str[i] != start[i]) return false;
		}
		return true;
	}

	ligand::ligand(const path& p) : p(p), num_atoms(0), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0), connector1(0), connector2(0), logp(0) // TODO: comment logp(0)
	{
		// Initialize necessary variables for constructing a ligand.
		frames.reserve(30); // A ligand typically consists of <= 30 frames.
		frames.push_back(frame(0)); // ROOT is also treated as a frame. The parent of ROOT frame is dummy.
		mutable_atoms.reserve(20); // A ligand typically consists of <= 20 mutation points.

		// Initialize helper variables for parsing.
		size_t current = 0; // Index of current frame, initialized to ROOT frame.
		size_t num_lines = 0; // Used to track line number for reporting parsing errors, if any.
		string line; // A line of ligand file in pdbqt format.
		line.reserve(79); // According to PDBQT specification, the last item AutoDock atom type locates at 1-based [78, 79].

		// Parse ATOM/HETATM, BRANCH, ENDBRANCH.
		ifstream in(p); // Parsing starts. Open the file stream as late as possible.
		while (getline(in, line))
		{
			++num_lines;
			if (starts_with(line, "ATOM") || starts_with(line, "HETATM"))
			{
				// Parse the ATOM/HETATM line into an atom, which belongs to the last frame.
				frame& f = frames.back();
				f.atoms.push_back(atom(line));				
				
				// Validate the AutoDock4 atom type.
				atom& a = f.atoms.back();
				if (a.ad == AD_TYPE_SIZE) throw parsing_error(p, num_lines, "Atom type " + line.substr(77, isspace(line[78]) ? 1 : 2) + " is not supported by igrow.");
				
				// Find the neighbors in the same frame.
				const atom_index idx = atom_index(frames.size() - 1, f.atoms.size() - 1); // The index to the current atom.
				for (size_t i = 0; i < idx.index; ++i) // Exclude the last atom which is a itself.
				{
					atom& b = f.atoms[i];
					if (a.is_neighbor(b))
					{
						a.neighbors.push_back(atom_index(idx.frame, i));
						BOOST_ASSERT(a.neighbors.size() <= 4);
						b.neighbors.push_back(idx);
						BOOST_ASSERT(b.neighbors.size() <= 4);
					}
				}
				
				// Add rotorX as a neighbor if the current atom is rotorY of a BRANCH frame.
				if ((current) && (!idx.index))
				{
					BOOST_ASSERT(f.atoms.size() == 1);
					a.neighbors.push_back(atom_index(f.parent, f.rotorX));
					BOOST_ASSERT(a.neighbors.size() <= 4);
					frames[f.parent].atoms[f.rotorX].neighbors.push_back(idx);
					BOOST_ASSERT(frames[f.parent].atoms[f.rotorX].neighbors.size() <= 4);
				}
				
				++num_atoms;
				if (a.is_mutable()) mutable_atoms.push_back(idx);
				if (!a.is_hydrogen()) ++num_heavy_atoms;
				if (a.is_hb_donor()) ++num_hb_donors;
				if (a.is_hb_acceptor()) ++num_hb_acceptors;
				mw += a.atomic_weight();
			}
			else if (starts_with(line, "BRANCH"))
			{
				// Insert a new frame whose parent is the current frame.
				frames.push_back(frame(current));

				// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
				// Y is not parsed because the atom immediately follows "BRANCH" must be Y in pdbqt files created by the prepare_ligand4.py script of MGLTools.
				// This assumption fails if pdbqt files are prepared by OpenBabel. In this case, the class frame should incorporate a new field rotorY to track the origin.
				const size_t x = right_cast<size_t>(line, 7, 10);

				// Find the corresponding heavy atom with X as its atom number in the current frame.
				frame& f = frames[current];
				for (size_t i = 0; i < f.atoms.size(); ++i)
				{
					if (f.atoms[i].number == x)
					{
						frames.back().rotorX = i;
						break;
					}
				}

				// Now the current frame is the newly inserted BRANCH frame.
				current = frames.size() - 1;

				// The parent frame has the current frame as one of its branches.
				f.branches.push_back(current);
			}
			else if (starts_with(line, "ENDBRANCH"))
			{
				// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
				// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
				if (frames.back().atoms.empty()) throw parsing_error(p, num_lines, "An empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

				// Now the parent of the following frame is the parent of current frame.
				current = frames[current].parent;
			}
		}
		in.close(); // Parsing finishes. Close the file stream as soon as possible.

		BOOST_ASSERT(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
		BOOST_ASSERT(num_atoms >= num_heavy_atoms);
		BOOST_ASSERT(num_atoms + 3 <= num_lines); // ROOT, ENDROOT, TORSDOF
		
		// Determine the number of rotatable bonds.
		num_rotatable_bonds = frames.size() - 1;

		// Determine if the current ligand is able to perform mutation or crossover.
		mutation_feasible = !mutable_atoms.empty();
		crossover_feasible = frames.size() > 1;
	}

	void ligand::save(const path& p)
	{
		// Update the path to the current ligand.
		this->p = p;

		using namespace std;
		ofstream out(p); // Dumping starts. Open the file stream as late as possible.
		out.setf(ios::fixed, ios::floatfield);
		out << setprecision(3);

		// Dump the ROOT frame.
		out << "ROOT\n";
		{
			const frame& f = frames.front();
			const size_t num_atoms = f.atoms.size();
			for (size_t i = 0; i < num_atoms; ++i)
			{
				const atom& a = f.atoms[i];
				out << "ATOM  " << setw(5) << a.number << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << endl;
			}
		}
		out << "ENDROOT\n";

		// Dump the BRANCH frames.
		const size_t num_frames = frames.size();
		vector<bool> dump_branches(num_frames); // dump_branches[0] is dummy. The ROOT frame has been dumped.
		vector<size_t> stack; // Stack to track the traversal sequence of frames in order to avoid recursion.
		stack.reserve(num_frames - 1); // The ROOT frame is excluded.
		{
			const frame& f = frames.front();
			for (size_t i = f.branches.size(); i > 0;)
			{
				stack.push_back(f.branches[--i]);
			}
		}
		while (!stack.empty())
		{
			const size_t fn = stack.back();
			const frame& f = frames[fn];
			if (!dump_branches[fn]) // This BRANCH frame has not been dumped.
			{
				out << "BRANCH"    << setw(4) << frames[f.parent].atoms[f.rotorX].number << setw(4) << f.atoms.front().number << endl;
				const size_t num_atoms = f.atoms.size();
				for (size_t i = 0; i < num_atoms; ++i)
				{
					const atom& a = f.atoms[i];
					out << "ATOM  " << setw(5) << a.number << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << endl;
				}
				dump_branches[fn] = true;
				for (size_t i = f.branches.size(); i > 0;)
				{
					stack.push_back(f.branches[--i]);
				}
			}
			else // This BRANCH frame has been dumped.
			{
				out << "ENDBRANCH" << setw(4) << frames[f.parent].atoms[f.rotorX].number << setw(4) << f.atoms.front().number << endl;
				stack.pop_back();
			}
		}
		out << "TORSDOF " << num_rotatable_bonds << std::endl;
		out.close();
	}

	ligand* ligand::mutate(const ligand& other, const mt19937eng& eng) const
	{
		// Initialize random number generators for obtaining two random mutable atoms.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t> > uniform_mutable_atom_gen_1(eng, uniform_int_distribution<size_t>(0, this->mutable_atoms.size() - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t> > uniform_mutable_atom_gen_2(eng, uniform_int_distribution<size_t>(0, other.mutable_atoms.size() - 1));

		// Obtain a random mutable atom from the current ligand and the other ligand respectively.
		const atom_index& ma1 = this->mutable_atoms[uniform_mutable_atom_gen_1()];
		const atom_index& ma2 = other.mutable_atoms[uniform_mutable_atom_gen_2()];

		// Obtain constant references to the mutable atoms.
		const frame& f1 = this->frames[ma1.frame];
		const frame& f2 = other.frames[ma2.frame];
		const atom& p1 = f1.atoms[ma1.index]; // Constant reference to the mutable atom of the current ligand.
		const atom& p2 = f2.atoms[ma2.index]; // Constant reference to the mutable atom of the other ligand.

		// The mutable atom should consist of only one non-rotatable single bond.
		BOOST_ASSERT(p1.neighbors.size() == 1);
		BOOST_ASSERT(p2.neighbors.size() == 1);

		// Both the connector and mutable atom should be in the same frame.
		BOOST_ASSERT(p1.neighbors.front().frame == ma1.frame);
		BOOST_ASSERT(p2.neighbors.front().frame == ma2.frame);

		// Obtain constant references to the connector atoms.
		const atom& c1 = f1.atoms[p1.neighbors.front().index];
		const atom& c2 = f2.atoms[p2.neighbors.front().index];

		// Initialize a child ligand.
		ligand child;
		child.parent1 = this->p;
		child.parent2 = other.p;
		child.connector1 = c1.number;
		child.connector2 = c2.number;
		
		// The number of rotatable bonds of child ligand is equal to the sum of its parent ligands plus 1.
		child.num_rotatable_bonds = this->num_rotatable_bonds + other.num_rotatable_bonds + 1;
		
		// The number of atoms of child ligand is equal to the sum of its parent ligands minus 2.
		child.num_atoms = this->num_atoms + other.num_atoms - 2;

		// The number of heavy atoms of child ligand is equal to the sum of its parent ligands minus the two mutable atoms if they are heavy atoms.
		child.num_heavy_atoms = this->num_heavy_atoms + other.num_heavy_atoms - ((p1.is_hydrogen() ? 0 : 1) + (p2.is_hydrogen() ? 0 : 1));

		// The number of hydrogen bond donors of child ligand is equal to the sum of its parent ligands minus the two mutable atoms if they are hydrogen bond donors.
		child.num_hb_donors = this->num_hb_donors + other.num_hb_donors - ((p1.is_hb_donor() ? 1 : 0) + (p2.is_hb_donor() ? 1 : 0));

		// The number of hydrogen bond acceptors of child ligand is equal to the sum of its parent ligands.
		child.num_hb_acceptors = this->num_hb_acceptors + other.num_hb_acceptors;

		// The molecular weight of child ligand is equal to the sum of its parent ligands minus the two mutable atoms.
		child.mw = this->mw + other.mw - (p1.atomic_weight() + p2.atomic_weight());

		// The logP of child ligand is equal to the sum of its parent ligands minus the two mutable atoms.
		child.logp = this->logp + other.logp; // TODO: comment this line.

		// The number of frames of child ligand is equal to the sum of its parent ligands.
		const size_t num_frames_1 = this->frames.size();
		const size_t num_frames_2 = other.frames.size();
		child.frames.reserve(num_frames_1 + num_frames_2);
		
		// Find the number of current ligand's frames that will be directly copied to the child ligand.
		size_t num_frames_split = 0; // Frames [0, num_frames_split) constitute part 1, and frames [num_frames_split, num_frames_1) constitute part 2.
		size_t insertion_index; // The index at which the new frame number of the other ligand's ma2 frame will be inserted into the branches of the current ligand's ma1 frame.
		const size_t num_branches_1 = f1.branches.size();
		for (size_t i = 0; i < num_branches_1; ++i)
		{
			if (c1.number < f1.atoms[this->frames[f1.branches[i]].rotorX].number)
			{
				num_frames_split = f1.branches[i];
				BOOST_ASSERT(num_frames_split); // At least one frame of the current ligand will be directly copied.
				insertion_index = i;				
				break;
			}
		}
		if (!num_frames_split) // The insertion index is at the end of the branches.
		{
			insertion_index = f1.branches.size();
			
			// Trace back the parent frames to the immediately next branch number, which is equal to num_frames_split.
			size_t c = ma1.frame; // Current frame.
			size_t p = f1.parent; // Parent frame.
			while (c && (c == this->frames[p].branches.back()))
			{
				c = p;
				p = this->frames[p].parent;
			}
			if (c) // The next branch number has been found in frame p.
			{
				const frame& pf = this->frames[p];
				const size_t num_branches_p = pf.branches.size();
				for (size_t i = 0; i < num_branches_p; ++i)
				{
					if (pf.branches[i] == c)
					{
						BOOST_ASSERT(i + 1 < num_branches_p);
						num_frames_split = pf.branches[i + 1];
						break;
					}
				}
			}
			else // The entire frame tree of the current ligand has been traversed.
			{
				num_frames_split = num_frames_1;
			}
		}
		BOOST_ASSERT(num_frames_split >= ma1.frame + 1);
		BOOST_ASSERT(num_frames_split <= num_frames_1);

		// Copy part 1 of the current ligand to the child ligand.
		// TODO: use std::copy instead.
		for (size_t i = 0; i < num_frames_split; ++i)
		{
			child.frames.push_back(this->frames[i]);
		}
		
		// Insert the new frame number of ma2 frame of the other ligand at the previously found index.
		frame& f3 = child.frames[ma1.frame];
		f3.branches.insert(f3.branches.begin() + insertion_index, num_frames_split);
		
		// Increment the branch numbers of part 2 frames by num_frames_2.
		const size_t num_branches_3 = f3.branches.size();
		BOOST_ASSERT(num_branches_3 == num_branches_1 + 1);
		for (size_t i = insertion_index + 1; i < num_branches_3; ++i)
		{
			f3.branches[i] += num_frames_2;
		}
		for (size_t k = 0; k < ma1.frame; ++k)
		{
			frame& f = child.frames[k];
			const size_t num_branches_f = f.branches.size();
			for (size_t i = 0; i < num_branches_f; ++i)
			{
				if (f.branches[i] >= num_frames_split) f.branches[i] += num_frames_2;
			}
		}
		
		// Remove mutable atom 1 from ma1 frame.
		f3.atoms.erase(f3.atoms.begin() + ma1.index);
		
		// Update a connector neighbor of c3 from mutable atom 1 to the rotorY of ma2 frame of the other ligand.
		atom& c3 = f3.atoms[p1.neighbors.front().index];
		for (size_t i = c3.neighbors.size(); i > 0;)
		{
			if (c3.neighbors[--i] == ma1)
			{
				c3.neighbors[i] = atom_index(num_frames_split, 0); // Assume rotorY is always 0.
				break;
			}
		}		
		
		// Find the traversal sequence of the other ligand starting from ma2 frame, as well as its reverse mapping.
		vector<size_t> traversal;
		traversal.reserve(num_frames_2);
		vector<size_t> mapping(num_frames_2);
		vector<size_t> stack;
		stack.reserve(num_frames_2);
		stack.push_back(ma2.frame);
		while (!stack.empty())
		{
			const size_t k = stack.back();
			mapping[k] = traversal.size();
			traversal.push_back(k);			
			const frame& f = other.frames[k];
			stack.pop_back();
			for (size_t i = f.branches.size(); i > 0;)
			{
				const size_t j = f.branches[--i];
				if (std::find(traversal.cbegin(), traversal.cend(), j) == traversal.end()) stack.push_back(j);
			}
			if (std::find(traversal.cbegin(), traversal.cend(), f.parent) == traversal.end()) stack.push_back(f.parent);
		}		
		
		// Copy the other ligand to the child ligand, and update the frame numbers of parent, branches, and atom neighbors;
		for (size_t k = 0; k < num_frames_2; ++k)
		{
			child.frames.push_back(other.frames[traversal[k]]);
			frame& f = child.frames.back();
			f.parent = mapping[f.parent] + num_frames_split;
			const size_t num_branches = f.branches.size();
			for (size_t i = 0; i < num_branches; ++i)
			{
				f.branches[i] = mapping[f.branches[i]] + num_frames_split;
			}
			const size_t num_atoms = f.atoms.size();
			for (size_t i = 0; i < num_atoms; ++i)
			{
				atom& a = f.atoms[i];
				const size_t num_neighbors = a.neighbors.size();
				for (size_t j = 0; j < num_neighbors; ++j)
				{
					a.neighbors[j].frame = mapping[a.neighbors[j].frame] + num_frames_split;
				}
			}
		}
		
		// The parent path from ma1 in the other ligand becomes the branch path in the child ligand.
		if (ma2.frame)
		{
			const size_t k0 = mapping.front();
			BOOST_ASSERT(k0 > 0);
			const size_t least_branch = k0 - 1; // The least frame number of num_frames_split + k0 frame of the child ligand is always equal to k0 - 1.
			frame& f = child.frames[num_frames_split + k0];
			f.parent = num_frames_split + least_branch;
			const size_t num_branches_k0 = f.branches.size();
			for (size_t i = 0; i < num_branches_k0; ++i)
			{
				if (f.branches[i] == f.parent)
				{
					f.branches.erase(f.branches.begin() + i);
					break;
				}
			}
			for (size_t i = least_branch; i > 0; --i)
			{
				frame& f = child.frames[num_frames_split + i];
				f.parent = num_frames_split + i - 1;
				std::sort(f.branches.begin(), f.branches.end());
				BOOST_ASSERT(f.branches.front() == f.parent);
				f.branches.front() += 2; // Previously this is the parent frame, now it becomes the first branch frame, so the frame number difference is +2.
			}
			frame& f0 = child.frames[num_frames_split];
			f0.branches.insert(f0.branches.begin(), num_frames_split + 1);
		}
		child.frames[num_frames_split].parent = ma1.frame;

		// Copy part 2 of the current ligand to the child ligand, and update the frame numbers of parent, branches, and atom neighbors.
		for (size_t k = num_frames_split; k < num_frames_1; ++k)
		{
			child.frames.push_back(this->frames[k]);
			frame& f = child.frames.back();
			if (f.parent >= num_frames_split) f.parent += num_frames_2;
			const size_t num_branches = f.branches.size();
			for (size_t i = 0; i < num_branches; ++i)
			{
				if (f.branches[i] >= num_frames_split) f.branches[i] += num_frames_2;
			}
			const size_t num_atoms = f.atoms.size();
			for (size_t i = 0; i < num_atoms; ++i)
			{
				atom& a = f.atoms[i];
				const size_t num_neighbors = a.neighbors.size();
				for (size_t j = 0; j < num_neighbors; ++j)
				{
					if (a.neighbors[j].frame >= num_frames_split) a.neighbors[j].frame += num_frames_2;
				}
			}
		}

		// Remove mutable atom 2 from ma2 frame.
		frame& f4 = child.frames[num_frames_split];
		f4.atoms.erase(f4.atoms.begin() + ma2.index);

		// The number of mutable atoms of child ligand is equal to the sum of its parent ligands minus 2.
		const size_t num_mutatable_atoms_1 = this->mutable_atoms.size();
		const size_t num_mutatable_atoms_2 = other.mutable_atoms.size();
		child.mutable_atoms.reserve(num_mutatable_atoms_1 + num_mutatable_atoms_2 - 2);
		
		// Copy the mutable atoms of the current ligand except ma1 to the child ligand.		
		for (size_t i = 0; i < num_mutatable_atoms_1; ++i)
		{
			if (this->mutable_atoms[i] != ma1)
			{
				child.mutable_atoms.push_back(this->mutable_atoms[i]);
				atom_index& ma = child.mutable_atoms.back();
				if (ma.frame > num_frames_split) ma.frame += num_frames_2;
			}
		}
		
		// Copy the mutable atoms of the other ligand except ma2 to the child ligand.
		for (size_t i = 0; i < num_mutatable_atoms_2; ++i)
		{
			if (this->mutable_atoms[i] != ma2)
			{
				child.mutable_atoms.push_back(other.mutable_atoms[i]);
				atom_index& ma = child.mutable_atoms.back();
				ma.frame = num_frames_split + mapping[ma.frame];
			}
		}
		
		// Determine if the child ligand is able to perform mutation or crossover.
		child.mutation_feasible = !child.mutable_atoms.empty();
		child.crossover_feasible = child.frames.size() > 1;

		return new ligand(child);
	}

	void ligand::evaluate_efficacy(const fl free_energy)
	{
		this->free_energy = free_energy;
		efficacy = free_energy / num_heavy_atoms;
	}

/*
	// add a fragment to the position of a randomly selected hydrogen

	void ligand::mutate(ligand fragment)
	{
		// select hydrogen
		int count(0), connectIndex, fragHydrogen(-1), fragIndex, linkerHydrogen(-1);
		// atom copy and atom reference
		atom curHydrogen, *fragHydrogenAtom, *fragConnectAtom;
		string element1, element2;
		Vec3d delta;

		linkerHydrogen = IndexOfRandomHydrogen();
		// obtain information on hydrogens
		curHydrogen = atoms[linkerHydrogen];
		// hydrogen has only 1 connection, get the first index in array
		connectIndex = *(curHydrogen.IndexArray.begin());
		const atom connectAtom = atoms[connectIndex];
		// do the same on the fragment, there must be so hydrogen in the library...
		fragHydrogen = fragment.IndexOfRandomHydrogen();
		// these atoms would be updated, get reference instead of copying
		fragHydrogenAtom = &fragment.atoms[fragHydrogen];
		fragIndex = *(fragHydrogenAtom->IndexArray.begin());
		fragConnectAtom = &fragment.atoms[fragIndex];

		// move fragment in place
		// pick the bond to be adjusted
		element1 = connectAtom.element;
		element2 = fragConnectAtom->element;
		delta = curHydrogen.coordinates - connectAtom.coordinates;
		delta.normalize();
		bond_library library;
		// scale it to reflect real molecular bond length
		delta *= library.length(element1, element2);
		// move fragment in place
		fragment.Translate(fragIndex, (connectAtom.coordinates + delta));

		// rotate fragment to minimise distance between fragment's hydrogen and connected atom
		// Rotating the fragment so it has the correct orientation...

		// use geometry to calculate the angle necessary to achieve minimum distance
		Vec3d A = connectAtom.coordinates - fragConnectAtom->coordinates;
		Vec3d B = fragHydrogenAtom->coordinates - fragConnectAtom->coordinates;
		double dist = (fragHydrogenAtom->coordinates - connectAtom.coordinates).length();
		// normal to the 3 points
		Vec3d normal = B^A;
		// angle between 3 points with respect to fragment connected atom
		double angle = acos((A * B) / sqrt(A.length2() * B.length2()));
		// rotate along the line using connected atom as pivot
		// prevent the situation when it is already in place in the beginning
		if (normal != Vec3d())
			fragment.RotateLine(normal, fragIndex, angle);
		// rotate the other side if the dist increased
		if ((fragHydrogenAtom->coordinates - connectAtom.coordinates).length() > dist)
			fragment.RotateLine(normal, fragIndex, -2 * angle);

		// SP2 bond pair appears in planar
		// The fragment bond is SP2-SP2, so planarizing...
		if (connectAtom.isSP2() && fragConnectAtom->isSP2())
		{
			int cur_dihedral(-1), frag_dihedral(-1);
			double dihedral, dist;
			set<int>::iterator it;
			map<double, double> best_angles;
			// find reference atoms for both connected atoms
			for (it = connectAtom.IndexArray.begin(); it != connectAtom.IndexArray.end(); ++it)
			{
				if (*it != linkerHydrogen)
					cur_dihedral = *it;
			}
			for (it = fragConnectAtom->IndexArray.begin(); it != fragConnectAtom->IndexArray.end(); ++it)
			{
				if (*it != fragHydrogen)
					frag_dihedral = *it;
			}
			// calculate dihedral angle on these four atoms
			dihedral = DihedralAngle(atoms[cur_dihedral].coordinates, connectAtom.coordinates, (*fragConnectAtom).coordinates, fragment.atoms[frag_dihedral].coordinates);
			// rotate so that the two planes are now parallel
			fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, -dihedral);
			// calculate distance in this orientation
			dist = MolecularDistance(fragment);
			best_angles.insert(pair<double, double>(dist, -dihedral));
			// try the opposite
			fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, pi);
			dist = MolecularDistance(fragment);
			best_angles.insert(pair<double, double>(dist, pi - dihedral));
			// try the other set of orientations
			fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, -pi + 2 * dihedral);
			dist = MolecularDistance(fragment);
			best_angles.insert(pair<double, double>(dist, dihedral));
			fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, pi);
			dist = MolecularDistance(fragment);
			best_angles.insert(pair<double, double>(dist, pi + dihedral));
			fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, -pi - dihedral);
			fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, best_angles.rbegin()->second); // maximize intra-molecular distance
		}
			// minimise hinderance for SP3 and SP bond
		else
		{
			// Rotating fragment around connecting bond to minimize steric hindrance...
			double BadContact, BestContact(0), BestAngle(0), Angle(0), delta(pi / 25);
			Vec3d normal = connectAtom.coordinates - fragConnectAtom->coordinates;
			// try a whole circle with definite steps
			while (Angle < 2 * pi)
			{
				Angle += delta;
				if (normal == Vec3d()) break; // impossible, useless
				fragment.RotateLine(normal, fragIndex, delta);
				BadContact = MolecularDistance(fragment);
				// maximize the distance
				if (BadContact > BestContact)
				{
					BestContact = BadContact;
					BestAngle = Angle;
				}
			}
			// rotate to least hindered orientation
			if (normal != Vec3d())
				fragment.RotateLine(normal, fragIndex, BestAngle);
		}

		// Merging the fragment with the original molecule...
		// remove hydrogen atom from both molecule
		DeleteAtom(linkerHydrogen);
		fragment.DeleteAtom(fragHydrogen);

		int updateIndex, cascadeIndex(MaxIndex());
		set<int> tempIndice;
		ostringstream output;

		for (map<int, atom>::iterator it = fragment.atoms.begin(); it != fragment.atoms.end(); ++it)
		{
			// copy the atom
			atom toAdd(it->second);
			updateIndex = atoi(toAdd.PDBIndex.c_str());
			// append to the largest index of this molecule
			updateIndex += cascadeIndex;
			output.str(string());
			output << updateIndex;
			toAdd.PDBIndex = output.str();
			tempIndice.clear();
			// produce a set of pending indice
			for (set<int>::iterator iter = toAdd.IndexArray.begin(); iter != toAdd.IndexArray.end(); ++iter)
				tempIndice.insert(*iter + cascadeIndex);
			// clear the original indice
			toAdd.IndexArray.clear();
			// move the indice to the original set
			for (set<int>::iterator iter = tempIndice.begin(); iter != tempIndice.end(); ++iter)
				toAdd.IndexArray.insert(*iter);
			// add atom to the current molecule
			atoms.insert(pair<int, atom > (updateIndex, toAdd));
		}

		// at last connect the fragment and molecule at the selected position
		atoms[connectIndex].IndexArray.insert(fragIndex + cascadeIndex);
		atoms[fragIndex + cascadeIndex].IndexArray.insert(connectIndex);
	}

	void ligand::Translate(int index, Vec3d origin)
	{
		// calculate real displacement to move
		Vec3d delta = origin - atoms[index].coordinates;
		// move all atoms in the molecule
		for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
			it->second.coordinates = it->second.coordinates + delta;
	}

	void ligand::RotateLine(Vec3d v1, Vec3d v2, int index, double radian)
	{
		Vec3d location, delta = atoms[index].coordinates;
		Mat4d rot;
		// create rotation matrix
		rot = rot.createRotation(radian, v1 - v2);
		for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
		{
			// shift to indexed atom position
			location = it->second.coordinates - delta;
			// perform rotation
			location = rot * location;
			// move back by delta amount
			it->second.coordinates = location + delta;
		}
	}

	// same method by calculating normal beforehand

	void ligand::RotateLine(Vec3d normal, int index, double radian)
	{
		Vec3d location, delta = atoms[index].coordinates;
		Mat4d rot;
		rot = rot.createRotation(radian, normal);
		for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
		{
			location = it->second.coordinates - delta;
			location = rot * location;
			it->second.coordinates = location + delta;
		}
	}

	// total minimum distance for all atoms, quadratic

	double ligand::MolecularDistance(ligand& other)
	{
		double min, dist, total_dist = 0;
		map<int, atom>::iterator it1, it2;
		for (it1 = atoms.begin(); it1 != atoms.end(); ++it1)
		{
			min = DBL_MAX;
			for (it2 = other.atoms.begin(); it2 != other.atoms.end(); ++it2)
			{
				dist = it1->second.DistanceTo(it2->second);
				if (dist < min)
					min = dist;
			}
			total_dist += dist;
		}
		return total_dist;
	}

	// in radian, angle between 2 planes (a1-a2-a3, a2-a3-a4) with respect to vector a2-a3

	double ligand::DihedralAngle(const Vec3d& a1, const Vec3d& a2, const Vec3d& a3, const Vec3d& a4)
	{
		Vec3d v1 = a2 - a1;
		Vec3d v2 = a3 - a2;
		Vec3d v3 = a4 - a3;
		return atan2((v1 * v2.length())*(v2^v3), (v1^v2)*(v2^v3));
	}
*/
}
