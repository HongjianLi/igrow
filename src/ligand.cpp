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

	ligand::ligand(const path& p) : p(p), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0), connector1(0), connector2(0), logp(0) // TODO: comment logp(0)
	{
		// Initialize necessary variables for constructing a ligand.
		frames.reserve(30); // A ligand typically consists of <= 30 frames.
		frames.push_back(frame(0, 0, 0)); // ROOT is also treated as a frame. The parent and rotorX of ROOT frame are dummy.
		mutable_atoms.reserve(20); // A ligand typically consists of <= 20 mutable atoms.

		// Initialize helper variables for parsing.
		size_t current = 0; // Index of current frame, initialized to ROOT frame.
		frame* f = &frames.front(); // Pointer to the current frame.
		size_t num_lines = 0; // Used to track line number for reporting parsing errors, if any.
		string line; // A line of ligand file in pdbqt format.
		line.reserve(79); // According to PDBQT specification, the last item AutoDock4 atom type locates at 1-based [78, 79].

		// Parse ATOM/HETATM, BRANCH, ENDBRANCH.
		ifstream in(p); // Parsing starts. Open the file stream as late as possible.
		while (getline(in, line))
		{
			++num_lines;
			if (starts_with(line, "ATOM") || starts_with(line, "HETATM"))
			{
				// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
				BOOST_ASSERT(current == frames.size() - 1);

				// Parse the ATOM/HETATM line into an atom, which belongs to the current frame.
				atoms.push_back(atom(line));
				
				// Validate the AutoDock4 atom type.
				atom& a = atoms.back();
				if (a.ad == AD_TYPE_SIZE) throw parsing_error(p, num_lines, "Atom type " + line.substr(77, isspace(line[78]) ? 1 : 2) + " is not supported by igrow.");
				
				// Find the neighbors in the same frame.
				const size_t idx = atoms.size() - 1; // The index to the current atom.
				BOOST_ASSERT(f->begin <= idx);
				for (size_t i = f->begin; i < idx; ++i) // Exclude the last atom which is the current atom itself.
				{
					atom& b = atoms[i];
					if (a.is_neighbor(b))
					{
						a.neighbors.push_back(i);
						b.neighbors.push_back(idx);
						BOOST_ASSERT(a.neighbors.size() <= 4);
						BOOST_ASSERT(b.neighbors.size() <= 4);
					}
				}
				
				// Add rotorX as a neighbor because rotorY is assumed to be the first atom of a BRANCH frame.
				if ((current) && (idx == f->begin))
				{
					a.neighbors.push_back(f->rotorX);
					atoms[f->rotorX].neighbors.push_back(idx);
					BOOST_ASSERT(a.neighbors.size() <= 4);
					BOOST_ASSERT(atoms[f->rotorX].neighbors.size() <= 4);
				}
				
				if (a.is_mutable()) mutable_atoms.push_back(idx);
				if (!a.is_hydrogen()) ++num_heavy_atoms;
				if (a.is_hb_donor()) ++num_hb_donors;
				if (a.is_hb_acceptor()) ++num_hb_acceptors;
				mw += a.atomic_weight();
			}
			else if (starts_with(line, "BRANCH"))
			{
				// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
				// Y is not parsed because the atom immediately follows "BRANCH" must be Y in pdbqt files created by the prepare_ligand4.py script of MGLTools.
				// This assumption fails if pdbqt files are prepared by OpenBabel. In this case, the class frame should incorporate a new field rotorY to track the origin.
				const size_t x = right_cast<size_t>(line, 7, 10);

				// Find the corresponding heavy atom with X as its atom number in the current frame.
				for (size_t i = f->begin; ; ++i)
				{
					if (atoms[i].number == x)
					{
						// Insert a new frame whose parent is the current frame.
						frames.push_back(frame(current, i, atoms.size()));
						break;
					}
				}

				// Now the current frame is the newly inserted BRANCH frame.
				current = frames.size() - 1;

				// The parent frame has the current frame as one of its branches.
				f->branches.push_back(current);

				// Update the pointer to the current frame.
				f = &frames[current];
			}
			else if (starts_with(line, "ENDBRANCH"))
			{
				// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
				// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
				if (f->begin == atoms.size()) throw parsing_error(p, num_lines, "An empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

				// Now the parent of the following frame is the parent of current frame.
				current = f->parent;

				// Update the pointer to the current frame.
				f = &frames[current];
			}
		}
		in.close(); // Parsing finishes. Close the file stream as soon as possible.

		BOOST_ASSERT(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
		BOOST_ASSERT(f == &frames.front()); // The frame pointer should point to the ROOT frame.
		
		// Determine the number of atoms.
		num_atoms = atoms.size();
		BOOST_ASSERT(num_atoms >= num_heavy_atoms);
		BOOST_ASSERT(num_atoms <= num_heavy_atoms + mutable_atoms.size());

		// Determin the number of rotatable bonds.
		num_rotatable_bonds = frames.size() - 1;
		BOOST_ASSERT(num_atoms + (num_rotatable_bonds << 1) + 3 <= num_lines); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines + REMARK lines (if any) == num_lines

		// Set frames[i].end = frames[i + 1].begin
		for (size_t i = 0; i < num_rotatable_bonds; ++i)
		{
			frames[i].end = frames[i + 1].begin;
		}
		frames.back().end = num_atoms;
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
			for (size_t i = f.begin; i < f.end; ++i)
			{
				const atom& a = atoms[i];
				out << "ATOM  " << setw(5) << a.number << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << '\n';
			}
		}
		out << "ENDROOT\n";

		// Dump the BRANCH frames.
		vector<bool> dump_branches(frames.size()); // dump_branches[0] is dummy. The ROOT frame has been dumped.
		vector<size_t> stack; // Stack to track the new_to_old_mapping sequence of frames in order to avoid recursion.
		stack.reserve(num_rotatable_bonds); // The ROOT frame is excluded.
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
				out << "BRANCH"    << setw(4) << atoms[f.rotorX].number << setw(4) << atoms[f.begin].number << '\n';
				for (size_t i = f.begin; i < f.end; ++i)
				{
					const atom& a = atoms[i];
					out << "ATOM  " << setw(5) << a.number << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << '\n';
				}
				dump_branches[fn] = true;
				for (size_t i = f.branches.size(); i > 0;)
				{
					stack.push_back(f.branches[--i]);
				}
			}
			else // This BRANCH frame has been dumped.
			{
				out << "ENDBRANCH" << setw(4) << atoms[f.rotorX].number << setw(4) << atoms[f.begin].number << '\n';
				stack.pop_back();
			}
		}
		out << "TORSDOF " << num_rotatable_bonds << '\n';
		out.close();
	}

	ligand::ligand(const ligand& l1, const ligand& l2, const mt19937eng& eng, const operation op) : parent1(l1.p), parent2(l2.p)
	{
		switch (op)
		{
		case operation_mutation:
			mutate(l1, l2, eng);
			break;
		case operation_crossover:
			crossover(l1, l2, eng);
			break;
		}
	}

	void ligand::mutate(const ligand& l1, const ligand& l2, const mt19937eng& eng)
	{
		// Initialize random number generators for obtaining two random mutable atoms.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_mutable_atom_gen_1(eng, uniform_int_distribution<size_t>(0, l1.mutable_atoms.size() - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_mutable_atom_gen_2(eng, uniform_int_distribution<size_t>(0, l2.mutable_atoms.size() - 1));

		// Obtain a random mutable atom from ligands 1 and 2 respectively.
		const size_t m1idx = l1.mutable_atoms[uniform_mutable_atom_gen_1()];
		const size_t m2idx = l2.mutable_atoms[uniform_mutable_atom_gen_2()];
		const atom& m1 = l1.atoms[m1idx]; // Constant reference to the mutable atom of ligand 1.
		const atom& m2 = l2.atoms[m2idx]; // Constant reference to the mutable atom of ligand 2.
		
		// The mutable atom should consist of only one non-rotatable single bond.
		BOOST_ASSERT(m1.neighbors.size() == 1);
		BOOST_ASSERT(m2.neighbors.size() == 1);

		// Obtain the frames of the two mutable atoms.
		const size_t f1idx = l1.get_frame(m1idx);
		const size_t f2idx = l2.get_frame(m2idx);
		const frame& f1 = l1.frames[f1idx];
		const frame& f2 = l2.frames[f2idx];

		// The mutable atom cannot be rotorY, but can be the last atom.
		BOOST_ASSERT(f1.begin <  m1idx);
		BOOST_ASSERT(f1.end   >= m1idx);
		BOOST_ASSERT(f2.begin <  m2idx);
		BOOST_ASSERT(f2.end   >= m2idx);

		// Obtain constant references to the connector atoms.
		const size_t c1idx = m1.neighbors.front();
		const size_t c2idx = m2.neighbors.front();
		const atom& c1 = l1.atoms[c1idx];
		const atom& c2 = l2.atoms[c2idx];
		
		// Both the connector and mutable atoms should be in the same frame.
		BOOST_ASSERT(f1idx == l1.get_frame(c1idx));
		BOOST_ASSERT(f2idx == l2.get_frame(c2idx));

		// Set the connector atoms.
		connector1 = c1.number;
		connector2 = c2.number;
		
		// The number of rotatable bonds of child ligand is equal to the sum of its parent ligands plus 1.
		num_rotatable_bonds = l1.num_rotatable_bonds + l2.num_rotatable_bonds + 1;
		
		// The number of atoms of child ligand is equal to the sum of its parent ligands minus 2.
		num_atoms = l1.num_atoms + l2.num_atoms - 2;

		// The number of heavy atoms of child ligand is equal to the sum of its parent ligands minus the two mutable atoms if they are heavy atoms.
		num_heavy_atoms = l1.num_heavy_atoms + l2.num_heavy_atoms - ((m1.is_hydrogen() ? 0 : 1) + (m2.is_hydrogen() ? 0 : 1));

		// The number of hydrogen bond donors of child ligand is equal to the sum of its parent ligands minus the two mutable atoms if they are hydrogen bond donors.
		num_hb_donors = l1.num_hb_donors + l2.num_hb_donors - ((m1.is_hb_donor() ? 1 : 0) + (m2.is_hb_donor() ? 1 : 0));

		// The number of hydrogen bond acceptors of child ligand is equal to the sum of its parent ligands.
		num_hb_acceptors = l1.num_hb_acceptors + l2.num_hb_acceptors;

		// The molecular weight of child ligand is equal to the sum of its parent ligands minus the two mutable atoms.
		mw = l1.mw + l2.mw - (m1.atomic_weight() + m2.atomic_weight());

		// The logP of child ligand is equal to the sum of its parent ligands minus the two mutable atoms.
		logp = l1.logp + l2.logp; // TODO: comment this line.

		// Determine the number of ligand 1's frames that will be directly copied to the child ligand. This is also the new frame number of ligand 2's f2 frame.
		size_t l1_num_frames_split = 0; // Frames [0, l1_num_frames_split) constitute part 1, and frames [l1_num_frames_split, l1_num_frames) constitute part 2.
		size_t insertion_index; // The index at which the new frame number of ligand 2's f2 frame will be inserted into the branches of child ligand's f1 frame.
		const size_t l1_num_frames = l1.frames.size();
		const size_t f1_num_branches = f1.branches.size();
		for (size_t i = 0; i < f1_num_branches; ++i)
		{
			if (c1.number < l1.atoms[l1.frames[f1.branches[i]].rotorX].number)
			{
				l1_num_frames_split = f1.branches[i];
				BOOST_ASSERT(l1_num_frames_split); // At least one frame of ligand 1 will be directly copied to the child ligand.
				insertion_index = i;				
				break;
			}
		}
		if (!l1_num_frames_split) // The insertion index is at the end of the branches of child ligand's f1 frame.
		{
			insertion_index = f1.branches.size();
			
			// Trace back the parents of f1 frame to the immediately next branch number, which is equal to l1_num_frames_split.
			size_t c = f1idx; // Current frame.
			size_t p = f1.parent; // Parent frame.
			while (c && (c == l1.frames[p].branches.back()))
			{
				c = p;
				p = l1.frames[p].parent;
			}
			if (c) // The immediately next branch number is found in frame p.
			{
				const frame& pf = l1.frames[p];
				const size_t pf_num_branches = pf.branches.size();
				for (size_t i = 0; i < pf_num_branches; ++i)
				{
					if (pf.branches[i] == c)
					{
						BOOST_ASSERT(i + 1 < pf_num_branches);
						l1_num_frames_split = pf.branches[i + 1];
						break;
					}
				}
			}
			else // The entire frame tree of ligand 1 has been traversed.
			{
				l1_num_frames_split = l1_num_frames;
			}
		}
		BOOST_ASSERT(l1_num_frames_split >= f1idx + 1);
		BOOST_ASSERT(l1_num_frames_split <= l1_num_frames);

		// Reserve enough capacity for storing frames and atoms.
		const size_t l2_num_frames = l2.frames.size();
		frames.reserve(l1_num_frames + l2_num_frames);
		atoms.reserve(num_atoms);

		// Copy part 1 of ligand 1 to the child ligand.
		// TODO: use std::copy instead.
		for (size_t i = 0; i < l1_num_frames_split; ++i)
		{
			frames.push_back(l1.frames[i]);
		}
		--frames[f1idx].end;
		for (size_t i = f1idx + 1; i < l1_num_frames_split; ++i)
		{
			frame& f = frames[i];
			--f.begin;
			--f.end;
			if (f.rotorX > m1idx) -- f.rotorX;
		}
		for (size_t i = 0; i < f1.begin; ++i)
		{
			atoms.push_back(l1.atoms[i]);
		}
		for (size_t i = f1.begin; i < l1.frames[l1_num_frames_split - 1].end; ++i)
		{
			if (i == m1idx) continue;
			atoms.push_back(l1.atoms[i]);
			atom& a = atoms.back();
			const size_t a_num_neighbors = a.neighbors.size();
			for (size_t j = 0; j < a_num_neighbors; ++j)
			{
				if (a.neighbors[j] > m1idx) --a.neighbors[j];
			}
		}

		// Insert the new frame number of ma2 frame of the other ligand at the previously found index.
		frame& f0 = frames[f1idx];
		f0.branches.insert(f0.branches.begin() + insertion_index, l1_num_frames_split);

		// Increment the branch numbers of part 2 frames by l2_num_frames.
		const size_t f0_num_branches = f0.branches.size();
		BOOST_ASSERT(f0_num_branches == f1_num_branches + 1);
		for (size_t i = insertion_index + 1; i < f0_num_branches; ++i)
		{
			f0.branches[i] += l2_num_frames;
		}
		for (size_t k = 0; k < f1idx; ++k)
		{
			frame& f = frames[k];
			const size_t f_num_branches = f.branches.size();
			for (size_t i = 0; i < f_num_branches; ++i)
			{
				if (f.branches[i] >= l1_num_frames_split) f.branches[i] += l2_num_frames;
			}
		}

		//// Update a connector neighbor of c3 from mutable atom 1 to the rotorY of ma2 frame of the other ligand.
		//atom& c3 = f3.atoms[m1.neighbors.front().index];
		//for (size_t i = c3.neighbors.size(); i > 0;)
		//{
		//	if (c3.neighbors[--i] == ma1)
		//	{
		//		c3.neighbors[i] = atom_index(l1_num_frames_split, 0); // Assume rotorY is always 0.
		//		break;
		//	}
		//}
		//
		//// Update atom neighbors because the indexes to all the atoms after ma2 are decremented by 1.
		//
		// Find the traversal sequence (i.e. new_to_old_mapping) of ligand 2 starting from f2 frame, as well as its reverse traversal sequence (i.e. old_to_new_mapping).
		vector<size_t> new_to_old_mapping;
		new_to_old_mapping.reserve(l2_num_frames);
		vector<size_t> old_to_new_mapping(l2_num_frames);
		vector<size_t> stack;
		stack.reserve(l2_num_frames);
		stack.push_back(f2idx);
		while (!stack.empty())
		{
			const size_t k = stack.back();
			old_to_new_mapping[k] = new_to_old_mapping.size();
			new_to_old_mapping.push_back(k);			
			const frame& f = l2.frames[k];
			stack.pop_back();
			for (size_t i = f.branches.size(); i > 0;)
			{
				const size_t j = f.branches[--i];
				if (std::find(new_to_old_mapping.cbegin(), new_to_old_mapping.cend(), j) == new_to_old_mapping.end()) stack.push_back(j);
			}
			if (std::find(new_to_old_mapping.cbegin(), new_to_old_mapping.cend(), f.parent) == new_to_old_mapping.end()) stack.push_back(f.parent);
		}		
		//
		//// Copy the other ligand to the child ligand, and update the frame numbers of parent, branches, and atom neighbors;
		//for (size_t k = 0; k < l2_num_frames; ++k)
		//{
		//	frames.push_back(l2.frames[new_to_old_mapping[k]]);
		//	frame& f = frames.back();
		//	f.parent = old_to_new_mapping[f.parent] + l1_num_frames_split;
		//	const size_t num_branches = f.branches.size();
		//	for (size_t i = 0; i < num_branches; ++i)
		//	{
		//		f.branches[i] = old_to_new_mapping[f.branches[i]] + l1_num_frames_split;
		//	}
		//	const size_t num_atoms = f.atoms.size();
		//	for (size_t i = 0; i < num_atoms; ++i)
		//	{
		//		atom& a = f.atoms[i];
		//		const size_t num_neighbors = a.neighbors.size();
		//		for (size_t j = 0; j < num_neighbors; ++j)
		//		{
		//			a.neighbors[j].frame = old_to_new_mapping[a.neighbors[j].frame] + l1_num_frames_split;
		//		}
		//	}
		//}
		//
		//// The parent path from ma1 in the other ligand becomes the branch path in the child ligand.
		//if (ma2.frame)
		//{
		//	const size_t k0 = old_to_new_mapping.front();
		//	BOOST_ASSERT(k0 > 0);
		//	const size_t least_branch = k0 - 1; // The least frame number of l1_num_frames_split + k0 frame of the child ligand is always equal to k0 - 1.
		//	frame& f = frames[l1_num_frames_split + k0];
		//	f.parent = l1_num_frames_split + least_branch;
		//	const size_t num_branches_k0 = f.branches.size();
		//	for (size_t i = 0; i < num_branches_k0; ++i)
		//	{
		//		if (f.branches[i] == f.parent)
		//		{
		//			f.branches.erase(f.branches.begin() + i);
		//			break;
		//		}
		//	}
		//	for (size_t i = least_branch; i > 0; --i)
		//	{
		//		frame& f = frames[l1_num_frames_split + i];
		//		f.parent = l1_num_frames_split + i - 1;
		//		std::sort(f.branches.begin(), f.branches.end());
		//		BOOST_ASSERT(f.branches.front() == f.parent);
		//		f.branches.front() += 2; // Previously this is the parent frame, now it becomes the first branch frame, so the frame number difference is +2.
		//	}
		//	frame& f0 = frames[l1_num_frames_split];
		//	f0.branches.insert(f0.branches.begin(), l1_num_frames_split + 1);
		//}
		//frames[l1_num_frames_split].parent = ma1.frame;

		//// Copy part 2 of the current ligand to the child ligand, and update the frame numbers of parent, branches, and atom neighbors.
		//for (size_t k = l1_num_frames_split; k < l1_num_frames; ++k)
		//{
		//	frames.push_back(l1.frames[k]);
		//	frame& f = frames.back();
		//	if (f.parent >= l1_num_frames_split) f.parent += l2_num_frames;
		//	const size_t num_branches = f.branches.size();
		//	for (size_t i = 0; i < num_branches; ++i)
		//	{
		//		if (f.branches[i] >= l1_num_frames_split) f.branches[i] += l2_num_frames;
		//	}
		//	const size_t num_atoms = f.atoms.size();
		//	for (size_t i = 0; i < num_atoms; ++i)
		//	{
		//		atom& a = f.atoms[i];
		//		const size_t num_neighbors = a.neighbors.size();
		//		for (size_t j = 0; j < num_neighbors; ++j)
		//		{
		//			if (a.neighbors[j].frame >= l1_num_frames_split) a.neighbors[j].frame += l2_num_frames;
		//		}
		//	}
		//}

		//// Remove mutable atom 2 from ma2 frame.
		//frame& f4 = frames[l1_num_frames_split];
		//f4.atoms.erase(f4.atoms.begin() + ma2.index);
		//
		//// Update c4 neighbor.
		//
		//// Update atom neighbors because the indexes to all the atoms after ma2 are decremented by 1.

		//// The number of mutable atoms of child ligand is equal to the sum of its parent ligands minus 2.
		//const size_t num_mutatable_atoms_1 = l1.mutable_atoms.size();
		//const size_t num_mutatable_atoms_2 = l2.mutable_atoms.size();
		//mutable_atoms.reserve(num_mutatable_atoms_1 + num_mutatable_atoms_2 - 2);
		//
		//// Copy the mutable atoms of the current ligand except ma1 to the child ligand.		
		//for (size_t i = 0; i < num_mutatable_atoms_1; ++i)
		//{
		//	if (l1.mutable_atoms[i] != ma1)
		//	{
		//		mutable_atoms.push_back(l1.mutable_atoms[i]);
		//		atom_index& ma = mutable_atoms.back();
		//		if (ma.frame > l1_num_frames_split) ma.frame += l2_num_frames;
		//	}
		//}
		//
		//// Copy the mutable atoms of the other ligand except ma2 to the child ligand.
		//for (size_t i = 0; i < num_mutatable_atoms_2; ++i)
		//{
		//	if (l2.mutable_atoms[i] != ma2)
		//	{
		//		mutable_atoms.push_back(l2.mutable_atoms[i]);
		//		atom_index& ma = mutable_atoms.back();
		//		ma.frame = l1_num_frames_split + old_to_new_mapping[ma.frame];
		//	}
		//}
	}
	
	void ligand::crossover(const ligand& l1, const ligand& l2, const mt19937eng& eng)
	{
		// Initialize random number generators for obtaining two random mutable atoms.
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_rotatable_bond_gen_1(eng, uniform_int_distribution<size_t>(0, l1.num_rotatable_bonds - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_rotatable_bond_gen_2(eng, uniform_int_distribution<size_t>(0, l2.num_rotatable_bonds - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t>> uniform_01_gen(eng, uniform_int_distribution<size_t>(0, 1));
		
		// Obtain a random rotatable bond from the current ligand and the other ligand respectively.
		const size_t rotatable_bond_1 = uniform_rotatable_bond_gen_1();
		const size_t rotatable_bond_2 = uniform_rotatable_bond_gen_2();

		// Initialize a child ligand.
		parent1 = l1.p;
		parent2 = l2.p;
		//connector1 = c1.number;
		//connector2 = c2.number;
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
