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
#include "fstream.hpp"
#include "mat3.hpp"
#include "ligand.hpp"

namespace igrow
{
	ligand::ligand(const path& p) : p(p), connector1(0), connector2(0), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0), logp(0) // TODO: comment logp(0)
	{
		// Initialize necessary variables for constructing a ligand.
		frames.reserve(30); // A ligand typically consists of <= 30 frames.
		frames.push_back(frame(0, 0, 0, 0)); // ROOT is also treated as a frame. The parent, rotorX, and rotorY of ROOT frame are dummy.
		frames.back().branches.reserve(4); // A frame typically consists of <= 4 BRANCH frames.
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

				// Validate the AutoDock4 atom type.
				const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
				const size_t ad = parse_ad_type_string(ad_type_string);
				if (ad == AD_TYPE_SIZE) throw parsing_error(p, num_lines, "Atom type " + ad_type_string + " is not supported by igrow.");

				// Parse the ATOM/HETATM line into an atom, which belongs to the current frame.
				atoms.push_back(atom(line.substr(12, 18), line.substr(54), right_cast<size_t>(line, 7, 11), vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54)), ad));

				// Update ligand properties.
				const atom& a = atoms.back();
				if (a.is_mutable()) mutable_atoms.push_back(a.number);
				if (!a.is_hydrogen()) ++num_heavy_atoms;
				if (a.is_hb_donor()) ++num_hb_donors;
				if (a.is_hb_acceptor()) ++num_hb_acceptors;
				mw += a.atomic_weight();
			}
			else if (starts_with(line, "BRANCH"))
			{
				// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
				frames.push_back(frame(current, right_cast<size_t>(line, 7, 10), right_cast<size_t>(line, 11, 14), atoms.size()));

				// Now the current frame is the newly inserted BRANCH frame.
				current = frames.size() - 1;

				// The parent frame has the current frame as one of its branches.
				f->branches.push_back(current);

				// Update the pointer to the current frame.
				f = &frames[current];

				// Reserve enough capacity for storing BRANCH frames.
				f->branches.reserve(4); // A frame typically consists of <= 4 BRANCH frames.
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

		// Determine the number of rotatable bonds.
		num_rotatable_bonds = frames.size() - 1;
		BOOST_ASSERT(num_atoms + (num_rotatable_bonds << 1) + 3 <= num_lines); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines + REMARK lines (if any) == num_lines

		// Determine the maximum atom serial number.
		max_atom_number = atoms.back().number;
		BOOST_ASSERT(max_atom_number >= num_atoms);

		// Set frames[i].end = frames[i + 1].begin
		for (size_t i = 0; i < num_rotatable_bonds; ++i)
		{
			frames[i].end = frames[i + 1].begin;
		}
		frames.back().end = num_atoms;
	}

	void ligand::save(const path& p) const
	{
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
		vector<size_t> stack; // Stack to track the l4_to_l2_mapping sequence of frames in order to avoid recursion.
		stack.reserve(num_rotatable_bonds); // The ROOT frame is excluded.
		{
			const frame& f = frames.front();
			for (auto i = f.branches.rbegin(); i < f.branches.rend(); ++i)
			{
				stack.push_back(*i);
			}
		}
		while (!stack.empty())
		{
			const size_t fn = stack.back();
			const frame& f = frames[fn];
			if (!dump_branches[fn]) // This BRANCH frame has not been dumped.
			{
				out << "BRANCH"    << setw(4) << f.rotorX << setw(4) << f.rotorY << '\n';
				for (size_t i = f.begin; i < f.end; ++i)
				{
					const atom& a = atoms[i];
					out << "ATOM  " << setw(5) << a.number << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << '\n';
				}
				dump_branches[fn] = true;
				for (auto i = f.branches.rbegin(); i < f.branches.rend(); ++i)
				{
					stack.push_back(*i);
				}
			}
			else // This BRANCH frame has been dumped.
			{
				out << "ENDBRANCH" << setw(4) << f.rotorX << setw(4) << f.rotorY << '\n';
				stack.pop_back();
			}
		}
		out << "TORSDOF " << num_rotatable_bonds << '\n';
		out.close();
	}

	std::pair<size_t, size_t> ligand::get_frame(const size_t srn) const
	{
		BOOST_ASSERT(num_rotatable_bonds == frames.size() - 1);
		for (size_t k = 0; k <= num_rotatable_bonds; ++k)
		{
			const frame& f = frames[k];
			const size_t srn_begin = atoms[f.begin].number;
			const size_t srn_end = atoms[f.end - 1].number;
			BOOST_ASSERT(srn_begin <= srn_end);
			if ((f.end - f.begin) == (srn_end - srn_begin + 1)) // The serial numbers are continuous, which is the most cases.
			{
				if ((srn_begin <= srn) && (srn <= srn_end)) return std::pair<size_t, size_t>(k, f.begin + srn - srn_begin);
			}
			else // The serial numbers are not continuous, but they are sorted. Binary search can be used.
			{
				// Linear search at the moment.
				for (size_t i = f.begin; i < f.end; ++i)
				{
					if (srn == atoms[i].number) return std::pair<size_t,  size_t>(k, i);
				}
			}
		}
		throw std::domain_error("Failed to find an atom with serial number " + lexical_cast<string>(srn));
	}

	ligand::ligand(const ligand& l1, const ligand& l2, const size_t g1, const size_t g2) : parent1(l1.p), parent2(l2.p)
	{
		BOOST_ASSERT(g1 < l1.mutable_atoms.size());
		BOOST_ASSERT(g2 < l2.mutable_atoms.size());
		const size_t m1srn = l1.mutable_atoms[g1];
		const size_t m2srn = l2.mutable_atoms[g2];
		BOOST_ASSERT(m1srn >= 1);
		BOOST_ASSERT(m1srn <= l1.max_atom_number);
		BOOST_ASSERT(m2srn >= 1);
		BOOST_ASSERT(m2srn <= l2.max_atom_number);

		// Obtain the frames and indices of the two mutable atoms.
		const std::pair<size_t, size_t> p1 = l1.get_frame(m1srn);
		const std::pair<size_t, size_t> p2 = l2.get_frame(m2srn);
		const size_t f1idx = p1.first;
		const size_t f2idx = p2.first;
		const frame& f1 = l1.frames[f1idx];
		const frame& f2 = l2.frames[f2idx];
		const size_t m1idx = p1.second;
		const size_t m2idx = p2.second;
		BOOST_ASSERT(f1.begin <= m1idx);
		BOOST_ASSERT(f1.end   >  m1idx);
		BOOST_ASSERT(f2.begin <= m2idx);
		BOOST_ASSERT(f2.end   >  m2idx);

		const atom& m1 = l1.atoms[m1idx]; // Constant reference to the mutable atom of ligand 1.
		const atom& m2 = l2.atoms[m2idx]; // Constant reference to the mutable atom of ligand 2.
		BOOST_ASSERT(m1.is_mutable());
		BOOST_ASSERT(m2.is_mutable());

		// Find the connector atom that is covalently bonded to the mutable atom for both ligands.
		size_t c1idx, c2idx;
		for (c1idx = f1.begin; (c1idx == m1idx) || (!m1.is_neighbor(l1.atoms[c1idx])); ++c1idx);
		for (c2idx = f2.begin; (c2idx == m2idx) || (!m2.is_neighbor(l2.atoms[c2idx])); ++c2idx);

		// The connector atom should be in the same frame as the mutable atom.
		BOOST_ASSERT(c1idx < f1.end);
		BOOST_ASSERT(c2idx < f2.end);

		// Obtain constant references to the connector atoms.
		const atom& c1 = l1.atoms[c1idx];
		const atom& c2 = l2.atoms[c2idx];

		// Both the connector and mutable atoms should be in the same frame.
		BOOST_ASSERT(f1idx == l1.get_frame(c1.number).first);
		BOOST_ASSERT(f2idx == l2.get_frame(c2.number).first);

		// Set the connector atoms.
		connector1 = c1.number;
		connector2 = c2.number;

		// The maximum atom serial number of child ligand is equal to the sum of its parent ligands.
		max_atom_number = l1.max_atom_number + l2.max_atom_number;

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
		logp = l1.logp + l2.logp;

		// Reserve enough capacity for storing atoms.
		atoms.reserve(num_atoms);

		// Reserve enough capacity for storing frames.
		const size_t l1_num_frames = l1.frames.size();
		const size_t l2_num_frames = l2.frames.size();
		frames.reserve(l1_num_frames + l2_num_frames);

		// Determine the number of ligand 1's frames up to f1.
		const size_t f1_num_frames = f1idx + 1;
		BOOST_ASSERT(f1_num_frames <= l1_num_frames);

		// Create new frames for ligand 1's frames that are before f1.
		for (size_t k = 0; k < f1idx; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				f.branches.push_back(b > f1idx ? l2_num_frames + b : b);
			}

			// Populate atoms.
			BOOST_ASSERT(f.begin == rf.begin);
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Create a new frame for ligand 1's f1 frame itself.
		{
			// The reference frame is f1.
			const size_t f1_num_branches = f1.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(f1.parent, f1.rotorX, f1.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(1 + f1_num_branches); // This frame exactly consists of 1 + f1_num_branches BRANCH frames.
			f.branches.push_back(f1_num_frames);
			for (size_t i = 0; i < f1_num_branches; ++i)
			{
				f.branches.push_back(l2_num_frames + f1.branches[i]);
			}

			// Populate atoms.
			BOOST_ASSERT(f.begin == f1.begin);
			for (size_t i = f1.begin; i < m1idx; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			for (size_t i = m1idx + 1; i < f1.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Find the traversal sequence (i.e. l4_to_l2_mapping) of ligand 2 starting from f2 frame, as well as its reverse traversal sequence (i.e. l2_to_l4_mapping).
		vector<size_t> l4_to_l2_mapping;
		l4_to_l2_mapping.reserve(l2_num_frames);
		vector<size_t> l2_to_l4_mapping(l2_num_frames);
		{
			vector<size_t> stack;
			stack.reserve(l2_num_frames);
			stack.push_back(f2idx);
			while (!stack.empty())
			{
				const size_t k = stack.back();
				stack.pop_back();
				l2_to_l4_mapping[k] = l4_to_l2_mapping.size();
				l4_to_l2_mapping.push_back(k);
				const frame& rf = l2.frames[k];
				for (auto i = rf.branches.rbegin(); i < rf.branches.rend(); ++i)
				{
					if (std::find(l4_to_l2_mapping.begin(), l4_to_l2_mapping.end(), *i) == l4_to_l2_mapping.end()) stack.push_back(*i);
				}
				if (std::find(l4_to_l2_mapping.begin(), l4_to_l2_mapping.end(), rf.parent) == l4_to_l2_mapping.end()) stack.push_back(rf.parent);
			}
		}
		BOOST_ASSERT(l4_to_l2_mapping.size() == l2_num_frames);
		BOOST_ASSERT(l4_to_l2_mapping[0] == f2idx);
		BOOST_ASSERT(l2_to_l4_mapping[f2idx] == 0);

		// Calculate the translation vector for moving ligand 2 to a nearby place of ligand 1.
		const vec3 c1_to_c2 = ((c1.covalent_radius() + c2.covalent_radius()) / (c1.covalent_radius() + m1.covalent_radius())) * (m1.coordinate - c1.coordinate); // Vector pointing from c1 to the new position of c2.
		const vec3 origin_to_c2 = c1.coordinate + c1_to_c2; // Translation vector to translate ligand 2 from origin to the new position of c2.
		const vec3 c2_to_c1_nd = (-1 * c1_to_c2).normalize(); // Normalized vector pointing from c2 to c1.
		const vec3 c2_to_m2_nd = (m2.coordinate - c2.coordinate).normalize(); // Normalized vector pointing from c2 to m2.
		const mat3 rot(cross_product(c2_to_m2_nd, c2_to_c1_nd).normalize(), c2_to_m2_nd * c2_to_c1_nd); // Rotation matrix to rotate m2 along the normal to the direction from the new position of c2 to c1.

		// Create a new frame for ligand 2's f2 frame itself. Its branches are separately considered, depending on whether f2 is the ROOT frame of ligand 2.
		{
			// The reference frame is f2.
			BOOST_ASSERT(&f2 == &l2.frames[l4_to_l2_mapping[0]]);

			// Create a new frame based on the reference frame.
			frames.push_back(frame(f1idx, connector1, l1.max_atom_number + connector2, atoms.size()));
			frame& f = frames.back();

			// Populate atoms.
			for (size_t i = f2.begin; i < m2idx; ++i)
			{
				const atom& ra = l2.atoms[i];
				atoms.push_back(atom(ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.number, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
			}
			for (size_t i = m2idx + 1; i < f2.end; ++i)
			{
				const atom& ra = l2.atoms[i];
				atoms.push_back(atom(ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.number, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		if (!f2idx) // f2 is the ROOT frame of ligand 2.
		{
			BOOST_ASSERT(l2_to_l4_mapping[0] == 0);
			const size_t f2_num_branches = f2.branches.size();
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(f2_num_branches); // This frame exactly consists of f2_num_branches BRANCH frames.
			for (size_t i = 0; i < f2_num_branches; ++i)
			{
				f.branches.push_back(f1_num_frames + l2_to_l4_mapping[f2.branches[i]]);
			}
		}
		else // f2 is not the ROOT frame of ligand 2.
		{
			{
				const size_t f2_num_branches = f2.branches.size();
				frame& f = frames.back();

				// Populate branches.
				f.branches.reserve(1 + f2_num_branches); // This frame exactly consists of 1 + f2_num_branches BRANCH frames.
				f.branches.push_back(f1_num_frames + l2_to_l4_mapping[f2.parent]);
				for (size_t i = 0; i < f2_num_branches; ++i)
				{
					f.branches.push_back(f1_num_frames + l2_to_l4_mapping[f2.branches[i]]);
				}
			}

			// Create new frames for ligand 2's frames that are parent frames of f2 except ROOT.
			BOOST_ASSERT(l2_to_l4_mapping[0] >= 1);
			for (size_t k = 1; k < l2_to_l4_mapping.front(); ++k)
			{
				// Obtain a constant reference to the corresponding frame of ligand 2.
				const frame& rf = l2.frames[l4_to_l2_mapping[k]];
				const size_t rf_num_branches = rf.branches.size();
				BOOST_ASSERT(rf_num_branches >= 1);

				// Create a new frame based on the reference frame.
				const frame& pf = l2.frames[l4_to_l2_mapping[k - 1]];
				BOOST_ASSERT(f1idx + k == frames.size() - 1);
				frames.push_back(frame(f1idx + k, l1.max_atom_number + pf.rotorY, l1.max_atom_number + pf.rotorX, atoms.size()));
				frame& f = frames.back();

				// Populate branches.
				f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
				f.branches.push_back(f1_num_frames + l2_to_l4_mapping[rf.parent]);
				const size_t b = l4_to_l2_mapping[k - 1];
				for (size_t i = 0; i < rf_num_branches; ++i)
				{
					if (rf.branches[i] == b) continue;
					f.branches.push_back(f1_num_frames + l2_to_l4_mapping[rf.branches[i]]);
				}

				// Populate atoms.
				for (size_t i = rf.begin; i < rf.end; ++i)
				{
					const atom& ra = l2.atoms[i];
					atoms.push_back(atom(ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.number, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
				}
				f.end = atoms.size();
				BOOST_ASSERT(f.begin < f.end);
			}

			// Create new frames for ligand 2's ROOT frame.
			{
				// Obtain a constant reference to the corresponding frame of ligand 2.
				const frame& rf = l2.frames.front();
				const size_t rf_num_branches = rf.branches.size();
				BOOST_ASSERT(rf_num_branches >= 1);

				// Create a new frame based on the reference frame.
				const frame& pf = l2.frames[l4_to_l2_mapping[l2_to_l4_mapping.front() - 1]];
				BOOST_ASSERT(f1idx + l2_to_l4_mapping.front() == frames.size() - 1);
				frames.push_back(frame(f1idx + l2_to_l4_mapping.front(), l1.max_atom_number + pf.rotorY, l1.max_atom_number + pf.rotorX, atoms.size()));
				frame& f = frames.back();

				// Populate branches.
				f.branches.reserve(rf_num_branches - 1); // This frame exactly consists of rf_num_branches - 1 BRANCH frames.
				const size_t b = l4_to_l2_mapping[l2_to_l4_mapping.front() - 1];
				for (size_t i = 0; i < rf_num_branches; ++i)
				{
					if (rf.branches[i] == b) continue;
					f.branches.push_back(f1_num_frames + l2_to_l4_mapping[rf.branches[i]]);
				}

				// Populate atoms.
				for (size_t i = rf.begin; i < rf.end; ++i)
				{
					const atom& ra = l2.atoms[i];
					atoms.push_back(atom(ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.number, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
				}
				f.end = atoms.size();
				BOOST_ASSERT(f.begin < f.end);
			}
		}

		// Create new frames for ligand 2's frames that are neither f2 nor f2's parent frames.
		for (size_t k = l2_to_l4_mapping.front() + 1; k < l2_num_frames; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 2.
			const frame& rf = l2.frames[l4_to_l2_mapping[k]];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(f1_num_frames + l2_to_l4_mapping[rf.parent], l1.max_atom_number + rf.rotorX, l1.max_atom_number + rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				f.branches.push_back(f1_num_frames + l2_to_l4_mapping[rf.branches[i]]);
			}

			// Populate atoms.
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				const atom& ra = l2.atoms[i];
				atoms.push_back(atom(ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.number, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Create new frames for ligand 1's frames that are after f1.
		for (size_t k = f1_num_frames; k < l1_num_frames; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent > f1idx ? l2_num_frames + rf.parent : rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				f.branches.push_back(b > f1idx ? l2_num_frames + b : b);
			}

			// Populate atoms.
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			f.end = atoms.size();
		}

		BOOST_ASSERT(frames.size() == l1_num_frames + l2_num_frames);
		BOOST_ASSERT(frames.size() == frames.capacity());
		BOOST_ASSERT(atoms.size() == num_atoms);
		BOOST_ASSERT(atoms.size() == atoms.capacity());

		// The number of mutable atoms of child ligand is equal to the sum of its parent ligands minus 2.
		const size_t l1_num_mutatable_atoms = l1.mutable_atoms.size();
		const size_t l2_num_mutatable_atoms = l2.mutable_atoms.size();
		mutable_atoms.reserve(l1_num_mutatable_atoms + l2_num_mutatable_atoms - 2);

		// Copy the mutable atoms of ligand 1 except m1 to the child ligand.
		for (size_t i = 0; i < g1; ++i)
		{
			mutable_atoms.push_back(l1.mutable_atoms[i]);
		}
		for (size_t i = g1 + 1; i < l1_num_mutatable_atoms; ++i)
		{
			mutable_atoms.push_back(l1.mutable_atoms[i]);
		}

		// Copy the mutable atoms of ligand 2 except m2 to the child ligand.
		for (size_t i = 0; i < g2; ++i)
		{
			mutable_atoms.push_back(l1.max_atom_number + l2.mutable_atoms[i]);
		}
		for (size_t i = g2 + 1; i < l2_num_mutatable_atoms; ++i)
		{
			mutable_atoms.push_back(l1.max_atom_number + l2.mutable_atoms[i]);
		}
		BOOST_ASSERT(mutable_atoms.size() == l1_num_mutatable_atoms + l2_num_mutatable_atoms - 2);
		BOOST_ASSERT(mutable_atoms.size() == mutable_atoms.capacity());
	}

	ligand::ligand(const ligand& l1, const ligand& l2, const size_t f1idx, const size_t f2idx, const size_t g1, const size_t g2) : parent1(l1.p), parent2(l2.p)
	{
		const frame& f1 = l1.frames[f1idx];
		const frame& f2 = l2.frames[f2idx];

		// Set the connector atoms.
		connector1 = (g1 ? f1.rotorY : f1.rotorX);
		connector2 = (g2 ? f2.rotorY : f2.rotorX);
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
