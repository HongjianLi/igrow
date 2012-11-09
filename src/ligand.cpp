/*

   Copyright (c) 2012, The Chinese University of Hong Kong

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
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include "mat3.hpp"
#include "ligand.hpp"

namespace igrow
{
	using namespace boost::filesystem;

	ligand::ligand(const path& p) : p(p), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0)
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
		string line; // A line of ligand file in PDBQT format.
		line.reserve(79); // According to PDBQT specification, the last item AutoDock4 atom type locates at 1-based [78, 79].

		// Parse ATOM/HETATM, BRANCH, ENDBRANCH.
		ifstream in(p); // Parsing starts. Open the file stream as late as possible.
		while (getline(in, line) && !starts_with(line, "TORSDOF"))
		{
			++num_lines;
			if (starts_with(line, "ATOM") || starts_with(line, "HETATM"))
			{
				// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
				BOOST_ASSERT(current == frames.size() - 1);
				BOOST_ASSERT(f == &frames.back());

				// Validate the AutoDock4 atom type.
				const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
				const size_t ad = parse_ad_type_string(ad_type_string);
				if (ad == AD_TYPE_SIZE) throw parsing_error(p, num_lines, "Atom type " + ad_type_string + " is not supported by igrow.");

				// Parse the ATOM/HETATM line into an atom, which belongs to the current frame.
				string name = line.substr(12, 4);
				boost::algorithm::trim(name);
				atoms.push_back(atom(name, line.substr(12, 18), line.substr(54), right_cast<size_t>(line, 7, 11), vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54)), ad));

				// Update ligand properties.
				const atom& a = atoms.back();
				if (a.is_mutable()) mutable_atoms.push_back(a.srn);
				if (!a.is_hydrogen()) ++num_heavy_atoms;
				if (a.is_hb_donor()) ++num_hb_donors; // TODO: use neighbor.
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

				// The ending index of atoms of previous frame is the starting index of atoms of current frame.
				frames[current - 1].end = f->begin;

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
		frames.back().end = num_atoms;

		// Determine the number of rotatable bonds.
		num_rotatable_bonds = frames.size() - 1;
		BOOST_ASSERT(num_atoms + (num_rotatable_bonds << 1) + 3 <= num_lines); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines + REMARK lines (if any) == num_lines

		// Determine the maximum atom serial number.
		max_atom_number = atoms.back().srn;
		BOOST_ASSERT(max_atom_number >= num_atoms);
	}

	void ligand::save() const
	{
		ofstream out(p); // Dumping starts. Open the file stream as late as possible.
		using namespace std;
		out.setf(ios::fixed, ios::floatfield);
		out << setprecision(3);

		// Dump the ROOT frame.
		out << "ROOT\n";
		{
			const frame& f = frames.front();
			for (size_t i = f.begin; i < f.end; ++i)
			{
				const atom& a = atoms[i];
				out << "ATOM  " << setw(5) << a.srn << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << '\n';
			}
		}
		out << "ENDROOT\n";

		// Dump the BRANCH frames.
		vector<bool> dump_branches(frames.size()); // dump_branches[0] is dummy. The ROOT frame has been dumped.
		vector<size_t> stack; // Stack to track the depth-first traversal sequence of frames in order to avoid recursion.
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
					out << "ATOM  " << setw(5) << a.srn << ' ' << a.columns_13_to_30 << setw(8) << a.coordinate[0] << setw(8) << a.coordinate[1] << setw(8) << a.coordinate[2] << a.columns_55_to_79 << '\n';
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

	void ligand::update(const path& p)
	{
		if (!exists(p))
		{
			fe = 0;
			le = 0;
			return;
		}
		string line;
		line.reserve(79);
		ifstream in(p);
		getline(in, line); // MODEL        1
		getline(in, line); // REMARK       NORMALIZED FREE ENERGY PREDICTED BY IDOCK:  -4.976 KCAL/MOL
		fe = right_cast<fl>(line, 56, 63);
		getline(in, line); // REMARK            TOTAL FREE ENERGY PREDICTED BY IDOCK:  -6.722 KCAL/MOL
		getline(in, line); // REMARK     INTER-LIGAND FREE ENERGY PREDICTED BY IDOCK:  -7.740 KCAL/MOL
		getline(in, line); // REMARK     INTRA-LIGAND FREE ENERGY PREDICTED BY IDOCK:   1.018 KCAL/MOL
		getline(in, line); // REMARK            LIGAND EFFICIENCY PREDICTED BY IDOCK:  -0.280 KCAL/MOL
		le = right_cast<fl>(line, 56, 63);
		for (size_t i = 0; getline(in, line) && !starts_with(line, "TORSDOF");)
		{
			if (starts_with(line, "ATOM"))
			{
				BOOST_ASSERT(atoms[i].srn == right_cast<size_t>(line, 7, 11));
				atoms[i++].coordinate = vec3(right_cast<fl>(line, 31, 38), right_cast<fl>(line, 39, 46), right_cast<fl>(line, 47, 54));
			}
		}
		in.close();
		save();
	}

	std::pair<size_t, size_t> ligand::get_frame(const size_t srn) const
	{
		BOOST_ASSERT(num_rotatable_bonds == frames.size() - 1);
		for (size_t k = 0; k <= num_rotatable_bonds; ++k)
		{
			const frame& f = frames[k];
			const size_t srn_begin = atoms[f.begin].srn;
			const size_t srn_end = atoms[f.end - 1].srn;
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
					if (srn == atoms[i].srn) return std::pair<size_t,  size_t>(k, i);
				}
			}
		}
		throw std::domain_error("Failed to find an atom with serial number " + lexical_cast<string>(srn));
	}

	ligand::ligand(const path& p, const ligand& l1, const ligand& l2, const size_t g1, const size_t g2) : p(p), parent1(l1.p), parent2(l2.p)
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
		BOOST_ASSERT(c1idx < f1.end);
		BOOST_ASSERT(c2idx < f2.end);

		// Obtain constant references to the connector atoms.
		const atom& c1 = l1.atoms[c1idx];
		const atom& c2 = l2.atoms[c2idx];
		BOOST_ASSERT(f1idx == l1.get_frame(c1.srn).first);
		BOOST_ASSERT(f2idx == l2.get_frame(c2.srn).first);

		// Set the connector bonds.
		connector1 = lexical_cast<string>(c1.srn) + ":" + c1.name + " - " + lexical_cast<string>(m1.srn) + ":" + m1.name;
		connector2 = lexical_cast<string>(c2.srn) + ":" + c2.name + " - " + lexical_cast<string>(m2.srn) + ":" + m2.name;

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
			frames.push_back(frame(f1idx, c1.srn, l1.max_atom_number + c2.srn, atoms.size()));
			frame& f = frames.back();

			// Populate atoms.
			for (size_t i = f2.begin; i < m2idx; ++i)
			{
				const atom& ra = l2.atoms[i];
				atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
			}
			for (size_t i = m2idx + 1; i < f2.end; ++i)
			{
				const atom& ra = l2.atoms[i];
				atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
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
					atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
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
					atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
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
				atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
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

	ligand::ligand(const path& p, const ligand& l1, const size_t f1idx) : p(p), parent1(l1.p), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0)
	{
		const frame& f1 = l1.frames[f1idx];

		// The maximum atom serial number of child ligand is equal to the sum of its parent ligands.
		max_atom_number = l1.max_atom_number + 1;

		// Reserve enough capacity for storing frames.
		const size_t l1_num_frames = l1.frames.size();
		frames.reserve(l1_num_frames - 1);

		// Reserve enough capacity for storing atoms.
		atoms.reserve(l1.num_atoms - 1);

		// Determine the number of frames of ligand 5. Here, ligand 5 = ligand 1 - ligand 3.
		size_t child;
		for (child = f1idx; l1.frames[child].branches.size(); child = l1.frames[child].branches.back());
		BOOST_ASSERT(child < l1_num_frames);
		const size_t l5_num_frames = child - f1idx + 1;
		BOOST_ASSERT(l5_num_frames < l1_num_frames);

		// Create new frames for ligand 1's frames that are before f1's parent frame.
		for (size_t k = 0; k < f1.parent; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches);
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				f.branches.push_back(b > f1idx ? b - l5_num_frames : b);
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


		// Create a new frame for ligand 1's f1's parent frame.
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[f1.parent];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches);
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				if (b == f1idx) continue;
				f.branches.push_back(b > f1idx ? b - l5_num_frames : b);
			}

			// Populate atoms.
			BOOST_ASSERT(f.begin == rf.begin);
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}

			// Obtain the frames and indices of the two connector atoms.
			const std::pair<size_t, size_t> p1 = l1.get_frame(f1.rotorX);
			const std::pair<size_t, size_t> p2 = l1.get_frame(f1.rotorY);
			BOOST_ASSERT(p1.first == f1.parent);
			BOOST_ASSERT(p2.first == f1idx);

			// Obtain constant references to the connector atoms.
			const atom& c1 = l1.atoms[p1.second];
			const atom& m1 = l1.atoms[p2.second];
			BOOST_ASSERT(c1.srn == f1.rotorX);
			BOOST_ASSERT(m1.srn == f1.rotorY);

			// Set the connector bonds.
			connector1 = lexical_cast<string>(c1.srn) + ":" + c1.name + " - " + lexical_cast<string>(m1.srn) + ":" + m1.name;

			// Add a hydrogen.
			const vec3 c1_to_c2 = ((c1.covalent_radius() + ad_covalent_radii[0]) / (c1.covalent_radius() + m1.covalent_radius())) * (m1.coordinate - c1.coordinate); // Vector pointing from c1 to the new position of c2.
			const vec3 origin_to_c2 = c1.coordinate + c1_to_c2; // Translation vector to translate ligand 2 from origin to the new position of c2.
			atoms.push_back(atom("H", " H   <0> d        ", "  0.00  0.00     0.085 H ", max_atom_number, origin_to_c2, 0)); // c2 is a hydrogen.
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Create new frames for ligand 1's frames that are after f1's parent frame and before f1.
		for (size_t k = f1.parent + 1; k < f1idx; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches);
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				f.branches.push_back(b);
			}

			// Populate atoms.
			BOOST_ASSERT(f.begin == rf.begin + 1);
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Create new frames for ligand 1's frames that are after f1idx + l5_num_frames.
		for (size_t k = f1idx + l5_num_frames; k < l1_num_frames; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent > f1idx ? rf.parent - l5_num_frames : rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				BOOST_ASSERT(b > f1idx);
				f.branches.push_back(b - l5_num_frames);
			}

			// Populate atoms.
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Refresh the number of atoms.
		num_atoms = atoms.size();
		BOOST_ASSERT(num_atoms >= 1);
		BOOST_ASSERT(num_atoms < l1.num_atoms);

		// Refresh the number of rotatable bonds.
		num_rotatable_bonds = frames.size() - 1;
		BOOST_ASSERT(num_rotatable_bonds < l1.num_rotatable_bonds);

		// Refresh mutable_atoms, num_heavy_atoms, num_hb_donors, num_hb_acceptors and mw.
		mutable_atoms.reserve(num_atoms);
		for (const auto& a : atoms)
		{
			if (a.is_mutable()) mutable_atoms.push_back(a.srn);
			if (!a.is_hydrogen()) ++num_heavy_atoms;
			if (a.is_hb_donor()) ++num_hb_donors; // TODO: use neighbor.
			if (a.is_hb_acceptor()) ++num_hb_acceptors;
			mw += a.atomic_weight();
		}
		BOOST_ASSERT(mutable_atoms.size() <= l1.mutable_atoms.size() + 1);
	}

	ligand::ligand(const path& p, const ligand& l1, const ligand& l2, const size_t f1idx, const size_t f2idx, const bool dummy) : p(p), parent1(l1.p), parent2(l2.p), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0)
	{
		const frame& f1 = l1.frames[f1idx];
		const frame& f2 = l2.frames[f2idx];

		// The maximum atom serial number of child ligand is equal to the sum of its parent ligands.
		max_atom_number = l1.max_atom_number + l2.max_atom_number;

		// Reserve enough capacity for storing frames.
		const size_t l1_num_frames = l1.frames.size();
		const size_t l2_num_frames = l2.frames.size();
		frames.reserve(l1_num_frames + l2_num_frames);

		// Reserve enough capacity for storing atoms.
		atoms.reserve(l1.num_atoms + l2.num_atoms);

		// Determine the number of frames of ligand 5 and ligand 4. Here, ligand 5 = ligand 1 - ligand 3.
		size_t child;
		for (child = f1idx; l1.frames[child].branches.size(); child = l1.frames[child].branches.back());
		BOOST_ASSERT(child < l1_num_frames);
		const size_t l5_num_frames = child - f1idx + 1;
		BOOST_ASSERT(l5_num_frames < l1_num_frames);
		for (child = f2idx; l2.frames[child].branches.size(); child = l2.frames[child].branches.back());
		BOOST_ASSERT(child < l2_num_frames);
		const size_t l4_num_frames = child - f2idx + 1;
		BOOST_ASSERT(l4_num_frames <= l2_num_frames);

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
				f.branches.push_back(b > f1idx ? l4_num_frames + b - l5_num_frames : b);
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

		// Obtain the frames and indices of the two connector atoms.
		const std::pair<size_t, size_t> p1 = l1.get_frame(f1.rotorX);
		const std::pair<size_t, size_t> p2 = l2.get_frame(f2.rotorY);
		BOOST_ASSERT(p1.first == f1.parent);
		BOOST_ASSERT(p2.first == f2idx);

		// Obtain constant references to the connector atoms.
		const atom& c1 = l1.atoms[p1.second];
		const atom& c2 = l2.atoms[p2.second];
		BOOST_ASSERT(c1.srn == f1.rotorX);
		BOOST_ASSERT(c2.srn == f2.rotorY);

		// Obtain the frames and indices of the two virtual mutable atoms.
		const std::pair<size_t, size_t> q1 = l1.get_frame(f1.rotorY);
		const std::pair<size_t, size_t> q2 = l2.get_frame(f2.rotorX);
		BOOST_ASSERT(q1.first == f1idx);
		BOOST_ASSERT(q2.first == f2.parent);

		// Obtain constant references to the virtual mutable atoms.
		const atom& m1 = l1.atoms[q1.second];
		const atom& m2 = l2.atoms[q2.second];
		BOOST_ASSERT(m1.srn == f1.rotorY);
		BOOST_ASSERT(m2.srn == f2.rotorX);

		// Set the connector bonds.
		connector1 = lexical_cast<string>(c1.srn) + ":" + c1.name + " - " + lexical_cast<string>(m1.srn) + ":" + m1.name;
		connector2 = lexical_cast<string>(c2.srn) + ":" + c2.name + " - " + lexical_cast<string>(m2.srn) + ":" + m2.name;

		// Calculate the translation vector for moving ligand 2 to a nearby place of ligand 1.
		const vec3 c1_to_c2 = ((c1.covalent_radius() + c2.covalent_radius()) / (c1.covalent_radius() + m1.covalent_radius())) * (m1.coordinate - c1.coordinate); // Vector pointing from c1 to the new position of c2.
		const vec3 origin_to_c2 = c1.coordinate + c1_to_c2; // Translation vector to translate ligand 2 from origin to the new position of c2.
		const vec3 c2_to_c1_nd = (-1 * c1_to_c2).normalize(); // Normalized vector pointing from c2 to c1.
		const vec3 c2_to_m2_nd = (m2.coordinate - c2.coordinate).normalize(); // Normalized vector pointing from c2 to m2.
		const mat3 rot(cross_product(c2_to_m2_nd, c2_to_c1_nd).normalize(), c2_to_m2_nd * c2_to_c1_nd); // Rotation matrix to rotate m2 along the normal to the direction from the new position of c2 to c1.

		// Create new frames for ligand 2's frames that are either f2 or f2's child frames.
		for (size_t k = 0; k < l4_num_frames; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 2.
			const frame& rf = l2.frames[f2idx + k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(k ? f1idx + rf.parent - f2idx : f1.parent, k ? l1.max_atom_number + rf.rotorX : f1.rotorX, l1.max_atom_number + rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				f.branches.push_back(f1idx + rf.branches[i] - f2idx);
			}

			// Populate atoms.
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				const atom& ra = l2.atoms[i];
				atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coordinate - c2.coordinate) + origin_to_c2, ra.ad));
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Create new frames for ligand 1's frames that are after f1idx + l5_num_frames.
		for (size_t k = f1idx + l5_num_frames; k < l1_num_frames; ++k)
		{
			// Obtain a constant reference to the corresponding frame of ligand 1.
			const frame& rf = l1.frames[k];
			const size_t rf_num_branches = rf.branches.size();

			// Create a new frame based on the reference frame.
			frames.push_back(frame(rf.parent > f1idx ? l4_num_frames + rf.parent - l5_num_frames : rf.parent, rf.rotorX, rf.rotorY, atoms.size()));
			frame& f = frames.back();

			// Populate branches.
			f.branches.reserve(rf_num_branches); // This frame exactly consists of rf_num_branches BRANCH frames.
			for (size_t i = 0; i < rf_num_branches; ++i)
			{
				const size_t b = rf.branches[i];
				f.branches.push_back(b > f1idx ? l4_num_frames + b - l5_num_frames : b);
			}

			// Populate atoms.
			for (size_t i = rf.begin; i < rf.end; ++i)
			{
				atoms.push_back(l1.atoms[i]);
			}
			f.end = atoms.size();
			BOOST_ASSERT(f.begin < f.end);
		}

		// Refresh the number of atoms.
		num_atoms = atoms.size();
		BOOST_ASSERT(num_atoms >= 1);
		BOOST_ASSERT(num_atoms < l1.num_atoms + l2.num_atoms);

		// Refresh the number of rotatable bonds.
		num_rotatable_bonds = frames.size() - 1;
		BOOST_ASSERT(num_rotatable_bonds >= 1);
		BOOST_ASSERT(num_rotatable_bonds <= l1.num_rotatable_bonds + l2.num_rotatable_bonds - 1);

		// Refresh mutable_atoms, num_heavy_atoms, num_hb_donors, num_hb_acceptors and mw.
		mutable_atoms.reserve(num_atoms);
		for (const auto& a : atoms)
		{
			if (a.is_mutable()) mutable_atoms.push_back(a.srn);
			if (!a.is_hydrogen()) ++num_heavy_atoms;
			if (a.is_hb_donor()) ++num_hb_donors; // TODO: use neighbor.
			if (a.is_hb_acceptor()) ++num_hb_acceptors;
			mw += a.atomic_weight();
		}
		BOOST_ASSERT(mutable_atoms.size() <= l1.mutable_atoms.size() + l2.mutable_atoms.size());
	}
}
