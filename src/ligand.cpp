#include <iomanip>
#include <boost/filesystem/fstream.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>
#include "array.hpp"
#include "ligand.hpp"
using namespace boost;
using namespace boost::filesystem;

ligand::ligand(const path& p) : p(p), num_hb_donors(0), num_hb_acceptors(0), ma(0)
{
	// Initialize necessary variables for constructing a ligand.
	frames.reserve(30); // A ligand typically consists of <= 30 frames.
	frames.push_back(frame(0, 0, 0, 0)); // ROOT is also treated as a frame. The parent, rotorX, and rotorY of ROOT frame are dummy.
	frames.back().branches.reserve(4); // A frame typically consists of <= 4 BRANCH frames.

	// Initialize helper variables for parsing.
	size_t current = 0; // Index of current frame, initialized to ROOT frame.
	frame* f = &frames.front(); // Pointer to the current frame.
	size_t num_lines = 0; // Used to track line number for reporting parsing errors, if any.
	string line; // A line of ligand file in PDBQT format.

	// Parse ATOM/HETATM, BRANCH, ENDBRANCH.
	for (boost::filesystem::ifstream ifs(p); getline(ifs, line);)
	{
		++num_lines;
		const string record = line.substr(0, 6);
		if (record == "TORSDO") break;
		if (record == "ATOM  " || record == "HETATM")
		{
			// Whenever an ATOM/HETATM line shows up, the current frame must be the last one.
			assert(current == frames.size() - 1);
			assert(f == &frames.back());

			// Validate the AutoDock4 atom type.
			const string ad_type_string = line.substr(77, isspace(line[78]) ? 1 : 2);
			const size_t ad = atom::parse_ad_string(ad_type_string);

			// Parse the ATOM/HETATM line into an atom, which belongs to the current frame.
			string name = line.substr(12, 4);
			boost::algorithm::trim(name);
			atoms.push_back(atom(name, line.substr(12, 18), line.substr(54), stoul(line.substr(6, 5)), {stod(line.substr(30, 8)), stod(line.substr(38, 8)), stod(line.substr(46, 8))}, ad));

			// Update ligand properties.
			const atom& a = atoms.back();
			if (a.is_hb_donor()) ++num_hb_donors;
			if (a.is_hb_acceptor()) ++num_hb_acceptors;
			ma += a.atomic_mass();
		}
		else if (record == "BRANCH")
		{
			// Parse "BRANCH   X   Y". X and Y are right-justified and 4 characters wide.
			frames.push_back(frame(current, stoul(line.substr(6, 4)), stoul(line.substr(10, 4)), atoms.size()));

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
		else if (record == "ENDBRA")
		{
			// A frame may be empty, e.g. "BRANCH   4   9" is immediately followed by "ENDBRANCH   4   9".
			// This emptiness is likely to be caused by invalid input structure, especially when all the atoms are located in the same plane.
			if (f->begin == atoms.size()) throw domain_error("Error parsing " + p.filename().string() + ": an empty BRANCH has been detected, indicating the input ligand structure is probably invalid.");

			// Now the parent of the following frame is the parent of current frame.
			current = f->parent;

			// Update the pointer to the current frame.
			f = &frames[current];
		}
	}

	assert(current == 0); // current should remain its original value if "BRANCH" and "ENDBRANCH" properly match each other.
	assert(f == &frames.front()); // The frame pointer should point to the ROOT frame.

	// Determine the number of atoms.
	num_atoms = atoms.size();
	frames.back().end = num_atoms;

	// Determine the number of rotatable bonds.
	num_rotatable_bonds = frames.size() - 1;
	assert(num_atoms + (num_rotatable_bonds << 1) + 3 <= num_lines); // ATOM/HETATM lines + BRANCH/ENDBRANCH lines + ROOT/ENDROOT/TORSDOF lines + REMARK lines (if any) == num_lines

	// Determine the maximum atom serial number.
	max_atom_number = atoms.back().srn;
	assert(max_atom_number >= num_atoms);
}

bool ligand::crossover_feasible() const
{
	return num_rotatable_bonds > 0;
}

bool ligand::operator<(const ligand& l) const
{
	return fe < l.fe;
}

void ligand::save() const
{
	boost::filesystem::ofstream ofs(p);
	ofs.setf(ios::fixed, ios::floatfield);
	ofs << setprecision(3);

	// Dump the ROOT frame.
	ofs << "ROOT\n";
	{
		const frame& f = frames.front();
		for (size_t i = f.begin; i < f.end; ++i)
		{
			const atom& a = atoms[i];
			ofs << "ATOM  " << setw(5) << a.srn << ' ' << a.columns_13_to_30 << setw(8) << a.coord[0] << setw(8) << a.coord[1] << setw(8) << a.coord[2] << a.columns_55_to_79 << '\n';
		}
	}
	ofs << "ENDROOT\n";

	// Dump the BRANCH frames.
	vector<bool> branches_written(frames.size()); // branches_written[0] is dummy. The ROOT frame has been written.
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
		if (!branches_written[fn]) // This BRANCH frame has not been written.
		{
			ofs << "BRANCH"    << setw(4) << f.rotorX << setw(4) << f.rotorY << '\n';
			for (size_t i = f.begin; i < f.end; ++i)
			{
				const atom& a = atoms[i];
				ofs << "ATOM  " << setw(5) << a.srn << ' ' << a.columns_13_to_30 << setw(8) << a.coord[0] << setw(8) << a.coord[1] << setw(8) << a.coord[2] << a.columns_55_to_79 << '\n';
			}
			branches_written[fn] = true;
			for (auto i = f.branches.rbegin(); i < f.branches.rend(); ++i)
			{
				stack.push_back(*i);
			}
		}
		else // This BRANCH frame has been written.
		{
			ofs << "ENDBRANCH" << setw(4) << f.rotorX << setw(4) << f.rotorY << '\n';
			stack.pop_back();
		}
	}
	ofs << "TORSDOF " << num_rotatable_bonds << '\n';
}

void ligand::update(const path& p)
{
	if (!exists(p))
	{
		fe = 0;
		return;
	}
	string line;
	boost::filesystem::ifstream ifs(p);
	getline(ifs, line); // MODEL        1
	getline(ifs, line); // REMARK       NORMALIZED FREE ENERGY PREDICTED BY IDOCK:  -4.976 KCAL/MOL
	fe = stod(line.substr(55, 8));
	getline(ifs, line); // REMARK            TOTAL FREE ENERGY PREDICTED BY IDOCK:  -6.722 KCAL/MOL
	getline(ifs, line); // REMARK     INTER-LIGAND FREE ENERGY PREDICTED BY IDOCK:  -7.740 KCAL/MOL
	getline(ifs, line); // REMARK     INTRA-LIGAND FREE ENERGY PREDICTED BY IDOCK:   1.018 KCAL/MOL
	getline(ifs, line); // REMARK    RF-SCORE BINDING AFFINITY PREDICTED BY IDOCK:   6.532 PKD
	for (size_t i = 0; getline(ifs, line);)
	{
		const string record = line.substr(0, 6);
		if (record == "TORSDO") break;
		if (record == "ATOM  ")
		{
			assert(atoms[i].srn == stoul(line.substr(6, 5)));
			atoms[i++].coord = {stod(line.substr(30, 8)), stod(line.substr(38, 8)), stod(line.substr(46, 8))};
		}
	}
	ifs.close();
	save();
}

pair<size_t, size_t> ligand::get_frame(const size_t srn) const
{
	assert(num_rotatable_bonds == frames.size() - 1);
	for (size_t k = 0; k <= num_rotatable_bonds; ++k)
	{
		const frame& f = frames[k];
		const size_t srn_begin = atoms[f.begin].srn;
		const size_t srn_end = atoms[f.end - 1].srn;
		assert(srn_begin <= srn_end);
		if (f.end - f.begin == srn_end - srn_begin + 1) // The serial numbers are continuous, which is the most cases.
		{
			if (srn_begin <= srn && srn <= srn_end) return pair<size_t, size_t>(k, f.begin + srn - srn_begin);
		}
		else // The serial numbers are not continuous, but they are sorted.
		{
			// Linear search at the moment. Binary search can be used.
			for (size_t i = f.begin; i < f.end; ++i)
			{
				if (srn == atoms[i].srn) return pair<size_t,  size_t>(k, i);
			}
		}
	}
	throw domain_error("Failed to find an atom with serial number " + to_string(srn));
}

ligand::ligand(const path& p, const ligand& l1, const ligand& l2, const size_t f1idx, const size_t f2idx) : p(p), parent1(l1.p), parent2(l2.p), num_hb_donors(0), num_hb_acceptors(0), ma(0)
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
	assert(child < l1_num_frames);
	const size_t l5_num_frames = child - f1idx + 1;
	assert(l5_num_frames < l1_num_frames);
	for (child = f2idx; l2.frames[child].branches.size(); child = l2.frames[child].branches.back());
	assert(child < l2_num_frames);
	const size_t l4_num_frames = child - f2idx + 1;
	assert(l4_num_frames <= l2_num_frames);

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
		assert(f.begin == rf.begin);
		for (size_t i = rf.begin; i < rf.end; ++i)
		{
			atoms.push_back(l1.atoms[i]);
		}
		f.end = atoms.size();
		assert(f.begin < f.end);
	}

	// Obtain the frames and indices of the two connector atoms.
	const pair<size_t, size_t> p1 = l1.get_frame(f1.rotorX);
	const pair<size_t, size_t> p2 = l2.get_frame(f2.rotorY);
	assert(p1.first == f1.parent);
	assert(p2.first == f2idx);

	// Obtain constant references to the connector atoms.
	const atom& c1 = l1.atoms[p1.second];
	const atom& c2 = l2.atoms[p2.second];
	assert(c1.srn == f1.rotorX);
	assert(c2.srn == f2.rotorY);

	// Obtain the frames and indices of the two virtual mutable atoms.
	const pair<size_t, size_t> q1 = l1.get_frame(f1.rotorY);
	const pair<size_t, size_t> q2 = l2.get_frame(f2.rotorX);
	assert(q1.first == f1idx);
	assert(q2.first == f2.parent);

	// Obtain constant references to the virtual mutable atoms.
	const atom& m1 = l1.atoms[q1.second];
	const atom& m2 = l2.atoms[q2.second];
	assert(m1.srn == f1.rotorY);
	assert(m2.srn == f2.rotorX);

	// Set the connector bonds.
	connector1 = to_string(c1.srn) + ":" + c1.name + " - " + to_string(m1.srn) + ":" + m1.name;
	connector2 = to_string(c2.srn) + ":" + c2.name + " - " + to_string(m2.srn) + ":" + m2.name;

	// Calculate the translation vector for moving ligand 2 to a nearby place of ligand 1.
	const array<double, 3> c1_to_c2 = ((c1.covalent_radius() + c2.covalent_radius()) / (c1.covalent_radius() + m1.covalent_radius())) * (m1.coord - c1.coord); // Vector pointing from c1 to the new position of c2.
	const array<double, 3> origin_to_c2 = c1.coord + c1_to_c2; // Translation vector to translate ligand 2 from origin to the new position of c2.
	const array<double, 3> c2_to_c1_nd = normalize(-1 * c1_to_c2); // Normalized vector pointing from c2 to c1.
	const array<double, 3> c2_to_m2_nd = normalize(m2.coord - c2.coord); // Normalized vector pointing from c2 to m2.
	const array<double, 9> rot = vec3_to_mat3(normalize(cross_product(c2_to_m2_nd, c2_to_c1_nd)), c2_to_m2_nd * c2_to_c1_nd); // Rotation matrix to rotate m2 along the normal to the direction from the new position of c2 to c1.

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
			atoms.push_back(atom(ra.name, ra.columns_13_to_30, ra.columns_55_to_79, l1.max_atom_number + ra.srn, rot * (ra.coord - c2.coord) + origin_to_c2, ra.ad));
		}
		f.end = atoms.size();
		assert(f.begin < f.end);
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
		assert(f.begin < f.end);
	}

	// Refresh the number of atoms.
	num_atoms = atoms.size();
	assert(num_atoms >= 1);
	assert(num_atoms < l1.num_atoms + l2.num_atoms);

	// Refresh the number of rotatable bonds.
	num_rotatable_bonds = frames.size() - 1;
	assert(num_rotatable_bonds >= 1);
	assert(num_rotatable_bonds <= l1.num_rotatable_bonds + l2.num_rotatable_bonds - 1);

	// Refresh num_hb_donors, num_hb_acceptors and mw.
	for (const auto& a : atoms)
	{
		if (a.is_hb_donor()) ++num_hb_donors;
		if (a.is_hb_acceptor()) ++num_hb_acceptors;
		ma += a.atomic_mass();
	}
}
