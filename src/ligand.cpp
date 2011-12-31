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

	using namespace boost;

	ligand::ligand(const path& p) : p(p), num_heavy_atoms(0), num_hb_donors(0), num_hb_acceptors(0), mw(0), connector1(0), connector2(0), logp(0) // TODO: comment logp(0)
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
				
				// Find the neighbors of a in the same frame.
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
						
						// Adjust num_hb_donors because AutoDock4 does not include explicit atom types for hydrogen bond donors.
						if (a.is_polar_hydrogen()) // b is a hydrogen bond donor.
						{
							// Increment num_hb_donors if b has not been counted.
							if (!f.atoms[b.neighbors.back().index].is_polar_hydrogen()) ++num_hb_donors;
							BOOST_ASSERT(a.neighbors.size() == 1);
							break; // A hydrogen can only have one neighbor.
						}
					}
				}
				
				if (a.is_mutable()) mutable_atoms.push_back(idx);
				if (!a.is_hydrogen()) ++ num_heavy_atoms;				
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
		BOOST_ASSERT(num_heavy_atoms + 2 < num_lines); // ROOT, ENDROOT

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
			const size_t num_branches = f.branches.size();
			for (size_t i = 0; i < num_branches; ++i)
			{
				stack.push_back(f.branches[i]);
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
		out << "TORSDOF " << (num_frames - 1);
		out.close();
	}

	ligand* ligand::mutate(const ligand& other, const mt19937eng& eng) const
	{
		using boost::random::variate_generator;
		using boost::random::uniform_int_distribution;
		variate_generator<mt19937eng, uniform_int_distribution<size_t> > uniform_mutation_point_gen_1(eng, uniform_int_distribution<size_t>(0, this->mutable_atoms.size() - 1));
		variate_generator<mt19937eng, uniform_int_distribution<size_t> > uniform_mutation_point_gen_2(eng, uniform_int_distribution<size_t>(0, other.mutable_atoms.size() - 1));
		const atom_index& mp1 = this->mutable_atoms[uniform_mutation_point_gen_1()];
		const atom_index& mp2 = other.mutable_atoms[uniform_mutation_point_gen_2()];

		const frame& f1 = this->frames[mp1.frame];
		const frame& f2 = other.frames[mp2.frame];
		const atom& p1 = f1.atoms[mp1.index]; // Mutation point 1
		const atom& p2 = f2.atoms[mp2.index]; // Mutation point 2
		BOOST_ASSERT(p1.neighbors.size() == 1); // The mutation point should consists of only one non-rotatable single bond.
		BOOST_ASSERT(p2.neighbors.size() == 1); // The mutation point should consists of only one non-rotatable single bond.
		BOOST_ASSERT(p1.neighbors.front().frame == mp1.frame); // Both the connector and mutation point should be in the same frame.
		BOOST_ASSERT(p2.neighbors.front().frame == mp2.frame); // Both the connector and mutation point should be in the same frame.
		const atom& c1 = f1.atoms[p1.neighbors.front().index]; // Connector 1
		const atom& c2 = f2.atoms[p2.neighbors.front().index]; // Connector 2

		ligand child;
		child.parent1 = this->p;
		child.parent2 = other.p;
		child.connector1 = c1.number;
		child.connector2 = c2.number;
		child.frames.reserve(this->frames.size() + other.frames.size());
		child.mutable_atoms.reserve(this->mutable_atoms.size() + other.mutable_atoms.size());
		child.num_heavy_atoms = this->num_heavy_atoms + other.num_heavy_atoms - ((p1.is_hydrogen() ? 0 : 1) + (p2.is_hydrogen() ? 0 : 1));
		child.num_hb_donors = this->num_hb_donors + other.num_hb_donors;
		child.num_hb_acceptors = this->num_hb_acceptors + other.num_hb_acceptors;
		child.mw = this->mw + other.mw - (p1.atomic_weight() + p2.atomic_weight());
		child.logp = this->logp + other.logp; // TODO: comment this line.

		// Check ligand validity, i.e. steric clash, rule of 5

		return new ligand(other);
	}

	void ligand::evaluate_efficacy()
	{
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
