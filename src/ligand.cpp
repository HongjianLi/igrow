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

#include "common.hpp"
#include "ligand.hpp"
#include "mat.hpp"
#include <cfloat>
#include <sstream>
#include <fstream>

const double M_PI = 3.14159265358979323846;

using namespace std;

// default to create bond when not found

Ligand::Ligand()
{
    SkipConnectGenerate = false;
}

Ligand* Ligand::split(const Ligand& ref)
{
    Ligand ref_copy(ref);
    Ligand* child = new Ligand();

    // to detect whether it is a ring structure
    list<int> ring;
    int count = DetectRing(ring);
    // clean overlapping and traversed set
    overlap.clear();
    scanned.clear();
    // replace with hydrogen only
    toAddAtoms.clear();
    map<int, atom>::iterator it1, it2;
    int indexToAdd;
    atom atomH;

    // loop through all atoms with the reference, assign overlaps by same type and position
    for (it1 = atoms.begin(); it1 != atoms.end(); ++it1)
        for (it2 = ref_copy.atoms.begin(); it2 != ref_copy.atoms.end(); ++it2)
            if ((it1->second.element == it2->second.element) && (it1->second.coordinates == it2->second.coordinates))
                overlap.insert(it1->first);

    // randomly get some fragments
    scan_recursive(MinIndex());

    // make a copy for the child
    for (it1 = atoms.begin(); it1 != atoms.end(); ++it1)
        child->atoms.insert(pair<int, atom > (it1->first, atom(it1->second)));

    // start from the end
    map<int, atom>::reverse_iterator r_iter = atoms.rbegin();
    set<int> toAddchild, toAddparent;
    toAddchild.clear();
    toAddparent.clear();
    while (r_iter != atoms.rend()) {
        // atoms is being rejected, see if replacement is needed
        if (scanned.find(r_iter->first) == scanned.end()) {
            // replacement found
            if (toAddAtoms.find(r_iter->first) != toAddAtoms.end())
                toAddchild.insert(r_iter->first);
            // delete the atom from the child
            child->DeleteAtom(r_iter->first);
            /*if (toAddAtoms.find(pendingIndex) != toAddAtoms.end()) {
                    // retrieve corresponding hydrogen atom
                    atomH = toAddAtoms[pendingIndex];
                    child->atoms.insert(pair<int,Atom>(pendingIndex,atomH));
                    // get its connected atom index
                    indexToAdd = *(atomH.IndexArray.begin());
                    // update connection
                    child->atoms[indexToAdd].IndexArray.insert(pendingIndex);
            }*/
            // atoms to be deleted in this group
        }
        else if (overlap.find(r_iter->first) == overlap.end()) {
            // replacement found
            if (toAddAtoms.find(r_iter->first) != toAddAtoms.end())
                toAddparent.insert(r_iter->first);
            // delete atom from this molecule
            r_iter = DeleteAtom(r_iter->first);
            /*if (toAddAtoms.find(pendingIndex) != toAddAtoms.end()) {
                    atomH = toAddAtoms[pendingIndex];
                    atoms.insert(pair<int,Atom>(pendingIndex,atomH));
                    indexToAdd = *(atomH.IndexArray.begin());
                    atoms[indexToAdd].IndexArray.insert(pendingIndex);
            }*/
            // iterator has already incremented in the removal process
            continue;
        }
        ++r_iter;
    }
    for (set<int>::iterator set_it = toAddchild.begin(); set_it != toAddchild.end(); ++set_it) {
        // retrieve corresponding hydrogen atom
        atomH = toAddAtoms[*set_it];
        child->atoms.insert(pair<int, atom > (*set_it, atomH));
        // get its connected atom index
        indexToAdd = *(atomH.IndexArray.begin());
        // update connection
        child->atoms[indexToAdd].IndexArray.insert(*set_it);
    }
    for (set<int>::iterator set_it = toAddparent.begin(); set_it != toAddparent.end(); ++set_it) {
        // retrieve corresponding hydrogen atom
        atomH = toAddAtoms[*set_it];
        atoms.insert(pair<int, atom > (*set_it, atomH));
        // get its connected atom index
        indexToAdd = *(atomH.IndexArray.begin());
        // update connection
        atoms[indexToAdd].IndexArray.insert(*set_it);
    }

#ifdef OBJECT_SPACE
    RotateLine(mole_axis, mole_angle);
    //Translate(translation);
    //child->translation = Vec3d(ref.translation);
    child->mole_axis = Vec3d(ref.mole_axis);
    child->mole_angle = ref.mole_angle;
    child->RotateLine(child->mole_axis, child->mole_angle);
    //child->Translate(child->translation);
#endif

    return child;
}

// add a fragment to the position of a randomly selected hydrogen

int Ligand::mutate(string FilenameOfFragment)
{
    // local molecule copy
    Ligand fragment;
    fragment.LoadPDB(FilenameOfFragment);

    // select hydrogen
    int count(0), connectIndex, fragHydrogen(-1), fragIndex, linkerHydrogen(-1);
    // atom copy and atom reference
    atom curHydrogen, connectAtom, *fragHydrogenAtom, *fragConnectAtom;
    string element1, element2;
    Vec3d delta;

    linkerHydrogen = IndexOfRandomHydrogen();
    // obtain information on hydrogens
    curHydrogen = atoms[linkerHydrogen];
    // hydrogen has only 1 connection, get the first index in array
    connectIndex = *(curHydrogen.IndexArray.begin());
    connectAtom = atoms[connectIndex];
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
    if (connectAtom.isSP2() && fragConnectAtom->isSP2()) {
        int cur_dihedral(-1), frag_dihedral(-1);
        double dihedral, dist;
        set<int>::iterator it;
        map<double, double> best_angles;
        // find reference atoms for both connected atoms
        for (it = connectAtom.IndexArray.begin(); it != connectAtom.IndexArray.end(); ++it) {
            if (*it != linkerHydrogen)
                cur_dihedral = *it;
        }
        for (it = fragConnectAtom->IndexArray.begin(); it != fragConnectAtom->IndexArray.end(); ++it) {
            if (*it != fragHydrogen)
                frag_dihedral = *it;
        }
        // calculate dihedral angle on these four atoms
        dihedral = DihedralAngle(atoms[cur_dihedral], connectAtom, *fragConnectAtom, fragment.atoms[frag_dihedral]);
        // rotate so that the two planes are now parallel
        fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, -dihedral);
        // calculate distance in this orientation
        dist = MolecularDistance(fragment);
        best_angles.insert(pair<double, double>(dist, -dihedral));
        // try the opposite
        fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, M_PI);
        dist = MolecularDistance(fragment);
        best_angles.insert(pair<double, double>(dist, M_PI - dihedral));
        // try the other set of orientations
        fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, -M_PI + 2 * dihedral);
        dist = MolecularDistance(fragment);
        best_angles.insert(pair<double, double>(dist, dihedral));
        fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, M_PI);
        dist = MolecularDistance(fragment);
        best_angles.insert(pair<double, double>(dist, M_PI + dihedral));
        fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, -M_PI - dihedral);
        fragment.RotateLine(connectAtom.coordinates, fragConnectAtom->coordinates, fragIndex, best_angles.rbegin()->second); // maximize intra-molecular distance
    }
        // minimise hinderance for SP3 and SP bond
    else {
        // Rotating fragment around connecting bond to minimize steric hindrance...
        double BadContact, BestContact(0), BestAngle(0), Angle(0), delta(M_PI / 25);
        Vec3d normal = connectAtom.coordinates - fragConnectAtom->coordinates;
        // try a whole circle with definite steps
        while (Angle < 2 * M_PI) {
            Angle += delta;
            if (normal == Vec3d()) break; // impossible, useless
            fragment.RotateLine(normal, fragIndex, delta);
            BadContact = MolecularDistance(fragment);
            // maximize the distance
            if (BadContact > BestContact) {
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

    for (map<int, atom>::iterator it = fragment.atoms.begin(); it != fragment.atoms.end(); ++it) {
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

    return 0;
}

int Ligand::IndexOfRandomHydrogen()
{
    set<int> index_of_hydrogen;
    set<int>::iterator it;
    // find out all hydrogen atoms in the molecule
    for (map<int, atom>::iterator iter = atoms.begin(); iter != atoms.end(); ++iter) {
        if (iter->second.element.at(0) == 'H' && iter->second.element.size() == 1)
            index_of_hydrogen.insert(iter->first);
    }
    // if there is no more hydrogen atom, return error
    if (index_of_hydrogen.empty()) return -1;
    // randomly choose an index
//    int index = generator_int() % index_of_hydrogen.size();
    int index = 0 % index_of_hydrogen.size();    
    it = index_of_hydrogen.begin();
    advance(it, index);
    // return index of chosen hydrogen
    return *it;
}

void Ligand::Translate(int index, Vec3d origin)
{
    // calculate real displacement to move
    Vec3d delta = origin - atoms[index].coordinates;
    // move all atoms in the molecule
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
        it->second.coordinates = it->second.coordinates + delta;
}

void Ligand::Translate(Vec3d origin)
{
    // determine centre of gravity
    Vec3d centre = CentreOfGravity();
    Vec3d delta = origin - centre;
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
        it->second.coordinates = it->second.coordinates + delta;
}

void Ligand::RotateLine(Vec3d v1, Vec3d v2, int index, double radian)
{
    Vec3d location, delta = atoms[index].coordinates;
    Mat4d rot;
    // create rotation matrix
    rot = rot.createRotation(radian, v1 - v2);
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        // shift to indexed atom position
        location = it->second.coordinates - delta;
        // perform rotation
        location = rot * location;
        // move back by delta amount
        it->second.coordinates = location + delta;
    }
}

// same method by calculating normal beforehand

void Ligand::RotateLine(Vec3d normal, int index, double radian)
{
    Vec3d location, delta = atoms[index].coordinates;
    Mat4d rot;
    rot = rot.createRotation(radian, normal);
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        location = it->second.coordinates - delta;
        location = rot * location;
        it->second.coordinates = location + delta;
    }
}

void Ligand::RotateLine(Vec3d normal, double radian)
{
    Vec3d location, delta = CentreOfGravity();
    Mat4d rot;
    rot = rot.createRotation(radian, normal);
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        location = it->second.coordinates - delta;
        location = rot * location;
        it->second.coordinates = location + delta;
    }
}

// centre of rotation given by indexed atom

void Ligand::Rotate(int index, Vec3d EulerAngles)
{
    Vec3d centre = atoms[index].coordinates;
    Rotate(centre, EulerAngles);
}

void Ligand::Rotate(Vec3d EulerAngles)
{
    // determine centre of gravity
    Vec3d centre = CentreOfGravity();
    Rotate(centre, EulerAngles);
}

int Ligand::SavePDB(string Filename)
{
    // try to open file, always overwrite
    ofstream out_file;
    out_file.open(Filename.c_str(), ios::out);
    if (!out_file) return -1;

    string connectData;
    bool hasHeader = false;
    // save comments
    for (list<string>::iterator it = comments.begin(); it != comments.end(); ++it) {
        if (!hasHeader) {
            out_file << "HEADER    " << *it << endl;
            hasHeader = true;
        }
        else {
            out_file << "REMARK    " << *it << endl;
        }
    }

    // save coordinates
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
        out_file << it->second.WritePDBLine(it->first) << endl;

    // connect data
    ostringstream output;
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        // append index
        output.str(string());
        output.fill(' ');
        output.width(5);
        output << right << it->first;
        connectData = "CONECT" + output.str();
        // append all connection data
        for (set<int>::iterator iter = it->second.IndexArray.begin(); iter != it->second.IndexArray.end(); ++iter) {
            output.str(string());
            output.fill(' ');
            output.width(5);
            output << *iter;
            connectData += output.str();
        }
        // fill the line to 82 characters
        if (connectData.length() < 82) {
            output.str(string());
            output.fill(' ');
            output.width(82);
            output << left << connectData;
        }
        // write to file where there is something
        if (!(it->second.IndexArray.empty()))
            out_file << output.str();
        out_file << endl;
    }
    /*// ID data
    connectData = "REMARK     ID ";
    for (map<int,Atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
            output.str(string());
            output << left << it->first << " " << it->second.ID;
            out_file << connectData << output.str();
            out_file << endl;
    }*/
#ifdef OBJECT_SPACE
    // orientation and translation
    out_file << "REMARK     ROTX " << mole_axis.n[0] << endl;
    out_file << "REMARK     ROTY " << mole_axis.n[1] << endl;
    out_file << "REMARK     ROTZ " << mole_axis.n[2] << endl;
    out_file << "REMARK     ANGLE " << mole_angle << endl;
#endif

    out_file.close();
    return 0;
}

#define LINE_LENGTH	512

int Ligand::LoadPDB(string Filename)
{
    // reset container
    comments.clear();
    scanned.clear();
    overlap.clear();
    atoms.clear();
    toAddAtoms.clear();

    // produce an ID for this molecule
//    unsigned int ID = generator_int();
    unsigned int ID = 0;

    ifstream in_file;
    in_file.open(Filename.c_str());
    if (!in_file) return -1;

    bool connectData = false;
    char buffer[LINE_LENGTH];
    string oneLine;
    int index, test;
    char *context;
    // read whole file
    while (!in_file.eof()) {
        // char based line
        in_file.getline(buffer, LINE_LENGTH);
        // string based line
        oneLine = string(buffer);
        if (oneLine.length() > 7) {

            if (oneLine.substr(0, 5) == string("ATOM ") || oneLine.substr(0, 7) == string("HETATM ")) {
                atom toAdd;
                // invoke method to read the line
                toAdd.ReadPDBLine(oneLine);
                toAdd.ID = ID;
                // read index from file
                index = atoi(oneLine.substr(6, 6).c_str());
                atoms.insert(pair<int, atom > (index, toAdd));
            }
            else if (oneLine.substr(0, 7) == string("CONECT ")) {
                // there is connection data in this PDB file
                connectData = true;
                memset(buffer, ' ', 6);
                index = strtol(buffer, &context, 10);
                test = strtol(context, &context, 10);
                while (test != 0) {
                    atoms[index].IndexArray.insert(test);
                    test = strtol(context, &context, 10);
                }
            }
            else if (oneLine.substr(0, 6) == string("HEADER")) {
                comments.push_back(oneLine.substr(10));
            }
            else if (oneLine.substr(0, 7) == string("REMARK ") && oneLine.substr(10, 3) == "FIX") {
                comments.push_back(oneLine.substr(10));
            }

            /*else if (oneLine.substr(0,7) == string("REMARK ")) {
                    // overwrite the ID generator if found
                    if (oneLine.substr(11, 3) == string("ID ")) {
                            memset(buffer, ' ', 14);
                            index = strtol(buffer, &context, 10);
                            ID = strtol(context, &context, 10);
                            atoms[index].ID = ID;
                    }
            #ifdef OBJECT_SPACE
                    if (oneLine.substr(11,5) == string("ROTX "))
                            mole_axis.n[0] = atof(oneLine.substr(18).c_str());
                    if (oneLine.substr(11,5) == string("ROTY "))
                            mole_axis.n[1] = atof(oneLine.substr(18).c_str());
                    if (oneLine.substr(11,5) == string("ROTZ "))
                            mole_axis.n[2] = atof(oneLine.substr(18).c_str());
                    if (oneLine.substr(11,7) == string("ANGLE "))
                            mole_angle = atof(oneLine.substr(19).c_str());
            #endif
            }*/
        }
    }
    in_file.close();

    // when there is no connection information, generate bonds
    if (!connectData && !SkipConnectGenerate) {
        CreateBonds();
    }
#ifdef OBJECT_SPACE
    if (mole_axis == Vec3d()) {
        mole_axis.n[1] = 1;
        mole_angle = 0;
    }
#endif
    return 0;
}

void Ligand::AddHydrogen()
{
    string carbon("C"), hydrogen("H"), nitrogen("N"), oxygen("O");

    map<int, atom>::iterator it;
    set<int>::iterator iter;

    // need to guess where the hydrogen lies...
    bond_library library;
    int nextIndex, ring_count, bonds = 0;
    Vec3d v1, v2, v3, normal, delta, delta2;
    Mat4d rot;
    // check whether there is aromatic structure in this molecule
    list<int> ring;
    set<int> elements_in_ring;
    ring_count = DetectAromatic(ring);
    if (ring_count) {
        for (int i = 0; i != ring_count; ++i) {
            DetectAromatic(ring, i);
            for (list<int>::iterator list_it = ring.begin(); list_it != ring.end(); ++list_it)
                elements_in_ring.insert(*list_it);
        }
    }
    // consider backbone atoms only (C,N,O)
    // loop through each atom to see if a hydrogen could be added
    for (it = atoms.begin(); it != atoms.end(); ++it) {
        bonds = 0;
        // template of a hydrogen
        atom toAdd;
        toAdd.element = hydrogen;
        toAdd.name = hydrogen + "  ";
        toAdd.Residue = it->second.Residue;
        nextIndex = MaxIndex() + 1;
        toAdd.PDBIndex = nextIndex;
        toAdd.ID = atoms.begin()->second.ID;
        toAdd.IndexArray.clear();
        // assume connection to the testing atom
        toAdd.IndexArray.insert(it->first);
        // get number of bonds
        for (iter = it->second.IndexArray.begin(); iter != it->second.IndexArray.end(); ++iter)
            bonds += library.type(it->second.element, string(atoms[*iter].element), it->second.DistanceTo(atoms[*iter]));
        // it is part of aromatic ring, only one way to add hydrogen atoms
        if (elements_in_ring.find(it->first) != elements_in_ring.end()) {
            // already saturated in this configuration
            if (it->second.IndexArray.size() > 2) continue;
            iter = it->second.IndexArray.begin();
            v1 = atoms[*iter].coordinates - it->second.coordinates;
            ++iter;
            v2 = atoms[*iter].coordinates - it->second.coordinates;
            normal = v2^v1;
            double angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
            angle = M_PI - angle / 2;
            // maximise atoms apart
            rot = rot.createRotation(angle, normal);
            v3 = rot*v1;
            v3.normalize();
            delta = library.length(it->second.element, hydrogen) * v3;
            toAdd.coordinates = it->second.coordinates + delta;
            // try the opposite
            rot = rot.createRotation(-angle, normal);
            v3 = rot*v1;
            v3.normalize();
            delta2 = library.length(it->second.element, hydrogen) * v3;
            v2.normalize();
            // check orientation, maximise distance
            if ((delta2 - v2).length2() > (delta - v2).length2())
                toAdd.coordinates = it->second.coordinates + delta2;
            atoms.insert(pair<int, atom > (nextIndex, toAdd));
            it->second.IndexArray.insert(nextIndex);
            continue;
        }
        // testing on carbon atom
        if (it->second.element == carbon) {
            // already saturated
            if (it->second.IndexArray.size() == 4) continue;
            // find out number of shared electron

            // maximum electron sharing reached
            if (bonds == 4) continue;
            if (bonds == 3) {
                // triple bond, add on the opposite
                if (it->second.IndexArray.size() == 1) {
                    // find displacement to move
                    v1 = atoms[*(it->second.IndexArray.begin())].coordinates - it->second.coordinates;
                    v1.normalize();
                    delta = library.length(carbon, hydrogen) * -v1;
                    // get it to place
                    toAdd.coordinates = it->second.coordinates + delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    // add back connection
                    it->second.IndexArray.insert(nextIndex);
                }
                // one connected with single bond and the other with double bound
                if (it->second.IndexArray.size() == 2) {
                    iter = it->second.IndexArray.begin();
                    v1 = atoms[*iter].coordinates - it->second.coordinates;
                    ++iter;
                    v2 = atoms[*iter].coordinates - it->second.coordinates;
                    normal = v2^v1;
                    double angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
                    angle = M_PI - angle / 2;
                    // atoms are 120 degrees apart on a plane
                    rot = rot.createRotation(angle, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta = library.length(carbon, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // try the opposite
                    rot = rot.createRotation(-angle, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta2 = library.length(it->second.element, hydrogen) * v3;
                    v2.normalize();
                    // check orientation, maximise distance
                    if ((delta2 - v2).length2() > (delta - v2).length2())
                        toAdd.coordinates = it->second.coordinates + delta2;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
                // all are connected through single bonds
                if (it->second.IndexArray.size() == 3) {
                    iter = it->second.IndexArray.begin();
                    v1 = atoms[*iter].coordinates - it->second.coordinates;
                    v1.normalize();
                    ++iter;
                    v2 = atoms[*iter].coordinates - it->second.coordinates;
                    v2.normalize();
                    ++iter;
                    v3 = atoms[*iter].coordinates - it->second.coordinates;
                    v3.normalize();
                    // average the 3 vectors to produce optimal direction
                    normal = (v1^v2) + (v2^v3) + (v3^v1);
                    normal.normalize();
                    delta = library.length(carbon, hydrogen) * normal;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // check orientation
                    if ((it->second.coordinates - delta - atoms[*iter].coordinates).length2() > (toAdd.coordinates - atoms[*iter].coordinates).length2())
                        toAdd.coordinates = it->second.coordinates - delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
            }
            if (bonds == 2) {
                // connected through a double bond
                if (it->second.IndexArray.size() == 1) {
                    // extract connected atom and find a reference atom
                    atom nearAtom = atoms[*(it->second.IndexArray.begin())];
                    int refPoint = -1;
                    for (iter = nearAtom.IndexArray.begin(); iter != nearAtom.IndexArray.end(); ++iter) {
                        if ((*iter) != it->first) {
                            refPoint = *iter;
                            break;
                        }
                    }
                    // calculate the normal to rotate
                    v1 = atoms[refPoint].coordinates - nearAtom.coordinates;
                    v2 = it->second.coordinates - nearAtom.coordinates;
                    normal = v2^v1;
                    // add a hydrogen to both +120 degree and -120 degree
                    // hydrogens lie parallel to the reference vector
                    rot = rot.createRotation(2 * M_PI / 3, normal);
                    v3 = rot*v2;
                    v3.normalize();
                    delta = library.length(carbon, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                    toAdd.PDBIndex = ++nextIndex;
                    rot = rot.createRotation(-2 * M_PI / 3, normal);
                    v3 = rot*v2;
                    v3.normalize();
                    delta = library.length(carbon, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
                // connected through two single bonds
                if (it->second.IndexArray.size() == 2) {
                    // find normal of two bonds
                    iter = it->second.IndexArray.begin();
                    v1 = atoms[*iter].coordinates - it->second.coordinates;
                    ++iter;
                    v2 = atoms[*iter].coordinates - it->second.coordinates;
                    normal = v2^v1;
                    double angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
                    angle = M_PI - angle / 2;
                    // get vector by rotating 120 degrees around normal, outer side
                    rot = rot.createRotation(angle, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta = v3;
                    // try the opposite
                    rot = rot.createRotation(-angle, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta2 = v3;
                    // make fair comparison
                    v2.normalize();
                    // check orientation, maximise distance
                    if ((delta2 - v2).length2() > (delta - v2).length2()) {
                        v3 = delta2;
                    }
                    else {
                        v3 = delta;
                    }
                    // find vector parallel to the vectors plane
                    normal = normal^v3;
                    // rotate vector up and down for 54.75 (half of 109.5) degrees to get placement
                    rot = rot.createRotation(54.75 / 180 * M_PI, normal);
                    v1 = rot*v3;
                    v1.normalize();
                    delta = library.length(carbon, hydrogen) * v1;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // check orientation, maximise distance
                    if ((it->second.coordinates - delta - atoms[*iter].coordinates).length2() > (toAdd.coordinates - atoms[*iter].coordinates).length2())
                        toAdd.coordinates = it->second.coordinates - delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                    toAdd.PDBIndex = ++nextIndex;
                    rot = rot.createRotation(-54.75 / 180 * M_PI, normal);
                    v2 = rot*v3;
                    v2.normalize();
                    delta = library.length(carbon, hydrogen) * v2;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // check orientation, maximise distance
                    if ((it->second.coordinates - delta - atoms[*iter].coordinates).length2() > (toAdd.coordinates - atoms[*iter].coordinates).length2())
                        toAdd.coordinates = it->second.coordinates - delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
            }
            if (bonds == 1) {
                // need to add three hydrogens
                // extract connected atom and find a reference atom
                atom nearAtom = atoms[*(it->second.IndexArray.begin())];
                int refPoint = -1;
                for (iter = nearAtom.IndexArray.begin(); iter != nearAtom.IndexArray.end(); ++iter) {
                    if ((*iter) != it->first) {
                        refPoint = *iter;
                        break;
                    }
                }
                if (refPoint == -1) continue;
                // calculate the normal to rotate
                v1 = atoms[refPoint].coordinates - nearAtom.coordinates;
                v2 = it->second.coordinates - nearAtom.coordinates;
                normal = -v2^v1;
                // produce 2 siblings, 120 degrees apart
                rot = rot.createRotation(2 * M_PI / 3, v2);
                v3 = rot*normal;
                rot = rot.createRotation(-2 * M_PI / 3, v2);
                delta = rot*normal;
                // rotate to correct orientation (19.5 degrees), produce 3
                rot = rot.createRotation(19.5 / 180 * M_PI, normal^v2);
                normal = rot*normal;
                normal.normalize();
                v1 = library.length(carbon, hydrogen) * normal;
                toAdd.coordinates = it->second.coordinates + v1;
                // check orientation
                if ((it->second.coordinates - v1 - nearAtom.coordinates).length2() > (toAdd.coordinates - nearAtom.coordinates).length2())
                    toAdd.coordinates = it->second.coordinates - v1;
                atoms.insert(pair<int, atom > (nextIndex, toAdd));
                it->second.IndexArray.insert(nextIndex);
                // second hydrogen
                toAdd.PDBIndex = ++nextIndex;
                rot = rot.createRotation(19.5 / 180 * M_PI, v3^v2);
                v3 = rot*v3;
                v3.normalize();
                v1 = library.length(carbon, hydrogen) * v3;
                toAdd.coordinates = it->second.coordinates + v1;
                // check orientation
                if ((it->second.coordinates - v1 - nearAtom.coordinates).length2() > (toAdd.coordinates - nearAtom.coordinates).length2())
                    toAdd.coordinates = it->second.coordinates - v1;
                atoms.insert(pair<int, atom > (nextIndex, toAdd));
                it->second.IndexArray.insert(nextIndex);
                // third hydrogen
                toAdd.PDBIndex = ++nextIndex;
                rot = rot.createRotation(19.5 / 180 * M_PI, delta^v2);
                delta = rot*delta;
                delta.normalize();
                v1 = library.length(carbon, hydrogen) * delta;
                toAdd.coordinates = it->second.coordinates + v1;
                // check orientation
                if ((it->second.coordinates - v1 - nearAtom.coordinates).length2() > (toAdd.coordinates - nearAtom.coordinates).length2())
                    toAdd.coordinates = it->second.coordinates - v1;
                atoms.insert(pair<int, atom > (nextIndex, toAdd));
                it->second.IndexArray.insert(nextIndex);
            }
        }
        if (it->second.element == nitrogen) {
            // saturated in a stable compound
            if (it->second.IndexArray.size() >= 3) continue;
            // find out number of shared electron
            if (it->second.IndexArray.size() == 2) {
                // SP2 orbital
                if (bonds == 3) {
                    iter = it->second.IndexArray.begin();
                    v1 = atoms[*iter].coordinates - it->second.coordinates;
                    ++iter;
                    v2 = atoms[*iter].coordinates - it->second.coordinates;
                    normal = v2^v1;
                    double angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
                    angle = M_PI - angle / 2;
                    // rotate 120 degrees
                    rot = rot.createRotation(angle, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta = library.length(nitrogen, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // try the opposite
                    rot = rot.createRotation(-angle, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta2 = library.length(it->second.element, hydrogen) * v3;
                    v2.normalize();
                    // check orientation, maximise distance
                    if ((delta2 - v2).length2() > (delta - v2).length2())
                        toAdd.coordinates = it->second.coordinates + delta2;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
                // trigonal pyramid, 107 degrees among bonds
                if (bonds == 2) {
                    iter = it->second.IndexArray.begin();
                    v1 = atoms[*iter].coordinates - it->second.coordinates;
                    ++iter;
                    v2 = atoms[*iter].coordinates - it->second.coordinates;
                    normal = v2^v1;
                    // project v1 along the expected vector/v1 plane
                    v3 = v1 * cos(107 / 180 * M_PI);
                    // project v1 along v2 onto normal/expected vector plane
                    delta = v1 * cos(53.5 / 180 * M_PI);
                    double angle = acos(v3.length2() / delta.length2());
                    rot = rot.createRotation(M_PI - angle, (normal^(v1 + v2)));
                    v3 = rot * (v1 + v2);
                    v3.normalize();
                    delta = library.length(nitrogen, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // check orientation, maximise distance
                    if ((it->second.coordinates - delta - atoms[*iter].coordinates).length2() > (toAdd.coordinates - atoms[*iter].coordinates).length2())
                        toAdd.coordinates = it->second.coordinates - delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
            }
            if (it->second.IndexArray.size() == 1) {
                // triple bonded nitrogen does not link to hydrogen
                // SP2 orbital
                if (bonds == 2) {
                    iter = it->second.IndexArray.begin();
                    v1 = atoms[*iter].coordinates - it->second.coordinates;
                    ++iter;
                    v2 = atoms[*iter].coordinates - it->second.coordinates;
                    normal = v2^v1;
                    // rotate 120 degrees
                    rot = rot.createRotation(2 * M_PI / 3, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta = library.length(nitrogen, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                    // rotate -120 degrees
                    toAdd.PDBIndex = ++nextIndex;
                    rot = rot.createRotation(-2 * M_PI / 3, normal);
                    v3 = rot*v1;
                    v3.normalize();
                    delta = library.length(nitrogen, hydrogen) * v3;
                    toAdd.coordinates = it->second.coordinates + delta;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
                // trigonal pyramid, 107 degrees among bonds
                if (bonds == 1) {
                    atom nearAtom = atoms[*(it->second.IndexArray.begin())];
                    // pick another point to make reference vector
                    int refPoint = -1;
                    for (iter = nearAtom.IndexArray.begin(); iter != nearAtom.IndexArray.end(); ++iter) {
                        if ((*iter) != it->first) {
                            refPoint = *iter;
                            break;
                        }
                    }
                    // calculate the normal to rotate
                    v1 = atoms[refPoint].coordinates - nearAtom.coordinates;
                    v2 = it->second.coordinates - nearAtom.coordinates;
                    normal = v2^v1;
                    // add a hydrogen to 107 degree
                    rot = rot.createRotation(107 / 180 * M_PI, normal);
                    v3 = rot*v2;
                    v3.normalize();
                    delta = library.length(nitrogen, hydrogen) * v3;
                    v1 = delta;
                    // make trigonal pyramid
                    v2 = -v2;
                    normal = v2^v1;
                    // project v1 along the expected vector/v1 plane
                    v3 = v1 * cos(107 / 180 * M_PI);
                    // project v1 along v2 onto normal/expected vector plane
                    delta = v1 * cos(53.5 / 180 * M_PI);
                    double angle = acos(v3.length2() / delta.length2());
                    rot = rot.createRotation(M_PI - angle, (normal^(v1 + v2)));
                    v3 = rot * (v1 + v2);
                    v3.normalize();
                    delta = library.length(nitrogen, hydrogen) * v3;
                    // when connected atom is SP2, maximize distance from the plane
                    if (nearAtom.isSP2()) {
                        // project onto the rotation plane
                        rot = rot.createRotation(17 / 180 * M_PI, (v1^v2));
                        v3 = rot*v1;
                        rot = rot.createRotation(17 / 180 * M_PI, (delta^v2));
                        normal = rot*delta;
                        // get angle from dot product
                        angle = acos((v3 * normal) / sqrt(v3.length2() * normal.length2()));
                        // rotate both vector to appropiate new position
                        rot = rot.createRotation((M_PI - angle) / 2, v2);
                        v1 = rot*v1;
                        delta = rot*delta;
                    }
                    else if (nearAtom.IndexArray.size() == 4) {
                        toAdd.coordinates = it->second.coordinates + v1;
                        angle = DihedralAngle(toAdd, it->second, nearAtom, atoms[refPoint]);
                        rot = rot.createRotation(M_PI / 3 - angle, v2);
                        v1 = rot*v1;
                        delta = rot*delta;
                    }
                    else if (nearAtom.IndexArray.size() == 3) {
                        toAdd.coordinates = it->second.coordinates + v1;
                        angle = DihedralAngle(toAdd, it->second, nearAtom, atoms[refPoint]);
                        rot = rot.createRotation(-angle, v2);
                        v1 = rot*v1;
                        delta = rot*delta;
                    }
                    // add atoms to the map
                    toAdd.coordinates = it->second.coordinates + v1;
                    // check orientation
                    if ((it->second.coordinates - v1 - nearAtom.coordinates).length2() > (toAdd.coordinates - nearAtom.coordinates).length2())
                        toAdd.coordinates = it->second.coordinates - v1;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                    toAdd.PDBIndex = ++nextIndex;
                    toAdd.coordinates = it->second.coordinates + delta;
                    // check orientation
                    if ((it->second.coordinates - v1 - nearAtom.coordinates).length2() > (toAdd.coordinates - nearAtom.coordinates).length2())
                        toAdd.coordinates = it->second.coordinates - v1;
                    atoms.insert(pair<int, atom > (nextIndex, toAdd));
                    it->second.IndexArray.insert(nextIndex);
                }
            }
        }
        if (it->second.element == oxygen) {
            // should have saturated already
            if (it->second.IndexArray.size() >= 2) continue;
            // carboxyal group, linear on opposite side
            v1 = atoms[*(it->second.IndexArray.begin())].coordinates - it->second.coordinates;
            v1.normalize();
            delta = library.length(oxygen, hydrogen) * -v1;
            toAdd.coordinates = it->second.coordinates + delta;
            atoms.insert(pair<int, atom > (nextIndex, toAdd));
            it->second.IndexArray.insert(nextIndex);
        }
    }
}

int Ligand::DetectAromatic(list<int>& ring, int index)
{
    string carbon("C");
    int ring_count = 0;
    // used to test whether some atoms lie on the some plane
    Vec3d normal;
    // copy feasible atoms from the molecule, work till none is left
    map<int, atom> possible_list;
    map<int, atom>::iterator it;
    // vector of identified ring
    static vector<list<int> > vector_of_rings;
    vector<list<int> >::iterator vector_it;
    // count number of occurence of each node
    multiset<int> set_counter;
    multiset<int>::iterator multiset_it;
    set<int>::iterator set_it;
    list<int>::iterator list_it, list_it2;
    atom* cur_atom;
    bool flag = false;
    // filter out element that cannot be in the ring
    for (it = atoms.begin(); it != atoms.end(); ++it)
        if (it->second.IndexArray.size() != 1)
            possible_list.insert(*it);
    // clear the vector when the expected ring is not indicated
    if (index == -1) {
        for (vector_it = vector_of_rings.begin(); vector_it != vector_of_rings.end(); ++vector_it)
            (*vector_it).clear();
        vector_of_rings.clear();
    }
    else {
        // retrieve ring indicated by index
        if (static_cast<unsigned int> (index) < vector_of_rings.size())
            ring = vector_of_rings[index];
        return vector_of_rings.size();
    }
    // scan over the possible list, put them into a multiset
    for (it = possible_list.begin(); it != possible_list.end(); ++it)
        for (set_it = it->second.IndexArray.begin(); set_it != it->second.IndexArray.end(); ++set_it)
            if (possible_list.find(*set_it) != possible_list.end())
                set_counter.insert(*set_it);
    // filter out outliners
    it = possible_list.begin();
    while (it != possible_list.end()) {
        // element that is referred once is not part of a cycle
        if (set_counter.count(it->first) == 1)
            possible_list.erase(it++);
        else
            ++it;
    }
    set_counter.clear();
    // retry the possible list
    for (it = possible_list.begin(); it != possible_list.end(); ++it)
        for (set_it = it->second.IndexArray.begin(); set_it != it->second.IndexArray.end(); ++set_it)
            if (possible_list.find(*set_it) != possible_list.end())
                set_counter.insert(*set_it);

    pair<int, int> ref_index;
    list<int> possible_ring;
    list<int> normal_type;
    set<int> ring_element;
    vector<Vec3d> vector_of_normal;
    double result;
    int best_index;
    // work on possible indice
    while (!possible_list.empty()) {
        if (possible_list.size() <= 2) {
            possible_list.clear();
            break;
        }
        // find an indexed atom with 2 connections only
        it = possible_list.begin();
        while (it != possible_list.end() && set_counter.count(it->first) != 2)
            ++it;
        if (it == possible_list.end()) {
            possible_list.clear();
            break;
        }
        cur_atom = &(it->second);
        possible_ring.clear();
        possible_ring.push_back(it->first);
        // traverse around the ring, obtain a list
        traverse_ring(possible_ring, possible_list, it->first);
        // cleanse the ring with unwanted endpoints
        ring_element.clear();
        for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it) {
            ring_count = 0;
            cur_atom = &atoms[*list_it];
            for (list_it2 = possible_ring.begin(); list_it2 != possible_ring.end(); ++list_it2) {
                if (cur_atom->IndexArray.find(*list_it2) != cur_atom->IndexArray.end())
                    ++ring_count;
            }
            // a traversed element has only 1 connection, remove it
            if (ring_count == 1) {
                ring_element.insert(*list_it);
                possible_list.erase(*list_it);
            }
        }
        list_it = possible_ring.begin();
        while (list_it != possible_ring.end()) {
            if (ring_element.find(*list_it) != ring_element.end())
                possible_ring.erase(list_it++);
            else
                ++list_it;
        }
        ring_element.clear();
        // debug
        /*cout << "The traversed ring has elements: ";
        for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it)
                cout << *list_it << " ";
        cout << endl;*/
        // debug end
        // initialise parameters
        vector_of_normal.clear();
        normal_type.clear();

        // if the ring contains only 2 element, remove it and restart testing
        if (possible_ring.size() <= 2) {
            while (!possible_ring.empty()) {
                possible_list.erase(*(possible_ring.begin()));
                possible_ring.pop_front();
            }
            continue;
        }

        // make a set of ring elements for fast comparison
        ring_element.clear();
        for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it)
            ring_element.insert(*list_it);

        // check their orientation
        list_it = possible_ring.begin();
        while (list_it != possible_ring.end()) {
            cur_atom = &atoms[*list_it];
            // need to check the connection in the list
            ref_index = pair<int, int>(-1, -1);
            for (set_it = cur_atom->IndexArray.begin(); set_it != cur_atom->IndexArray.end(); ++set_it) {
                flag = (ring_element.find(*set_it) != ring_element.end());
                if (flag && ref_index.first == -1)
                    ref_index.first = *set_it;
                else if (flag && ref_index.second == -1)
                    ref_index.second = *set_it;
            }
            // get 2 adjacent indice that are in the list
            if (ref_index.second == -1) {
                possible_ring.erase(list_it++);
                continue;
            }
            // the normal of this atom
            normal = (atoms[ref_index.first].coordinates - cur_atom->coordinates)^(atoms[ref_index.second].coordinates - cur_atom->coordinates);
            normal.normalize();
            best_index = -1;
            if (!vector_of_normal.empty()) {
                // compare with available normals
                for (int i = 0; i != vector_of_normal.size(); ++i) {
                    result = normal * vector_of_normal[i];
                    if (result > 0.95 || result < -0.95)
                        best_index = i;
                }
            }
            // a close match found, otherwise push a new vector
            if (best_index != -1) {
                normal_type.push_back(best_index);
            }
            else {
                normal_type.push_back(vector_of_normal.size());
                vector_of_normal.push_back(normal);
            }
            ++list_it;
        }

        ring_element.clear();
        // only accept ring with 7 or less members
        for (int i = 0; i != vector_of_normal.size(); ++i) {
            ring_count = (int) count(normal_type.begin(), normal_type.end(), i);
            if (ring_count <= 2) {
                list_it = normal_type.begin();
                list_it2 = possible_ring.begin();
                // remove elements that has 2 or less in same orientation
                while (list_it != normal_type.end()) {
                    if (*list_it == i) {
                        normal_type.erase(list_it++);
                        // remove element from outermost list
                        possible_list.erase(*list_it2);
                        possible_ring.erase(list_it2++);
                    }
                    else {
                        ++list_it;
                        ++list_it2;
                    }
                }
            }
            else {
                // consider this normal has enough supporting nodes
                ring_element.insert(i);
            }
        }
        for (set_it = ring_element.begin(); set_it != ring_element.end(); ++set_it) {
            // count number of elements of the valid normal
            ring_count = (int) count(normal_type.begin(), normal_type.end(), *set_it);
            list<int> final_ring;
            list_it = normal_type.begin();
            list_it2 = possible_ring.begin();
            // remove elements that has moved to final ring
            while (list_it != normal_type.end()) {
                if (*list_it == *set_it) {
                    normal_type.erase(list_it++);
                    final_ring.push_back(*list_it2);
                    possible_ring.erase(list_it2++);
                }
                else {
                    ++list_it;
                    ++list_it2;
                }
            }
            // it is probably a double ring
            if (ring_count > 7) {
                list<int> local_ring;
                ring_count = split_ring(final_ring, local_ring);
                if (ring_count) {
                    for (int i = 0; i != ring_count; ++i) {
                        split_ring(final_ring, local_ring, i);
                        flag = checkHeteroAtom(local_ring);
                        if (flag) {
                            vector_of_rings.push_back(local_ring);
                        }
                        else {
                            for (list_it = local_ring.begin(); list_it != local_ring.end(); ++list_it)
                                possible_list.erase(*list_it);
                            local_ring.clear();
                        }
                    }
                }
                for (list_it = final_ring.begin(); list_it != final_ring.end(); ++list_it)
                    possible_list.erase(*list_it);
                if (ring_count) {
                    for (int i = 0; i != vector_of_rings.size(); ++i) {
                        for (list_it = vector_of_rings[i].begin(); list_it != vector_of_rings[i].end(); ++list_it) {
                            possible_list.erase(*list_it);
                        }
                    }
                }
            }
            else {
                flag = checkHeteroAtom(final_ring);
                for (list_it = final_ring.begin(); list_it != final_ring.end(); ++list_it) {
                    possible_list.erase(*list_it);
                }
                if (flag)
                    vector_of_rings.push_back(final_ring);
                else
                    final_ring.clear();
            }
        }
    }
    if (!vector_of_rings.empty())
        ring = vector_of_rings[0];
    return vector_of_rings.size();
}

int Ligand::DetectRing(list<int>& ring, int index)
{
    int ring_count = 0;
    // copy feasible atoms from the molecule, work till none is left
    map<int, atom> possible_list;
    map<int, atom>::iterator it;
    // vector of identified ring
    static vector<list<int> > vector_of_rings;
    vector<list<int> >::iterator vector_it;
    // count number of occurence of each node
    multiset<int> set_counter;
    multiset<int>::iterator multiset_it;
    set<int>::iterator set_it;
    list<int>::iterator list_it, list_it2;
    atom* cur_atom;
    bool flag = false;
    // filter out element that cannot be in the ring
    for (it = atoms.begin(); it != atoms.end(); ++it)
        if (it->second.IndexArray.size() != 1)
            possible_list.insert(*it);
    // clear the vector when the expected ring is not indicated
    if (index == -1) {
        for (vector_it = vector_of_rings.begin(); vector_it != vector_of_rings.end(); ++vector_it)
            (*vector_it).clear();
        vector_of_rings.clear();
    }
    else {
        // retrieve ring indicated by index
        if (static_cast<unsigned int> (index) < vector_of_rings.size())
            ring = vector_of_rings[index];
        return vector_of_rings.size();
    }
    // scan over the possible list, put them into a multiset
    for (it = possible_list.begin(); it != possible_list.end(); ++it)
        for (set_it = it->second.IndexArray.begin(); set_it != it->second.IndexArray.end(); ++set_it)
            if (possible_list.find(*set_it) != possible_list.end())
                set_counter.insert(*set_it);
    // filter out outliners
    it = possible_list.begin();
    while (it != possible_list.end()) {
        // element that is referred once is not part of a cycle
        if (set_counter.count(it->first) == 1)
            possible_list.erase(it++);
        else
            ++it;
    }
    set_counter.clear();
    // retry the possible list
    for (it = possible_list.begin(); it != possible_list.end(); ++it)
        for (set_it = it->second.IndexArray.begin(); set_it != it->second.IndexArray.end(); ++set_it)
            if (possible_list.find(*set_it) != possible_list.end())
                set_counter.insert(*set_it);

    list<int> possible_ring;
    set<int> ring_element;
    // work on possible indice
    while (!possible_list.empty()) {
        if (possible_list.size() <= 2) {
            possible_list.clear();
            break;
        }
        // find an indexed atom with 2 connections only
        it = possible_list.begin();
        while (it != possible_list.end() && set_counter.count(it->first) != 2)
            ++it;
        if (it == possible_list.end()) {
            possible_list.clear();
            break;
        }
        cur_atom = &(it->second);
        possible_ring.clear();
        possible_ring.push_back(it->first);
        // traverse around the ring, obtain a list
        traverse_ring(possible_ring, possible_list, it->first);
        // cleanse the ring with unwanted endpoints
        ring_element.clear();
        for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it) {
            ring_count = 0;
            cur_atom = &atoms[*list_it];
            for (list_it2 = possible_ring.begin(); list_it2 != possible_ring.end(); ++list_it2) {
                if (cur_atom->IndexArray.find(*list_it2) != cur_atom->IndexArray.end())
                    ++ring_count;
            }
            // a traversed element has only 1 connection, remove it
            if (ring_count == 1) {
                ring_element.insert(*list_it);
                possible_list.erase(*list_it);
            }
        }
        list_it = possible_ring.begin();
        while (list_it != possible_ring.end()) {
            if (ring_element.find(*list_it) != ring_element.end())
                possible_ring.erase(list_it++);
            else
                ++list_it;
        }
        ring_element.clear();
        // initialise parameters
        // if the ring contains only 2 element, remove it and restart testing
        if (possible_ring.size() <= 2) {
            while (!possible_ring.empty()) {
                possible_list.erase(*(possible_ring.begin()));
                possible_ring.pop_front();
            }
            continue;
        }

        // it is probably a double ring
        if (possible_ring.size() > 7) {
            list<int> local_ring;
            ring_count = split_ring(possible_ring, local_ring);
            if (ring_count) {
                for (int i = 0; i != ring_count; ++i) {
                    split_ring(possible_ring, local_ring, i);
                    vector_of_rings.push_back(local_ring);
                }
            }
            for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it)
                possible_list.erase(*list_it);
            if (ring_count) {
                for (int i = 0; i != vector_of_rings.size(); ++i) {
                    for (list_it = vector_of_rings[i].begin(); list_it != vector_of_rings[i].end(); ++list_it) {
                        possible_list.erase(*list_it);
                    }
                }
            }
        }
        else {
            // ensure the whole thing is connected
            flag = true;
            for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it) {
                ring_count = 0;
                cur_atom = &atoms[*list_it];
                // a heterocyclic ring will not contains too many heteroatom
                // a carbon atom should be in its unstaturated state to be aromatic
                cur_atom = &atoms[*list_it];
                for (list_it2 = possible_ring.begin(); list_it2 != possible_ring.end(); ++list_it2)
                    if (cur_atom->IndexArray.find(*list_it2) != cur_atom->IndexArray.end())
                        ++ring_count;
                // an atom is connected to exactly 2 others in this ring
                if (ring_count != 2)
                    flag = false;
            }
            for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it) {
                possible_list.erase(*list_it);
            }
            if (flag)
                vector_of_rings.push_back(possible_ring);
            else {
                // call split ring
                list<int> local_ring;
                ring_count = split_ring(possible_ring, local_ring);
                if (ring_count) {
                    for (int i = 0; i != ring_count; ++i) {
                        split_ring(possible_ring, local_ring, i);
                        vector_of_rings.push_back(local_ring);
                    }
                }
                for (list_it = possible_ring.begin(); list_it != possible_ring.end(); ++list_it)
                    possible_list.erase(*list_it);
                if (ring_count) {
                    for (int i = 0; i != vector_of_rings.size(); ++i) {
                        for (list_it = vector_of_rings[i].begin(); list_it != vector_of_rings[i].end(); ++list_it) {
                            possible_list.erase(*list_it);
                        }
                    }
                }
            }
        }
    }
    if (!vector_of_rings.empty())
        ring = vector_of_rings[0];
    return vector_of_rings.size();
}

int Ligand::synthesis(string FilenameOfFragment)
{
    // local molecule copy
    Ligand fragment;
    fragment.LoadPDB(FilenameOfFragment);

    // test if ring joining is successful, 0 is returned if successfully joined
    if (JoinRing(fragment) == 0) return 2;

    return -1;

    bond_library library;
    set<int> SP2atoms;
    map<int, int> double_bond_pairs;
    set<int>::iterator set_it;
    atom* cur_atom;
    int count(0), ringCount;
    list<int> ring;
    bool flag = false;

    // find out the set of atoms that are SP2
    for (map<int, atom>::iterator it = fragment.atoms.begin(); it != fragment.atoms.end(); ++it) {
        if (it->second.isSP2())
            SP2atoms.insert(it->first);
    }
    // detect number of aromatic ring in the fragment
    ringCount = DetectAromatic(ring);
    // remove all those SP2 atoms that are found in an aromatic ring structure
    if (ringCount) {
        for (int i = 0; i != ringCount; ++i) {
            DetectAromatic(ring, i);
            for (list<int>::iterator it = ring.begin(); it != ring.end(); ++it)
                SP2atoms.erase(*it);
        }
    }

    // there are enough atoms to have double bond
    if (SP2atoms.size() >= 2) {
        set_it = SP2atoms.begin();
        while (set_it != SP2atoms.end()) {
            cur_atom = &fragment.atoms[*set_it];
            flag = false;
            for (set<int>::iterator it = cur_atom->IndexArray.begin(); it != cur_atom->IndexArray.end(); ++it) {
                if (SP2atoms.find(*it) != SP2atoms.end()) {
                    double_bond_pairs.insert(pair<int, int>(*set_it, *it));
                    flag = true;
                }
            }
            // erase checked index
            if (flag)
                SP2atoms.erase(set_it++);
            else
                ++set_it;
        }
    }
    else {
        // indicate the failure of using this mode, use normal joining method
        return -1;
    }
    SP2atoms.clear();

    if (double_bond_pairs.empty())
        return -1;
    // give some chance that even this structure is eligible to synthesis we'll still reject it
//    if (generator() > 0.5)
//        return -1;

    // extract a pair to work on
    pair<int, int> chosen_pair;
    if (double_bond_pairs.size() > 1) {
        map<int, int>::iterator map_it = double_bond_pairs.begin();
//        ringCount = generator_int() % double_bond_pairs.size();
//        ringCount = generator_int() % double_bond_pairs.size();	
	ringCount = 0 % double_bond_pairs.size();	
        advance(map_it, ringCount);
        chosen_pair = *map_it;
    }
    else {
        chosen_pair = *(double_bond_pairs.begin());
    }

    // scan number of hydrogen on the pair, choose the one with more hydrogen
    count = 0;
    ringCount = 0;
    cur_atom = &fragment.atoms[chosen_pair.first];
    for (set_it = cur_atom->IndexArray.begin(); set_it != cur_atom->IndexArray.end(); ++set_it)
        if (fragment.atoms[*set_it].element == "H")
            ++count;
    cur_atom = &fragment.atoms[chosen_pair.second];
    for (set_it = cur_atom->IndexArray.begin(); set_it != cur_atom->IndexArray.end(); ++set_it)
        if (fragment.atoms[*set_it].element == "H")
            ++ringCount;

    // ensure the first element has more hydrogens
    if (ringCount > count)
        chosen_pair = pair<int, int>(chosen_pair.second, chosen_pair.first);
    // split the fragment into 2 smaller parts
    if (fragment.replace_bond(chosen_pair) == -1)
        return -1;

    // select hydrogen
    int connectIndex, fragHydrogen, fragIndex, linkerHydrogen(-1);
    // atom copy and atom reference
    atom curHydrogen, connectAtom, *fragHydrogenAtom, *fragConnectAtom;
    string element1, element2;
    Vec3d delta;

    // pick hydrogen of ligand
    linkerHydrogen = IndexOfRandomHydrogen();
//
    // obtain information on hydrogens
    curHydrogen = atoms[linkerHydrogen];
    // hydrogen has only 1 connection, get the first index in array
    connectIndex = *(curHydrogen.IndexArray.begin());
    connectAtom = atoms[connectIndex];
    // find hydrogen on the reacted bond
    set<int> pending_hydrogen;
    cur_atom = &fragment.atoms[chosen_pair.first];
    for (set_it = cur_atom->IndexArray.begin(); set_it != cur_atom->IndexArray.end(); ++set_it)
        if (fragment.atoms[*set_it].element == string("H"))
            pending_hydrogen.insert(*set_it);
    set_it = pending_hydrogen.begin();
//    advance(set_it, generator_int() % pending_hydrogen.size());
    advance(set_it, 0 % pending_hydrogen.size());        
    fragHydrogen = *set_it;
    // these atoms would be updated, get reference instead of copying
    fragHydrogenAtom = &fragment.atoms[fragHydrogen];
    fragIndex = *(fragHydrogenAtom->IndexArray.begin());
    fragConnectAtom = &fragment.atoms[fragIndex];

    // move fragment in place
    // Moving the fragment near the atom to which it will bind.
    // pick the bond to be adjusted
    element1 = connectAtom.element;
    element2 = fragConnectAtom->element;
    delta = curHydrogen.coordinates - connectAtom.coordinates;
    delta.normalize();
    // scale it to reflect real molecular bond length
    delta *= library.length(element1, element2);
    // move fragment in place
    fragment.Translate(fragIndex, (connectAtom.coordinates + delta));

    // rotate fragment to minimise distance between fragment's hydrogen and connected atom
    Vec3d A = connectAtom.coordinates - fragConnectAtom->coordinates;
    Vec3d B = fragHydrogenAtom->coordinates - fragConnectAtom->coordinates;
    // normal to the 3 points
    Vec3d normal = B^A;
    // angle between 3 points with respect to fragment connected atom
    double angle = acos((A * B) / sqrt(A.length2() * B.length2()));
    // rotate along the line using connected atom as pivot
    // prevent the situation when it is already in place in the beginning
    if (normal != Vec3d())
        fragment.RotateLine(normal, fragIndex, angle);

    // minimise hinderance for SP3 and SP bond
    // Rotate fragment around connecting bond to minimize steric hindrance...
    double BadContact, BestContact(0), BestAngle(0), Angle(0), step(M_PI / 25);
    normal = connectAtom.coordinates - fragConnectAtom->coordinates;
    // try a whole circle with definite steps
    while (Angle < 2 * M_PI) {
        Angle += step;
        if (normal == Vec3d()) break;
        fragment.RotateLine(normal, fragIndex, step);
        BadContact = MolecularDistance(fragment);
        // maximize the distance
        if (BadContact > BestContact) {
            BestContact = BadContact;
            BestAngle = Angle;
        }
    }
    // rotate to least hindered orientation
    if (normal != Vec3d())
        fragment.RotateLine(normal, fragIndex, BestAngle);


    // erging the fragment with the original molecule
    // remove hydrogen atom from both molecule
    DeleteAtom(linkerHydrogen);
    fragment.DeleteAtom(fragHydrogen);

    int updateIndex, cascadeIndex(MaxIndex());
    set<int> tempIndice;
    ostringstream output;

    for (map<int, atom>::iterator it = fragment.atoms.begin(); it != fragment.atoms.end(); ++it) {
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

    return 0;
}

bool Ligand::valid()
{    
    if (atoms.begin()->first == 0)
	atoms.erase(atoms.begin());
    if (atoms.size() > 80) return false; // 80 should be replaced by max_atoms
//    int donor, acceptor;
//    double weight, logP;
//    if (!Lipinski5(ref, donor, acceptor, weight, logP)) return false;
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
	if (it->second.coordinates[0] != it->second.coordinates[0]) return false; // Elminate NaN	    
    return true;
}

double Ligand::mwt()
{
    double MW(0);
    bond_library library;
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it)
    {
	MW += library.weight(it->second.element);
    }
    return MW;
}

// protected methods

bool Ligand::hasBadBonds()
{
    bool badContact(false), badLength(false);
    double dist;
    map<int, atom>::iterator iter, it;
    // account number of connection, it has maximum connected atom based on orbital
    for (iter = atoms.begin(); iter != atoms.end(); ++iter) {
        if (iter->second.element[0] == 'C' && iter->second.element.size() == 1 && iter->second.IndexArray.size() > 4)
            badContact = true;
        if (iter->second.element[0] == 'O' && iter->second.IndexArray.size() > 2)
            badContact = true;
        if (iter->second.element[0] == 'H' && iter->second.IndexArray.size() > 1)
            badContact = true;
        if (iter->second.element[0] == 'N' && iter->second.IndexArray.size() > 4)
            badContact = true;
    }
    // atoms cannot be too close which would cause steric clash
    for (iter = atoms.begin(); iter != atoms.end(); ++iter) {
        it = iter;
        ++it;
        while (it != atoms.end()) {
            dist = iter->second.DistanceTo(it->second);
            if (dist < 1.15) {
                if (iter->second.IndexArray.find(it->first) == iter->second.IndexArray.end())
                    badLength = true;
            }
            ++it;
        }
    }
    return (badContact || badLength);
}

// as it is a map, largest indexed element is at the back

int Ligand::MaxIndex()
{
    map<int, atom>::reverse_iterator rit;
    rit = atoms.rbegin();
    return rit->first;
}

// as it is a map, smallest indexed element is at the front

int Ligand::MinIndex()
{
    map<int, atom>::iterator it;
    it = atoms.begin();
    return it->first;
}

// will invalidate iterator, return updated one

map<int, atom>::reverse_iterator Ligand::DeleteAtom(int index)
{
    // find iterator to the atom using index
    map<int, atom>::iterator it = atoms.find(index);
    set<int>::iterator index_it;
    // return error when not found
    if (it == atoms.end()) return atoms.rend();
    // loop through connected indices
    for (index_it = it->second.IndexArray.begin(); index_it != it->second.IndexArray.end(); ++index_it) {
        // remove index from the indicated atoms
        atoms[*index_it].IndexArray.erase(index);
    }
    // removed indexed atom, can only erase using forward iterator
    atoms.erase(it++);
    map<int, atom>::reverse_iterator rit(it);
    return rit;
}

// total minimum distance for all atoms, quadratic

double Ligand::MolecularDistance(Ligand &other)
{
    double min, dist, total_dist = 0;
    map<int, atom>::iterator it1, it2;
    for (it1 = atoms.begin(); it1 != atoms.end(); ++it1) {
        min = DBL_MAX;
        for (it2 = other.atoms.begin(); it2 != other.atoms.end(); ++it2) {
            dist = it1->second.DistanceTo(it2->second);
            if (dist < min)
                min = dist;
        }
        total_dist += dist;
    }
    return total_dist;
}

// intramolecular distance, dRMSD quadratic

double Ligand::IntraMolecularDistance(Ligand &other)
{
    // needs to work on same number of atoms
    //if (atoms.size() != other.atoms.size()) return DBL_MAX;
    string carbon("C"), nitrogen("N"), oxygen("O"), hydrogen("H");
    int c1(0), o1(0), n1(0), other1(0), c2(0), o2(0), n2(0), other2(0), total;
    double pij, qij, dist = 0;
    map<int, atom>::iterator it1, it2, it3, it4;
    // check identity of molecules
    for (it1 = atoms.begin(); it1 != atoms.end(); ++it1) {
        if (it1->second.element == hydrogen)
            continue;
        if (it1->second.element == carbon)
            ++c1;
        else if (it1->second.element == nitrogen)
            ++n1;
        else if (it1->second.element == oxygen)
            ++o1;
        else
            ++other1;
    }
    for (it2 = other.atoms.begin(); it2 != other.atoms.end(); ++it2) {
        if (it2->second.element == hydrogen)
            continue;
        if (it2->second.element == carbon)
            ++c2;
        else if (it2->second.element == nitrogen)
            ++n2;
        else if (it2->second.element == oxygen)
            ++o2;
        else
            ++other2;
    }
    // must match non-hydrogen atoms
    if (c1 != c2 || n1 != n2 || o1 != o2 || other1 != other2) return DBL_MAX;
    total = c1 + n1 + o1 + other1;
    // as index i
    it1 = atoms.begin();
    it2 = other.atoms.begin();
    while (it1 != atoms.end()) {
        if (it1->second.element == hydrogen) {
            ++it1;
            continue;
        }
        if (it2->second.element == hydrogen) {
            ++it2;
            continue;
        }
        // as index j, ensure i < j condition
        it3 = it1;
        ++it3;
        it4 = it2;
        ++it4;
        while (it3 != atoms.end()) {
            if (it3->second.element == hydrogen) {
                ++it3;
                continue;
            }
            if (it4->second.element == hydrogen) {
                ++it4;
                continue;
            }
            // distance between atoms in both conformation
            pij = it3->second.DistanceTo(it1->second);
            qij = it4->second.DistanceTo(it2->second);
            dist += (pij - qij)*(pij - qij);
            ++it3;
            ++it4;
        }
        ++it1;
        ++it2;
    }
    // rectified by N*(N+1)
    dist *= 1 / (total * (total + 1));
    // take square root
    dist = sqrt(dist);
    return dist;
}

// in radian, angle between 2 planes (a1-a2-a3, a2-a3-a4) with respect to vector a2-a3

double Ligand::DihedralAngle(atom a1, atom a2, atom a3, atom a4)
{
    Vec3d v1 = a2.coordinates - a1.coordinates;
    Vec3d v2 = a3.coordinates - a2.coordinates;
    Vec3d v3 = a4.coordinates - a3.coordinates;
    return atan2((v1 * v2.length())*(v2^v3), (v1^v2)*(v2^v3));
}

// angle between v1-v2-v3

double Ligand::Angle(Vec3d v1, Vec3d v2, Vec3d v3)
{
    Vec3d b1 = v1 - v2;
    Vec3d b2 = v3 - v2;
    return acos((b1 * b2) / sqrt(b1.length2() * b2.length2()));
}

// each atom contributes a weight based on molecular weight

Vec3d Ligand::CentreOfGravity()
{
    double total_weight = 0;
    Vec3d centre = Vec3d();
    bond_library library;
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        centre += it->second.coordinates * library.weight(it->second.element);
        total_weight += library.weight(it->second.element);
    }
    centre = centre / (atoms.size() * total_weight);
    return centre;
}

// produce bond based on molecular distance when connection data is not available

void Ligand::CreateBonds()
{
    int i = 0;
    map<int, atom>::iterator it1, it2;
    double dist;
    bond_library library;
    for (it1 = atoms.begin(); it1 != atoms.end(); ++it1) {
        ++i;
        it2 = atoms.begin();
        advance(it2, i);
        while (it2 != atoms.end()) {
            dist = it1->second.DistanceTo(it2->second);
            if (dist < library.length(it1->second.element, it2->second.element)*1.2) {
                it1->second.IndexArray.insert(it2->first);
                it2->second.IndexArray.insert(it1->first);
            }
            ++it2;
        }
        // skip last
        if (i == atoms.size() - 1) break;
    }
}

// traverse molecule to determine some parts to retain

void Ligand::scan_recursive(int index)
{
    // since the ring is stored in a static vector, can be retrieved after computing once
    list<int> ring;
    bool sameRing = false;
    // obtain number of ring in this molecule
    int count, ringCount = DetectRing(ring, 0);

    scanned.insert(index);
    set<int>::iterator it;
    atom toAdd, connect_atom, cur_atom = atoms[index];
    bond_library library;
    Vec3d delta;
    string cur_element;
    ostringstream convertor;

    // default replacement atom as hydrogen
    toAdd.name = "H  ";
    toAdd.element = "H";
    toAdd.ID = cur_atom.ID;

    for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it) {
        // this atom has not been visited
        if (scanned.find(*it) == scanned.end()) {
            sameRing = false;
            connect_atom = atoms[*it];
            // the current atom is in the overlap set, the connected atom is not in the overlap set
            if (overlap.find(index) != overlap.end() && overlap.find(*it) == overlap.end()) {
                // go through all rings
                for (int i = 0; i != ringCount; ++i) {
                    DetectRing(ring, i);
                    count = 0;
                    for (list<int>::iterator list_it = ring.begin(); list_it != ring.end(); ++list_it)
                        if (*list_it == index || *list_it == *it)
                            ++count;
                    // both current atom and pending atom are in the same ring, cannot split them
                    if (count == 2)
                        sameRing = true;
                }
                if (sameRing) {
                    scan_recursive(*it);
                    continue;
                }
                // make a replacement atom
                cur_element = cur_atom.element;
                toAdd.IndexArray.insert(index);
                convertor.str(string());
                convertor << index;
                toAdd.PDBIndex = convertor.str();
                toAdd.Residue = cur_atom.Residue;
                delta = connect_atom.coordinates - cur_atom.coordinates;
                delta.normalize();
                delta *= library.length("H", cur_element);
                toAdd.coordinates = cur_atom.coordinates + delta;
                toAddAtoms.insert(pair<int, atom > (*it, toAdd));
                // 50% chance to include this connected atom
//                if (generator() > 0.5)
//                    scan_recursive(*it);
                /*} else if (overlap.find(index) == overlap.end() && cur_atom.Residue != connect_atom.Residue) {
                        // the current atom is not in the overlap set, the connected atom is from different residue
                        // make a replacement atom
                        cur_element = cur_atom.element;
                        toAdd.IndexArray.insert(index);
                        convertor.str(string());
                        convertor << index;
                        toAdd.PDBIndex = convertor.str();
                        toAdd.Residue = cur_atom.Residue;
                        delta = connect_atom.coordinates - cur_atom.coordinates;
                        delta.normalize();
                        delta *= library.length("H", cur_element);
                        toAdd.coordinates = cur_atom.coordinates + delta;
                        toAddAtoms.insert(pair<int,Atom>(*it,toAdd));
                        // 50% chance to include this connected atom
                        if (generator() > 0.5)
                                scan_recursive(*it);
                } else if (overlap.find(index) == overlap.end() && cur_atom.ID != connect_atom.ID) {
                        // the two atoms are originated from different fragment
                        // make a replacement atom
                        cur_element = cur_atom.element;
                        toAdd.IndexArray.insert(index);
                        convertor.str(string());
                        convertor << index;
                        toAdd.PDBIndex = convertor.str();
                        toAdd.Residue = cur_atom.Residue;
                        delta = connect_atom.coordinates - cur_atom.coordinates;
                        delta.normalize();
                        delta *= library.length("H", cur_element);
                        toAdd.coordinates = cur_atom.coordinates + delta;
                        toAddAtoms.insert(pair<int,Atom>(*it,toAdd));
                        // 50% chance to include this connected atom
                        if (generator() > 0.5)
                                scan_recursive(*it);*/
            }
            else {
                scan_recursive(*it);
            }
        }
    }
}

// produce rotation matrix

void Ligand::Rotate(Vec3d centre, Vec3d EulerAngles)
{
    // produce rotation matrix
    Vec3d location;
    Mat3d rotx, roty, rotz;
    double cos_angle, sin_angle;
    rotx = Mat3d();
    cos_angle = cos(EulerAngles.n[0]);
    sin_angle = sin(EulerAngles.n[0]);
    rotx.n[4] = cos_angle;
    rotx.n[5] = -sin_angle;
    rotx.n[7] = sin_angle;
    rotx.n[8] = cos_angle;
    roty = Mat3d();
    cos_angle = cos(EulerAngles.n[1]);
    sin_angle = sin(EulerAngles.n[1]);
    roty.n[0] = cos_angle;
    roty.n[2] = sin_angle;
    roty.n[6] = -sin_angle;
    roty.n[8] = cos_angle;
    rotz = Mat3d();
    cos_angle = cos(EulerAngles.n[2]);
    sin_angle = sin(EulerAngles.n[2]);
    rotz.n[0] = cos_angle;
    rotz.n[1] = -sin_angle;
    rotz.n[3] = sin_angle;
    rotz.n[4] = cos_angle;
    // rotate ligand
    for (map<int, atom>::iterator it = atoms.begin(); it != atoms.end(); ++it) {
        location = it->second.coordinates - centre;
        location = rotx * roty * rotz*location;
        it->second.coordinates = location + centre;
    }
}

void Ligand::traverse_ring(list<int>& possible_ring, map<int, atom>& possible_list, int index)
{
    atom* cur_atom = &atoms[index];
    bool flag;
    // for all its neighbour
    for (set<int>::iterator it = cur_atom->IndexArray.begin(); it != cur_atom->IndexArray.end(); ++it) {
        flag = false;
        // find an index that is a possible element in ring but not yet covered
        if (possible_list.find(*it) != possible_list.end()) {
            for (list<int>::iterator lit = possible_ring.begin(); lit != possible_ring.end(); ++lit)
                if ((*it) == (*lit))
                    flag = true;
            // no such element in the ring yet
            if (!flag) {
                possible_ring.push_back(*it);
                traverse_ring(possible_ring, possible_list, *it);
            }
        }
    }
}

int Ligand::split_ring(list<int>& input, list<int>& output, int index)
{
    // ring information are stored in this method
    static vector<list<int> > rings;
    list<int> local_ring;
    if (index != -1) {
        if (static_cast<unsigned int> (index) < rings.size())
            output = rings[index];
        return rings.size();
    }
    else {
        rings.clear();
    }
    // declare variable and initialisation
    bool flag, hasNext;
    int count(0);
    atom* cur_atom;
    list<int>::iterator it1, it2;
    // connect2 to hold a set of indice that have 2 connections only, similarly for connect3, connect4 hold the rest of the case
    // any index that has a single connection has been rejected before calling this method
    // the indice in connect3 and connect4 are retained in further splitting
    set<int> connect2, connect3, connect4, scanned;
    set<int>::iterator set_it;
    // for each atom, extract the number of connections which have been regarded as candidate
    for (it1 = input.begin(); it1 != input.end(); ++it1) {
        count = 0;
        cur_atom = &atoms[*it1];
        // find the number of connected atoms to be a possible element of a ring
        for (it2 = input.begin(); it2 != input.end(); ++it2) {
            if (cur_atom->IndexArray.find(*it2) != cur_atom->IndexArray.end())
                ++count;
        }
        if (count == 2)
            connect2.insert(*it1);
        else if (count == 3)
            connect3.insert(*it1);
        else
            connect4.insert(*it1);
    }
    // all element are in the same ring
    if (connect3.empty() && connect4.empty()) {
        rings.push_back(input);
        return 1;
    }
    //cout << "Connections detected 2: " << connect2.size() << " 3: " << connect3.size() << " >4: " << connect4.size() << endl;

    // connect all 2-connected edges
    while (!connect2.empty()) {
        // start from a 2-connected atom
        local_ring.push_back(*connect2.begin());
        // remove this element
        connect2.erase(connect2.begin());
        // maintain a scanned list
        scanned.clear();
        scanned.insert(local_ring.back());
        // indicate whether a 2-connected atom is being considered
        flag = true;
        // find trailing atom
        while (flag) {
            // get the tailing element
            cur_atom = &atoms[local_ring.back()];
            hasNext = false;
            // from all possible element, find those to which this atom connects
            for (it1 = input.begin(); it1 != input.end(); ++it1) {
                if (cur_atom->IndexArray.find(*it1) != cur_atom->IndexArray.end() && scanned.find(*it1) == scanned.end()) {
                    // update scanned
                    scanned.insert(*it1);
                    // push entry to tail
                    local_ring.push_back(*it1);
                    // check this new element is 2-connected
                    if (connect2.find(*it1) == connect2.end()) flag = false;
                    // remove used entry
                    connect2.erase(*it1);
                    // avoid tracing the same element
                    hasNext = true;
                    // only trace one branch at a time
                    break;
                }
            }
            if (!hasNext) flag = false;
        }
        flag = true;
        // find leading atom
        while (flag) {
            cur_atom = &atoms[local_ring.front()];
            hasNext = false;
            // find a neighbour that is 
            for (it1 = input.begin(); it1 != input.end(); ++it1) {
                if (cur_atom->IndexArray.find(*it1) != cur_atom->IndexArray.end() && scanned.find(*it1) == scanned.end()) {
                    scanned.insert(*it1);
                    local_ring.push_front(*it1);
                    if (connect2.find(*it1) == connect2.end()) flag = false;
                    connect2.erase(*it1);
                    hasNext = true;
                    break;
                }
            }
            if (!hasNext) flag = false;
        }
        rings.push_back(local_ring);
        local_ring.clear();
    }
    merge_list(connect3, rings, input, 3);
    merge_list(connect4, rings, input, 4);
    return rings.size();
}

void Ligand::merge_list(set<int>& connect, vector<list<int> >& rings, list<int>& input, int limit)
{
    unsigned int count, ref_index;
    set<int>::iterator set_it;
    list<int>::iterator it1;
    set<int> scanned;
    // scan repetition in all the lists, merge them based on 3-connected and 4-connected
    while (!connect.empty()) {
        scanned.clear();
        // find out the occurrence of the first atom in connect3 across all list of rings
        for (int i = 0; i != rings.size(); ++i) {
            for (it1 = rings[i].begin(); it1 != rings[i].end(); ++it1)
                if (*it1 == *connect.begin())
                    scanned.insert(i);
        }
        // need to merge to form a complete ring
        if (scanned.size() == limit) {
            count = input.size();
            ref_index = 0;
            // lookup the shortest list
            for (set_it = scanned.begin(); set_it != scanned.end(); ++set_it) {
                if (rings[*set_it].size() < count) {
                    count = rings[*set_it].size();
                    ref_index = *set_it;
                }
            }
            scanned.erase(ref_index);
            for (set_it = scanned.begin(); set_it != scanned.end(); ++set_it) {
                // 4 cases to join the lists
                if (rings[*set_it].front() == *connect.begin()) {
                    if (rings[ref_index].front() == *connect.begin()) {
                        for (it1 = rings[ref_index].begin(); it1 != rings[ref_index].end(); ++it1)
                            rings[*set_it].push_front(*it1);
                    }
                    else {
                        for (list<int>::reverse_iterator it = rings[ref_index].rbegin(); it != rings[ref_index].rend(); ++it)
                            rings[*set_it].push_front(*it);
                    }
                }
                else {
                    if (rings[ref_index].front() == *connect.begin()) {
                        for (it1 = rings[ref_index].begin(); it1 != rings[ref_index].end(); ++it1)
                            rings[*set_it].push_back(*it1);
                    }
                    else {
                        for (list<int>::reverse_iterator it = rings[ref_index].rbegin(); it != rings[ref_index].rend(); ++it)
                            rings[*set_it].push_back(*it);
                    }
                }
            }
            rings[ref_index].clear();
            vector<list<int> >::iterator it = rings.begin();
            advance(it, ref_index);
            rings.erase(it);
        }
        else {
            // the connection is not fully utilised
            // go through the rings
            for (int i = 0; i != rings.size(); ++i) {
                // skip the ring that has found this element
                if (scanned.find(i) != scanned.end()) continue;
                // the ring is already completed, skip it
                if (atoms[rings[i].front()].IndexArray.find(rings[i].back()) != atoms[rings[i].front()].IndexArray.end()) continue;
                // from all its connection, find which to add
                for (set_it = atoms[*connect.begin()].IndexArray.begin(); set_it != atoms[*connect.begin()].IndexArray.end(); ++set_it)
                    // a connection is found to be the front of a ring
                    if (rings[i].front() == *set_it) {
                        rings[i].push_front(*connect.begin());
                        scanned.insert(i);
                        break;
                    }
                    else if (rings[i].back() == *set_it) {
                        rings[i].push_back(*connect.begin());
                        scanned.insert(i);
                        break;
                    }
            }
        }
        connect.erase(connect.begin());
    }
}

// assume a double bond is reacted into a single bond

int Ligand::replace_bond(pair<int, int> bond)
{
    Ligand frag1, frag2, frag3, frag4;
    int connect1(-1), connect2(-1), connect3(-1), connect4(-1), nextIndex = MaxIndex() + 1;
    atom* cur_atom;
    Vec3d v1, v2, v3, normal, delta, bond_vector;
    Mat4d rot;
    bond_library library;
    string hydrogen("H");
    map<int, atom> side1, side2;
    double angle, dist, max_dist;
    ostringstream convertor;
    // retrieve all the connected atom of the affected bond
    cur_atom = &atoms[bond.first];
    for (set<int>::iterator it = cur_atom->IndexArray.begin(); it != cur_atom->IndexArray.end(); ++it) {
        if (connect1 == -1 && *it != bond.second)
            connect1 = *it;
        else if (connect2 == -1 && *it != bond.second)
            connect2 = *it;
    }
    cur_atom = &atoms[bond.second];
    for (set<int>::iterator it = cur_atom->IndexArray.begin(); it != cur_atom->IndexArray.end(); ++it) {
        if (connect3 == -1 && *it != bond.first)
            connect3 = *it;
        else if (connect4 == -1 && *it != bond.first)
            connect4 = *it;
    }
    split_fragment(frag1, connect1, bond.first);
    split_fragment(frag2, connect2, bond.first);
    split_fragment(frag3, connect3, bond.second);
    split_fragment(frag4, connect4, bond.second);

    // the 4 fragments are not mutually exclusive
    if (frag1.atoms.size() + frag2.atoms.size() + frag3.atoms.size() + frag4.atoms.size() > atoms.size())
        return -1;

    /*cout << "There 4 fragments are:" << endl;
    for (map<int,Atom>::iterator it = frag1.atoms.begin(); it != frag1.atoms.end(); ++it)
            cout << it->first << " ";
    cout << endl;
    for (map<int,Atom>::iterator it = frag2.atoms.begin(); it != frag2.atoms.end(); ++it)
            cout << it->first << " ";
    cout << endl;
    for (map<int,Atom>::iterator it = frag3.atoms.begin(); it != frag3.atoms.end(); ++it)
            cout << it->first << " ";
    cout << endl;
    for (map<int,Atom>::iterator it = frag4.atoms.begin(); it != frag4.atoms.end(); ++it)
            cout << it->first << " ";
    cout << endl;*/

    // adjust bond length to reflect single bond characteristic
    bond_vector = atoms[bond.second].coordinates - atoms[bond.first].coordinates;
    bond_vector.normalize();
    delta = library.length(atoms[bond.second].element, atoms[bond.first].element) * bond_vector;
    atoms[bond.second].coordinates = atoms[bond.first].coordinates + delta;
    // template of a hydrogen
    atom toAdd;
    toAdd.element = hydrogen;
    toAdd.name = hydrogen + "  ";
    toAdd.ID = cur_atom->ID;

    // need to produce 3 bonds (indexed bond.first)
    // extract connected atom and find a reference atom
    cur_atom = &atoms[bond.first];
    atom nearAtom = atoms[bond.second];
    // calculate the normal to rotate
    v1 = atoms[connect3].coordinates - nearAtom.coordinates;
    v2 = cur_atom->coordinates - nearAtom.coordinates;
    normal = -v2^v1;
    if (normal == Vec3d()) return -1;
    // produce 2 siblings, 120 degrees apart
    rot = rot.createRotation(2 * M_PI / 3, v2);
    v3 = rot*normal;
    rot = rot.createRotation(-2 * M_PI / 3, v2);
    delta = rot*normal;
    // rotate to correct orientation (19.5 degrees), produce 3
    rot = rot.createRotation(19.5 / 180 * M_PI, normal^v2);
    normal = rot*normal;
    normal.normalize();
    // calculate the new delta to the fragment
    v1 = library.length(cur_atom->element, atoms[connect1].element) * normal;
    frag1.Translate(connect1, cur_atom->coordinates + v1);
    // second bond
    rot = rot.createRotation(19.5 / 180 * M_PI, v3^v2);
    v3 = rot*v3;
    v3.normalize();
    v1 = library.length(cur_atom->element, atoms[connect2].element) * v3;
    frag2.Translate(connect2, cur_atom->coordinates + v1);
    // third bond (hydrogen)
    rot = rot.createRotation(19.5 / 180 * M_PI, delta^v2);
    delta = rot*delta;
    delta.normalize();
    v1 = library.length(cur_atom->element, hydrogen) * delta;
    toAdd.Residue = cur_atom->Residue;
    convertor.str(string());
    convertor << nextIndex;
    toAdd.PDBIndex = convertor.str();
    toAdd.IndexArray.clear();
    toAdd.IndexArray.insert(bond.first);
    toAdd.coordinates = cur_atom->coordinates + v1;
    side1.insert(pair<int, atom > (nextIndex, toAdd));
    cur_atom->IndexArray.insert(nextIndex);

    ++nextIndex;

    // need to produce 3 bonds (indexed bond.second)
    // extract connected atom and find a reference atom
    cur_atom = &atoms[bond.second];
    nearAtom = atoms[bond.first];
    // calculate the normal to rotate
    v1 = atoms[connect1].coordinates - nearAtom.coordinates;
    v2 = cur_atom->coordinates - nearAtom.coordinates;
    normal = -v2^v1;
    // produce 2 siblings, 120 degrees apart
    rot = rot.createRotation(2 * M_PI / 3, v2);
    v3 = rot*normal;
    rot = rot.createRotation(-2 * M_PI / 3, v2);
    delta = rot*normal;
    // rotate to correct orientation (19.5 degrees), produce 3
    rot = rot.createRotation(19.5 / 180 * M_PI, normal^v2);
    normal = rot*normal;
    normal.normalize();
    // calculate the new delta to the fragment
    v1 = library.length(cur_atom->element, atoms[connect3].element) * normal;
    frag3.Translate(connect3, cur_atom->coordinates + v1);
    // second bond
    rot = rot.createRotation(19.5 / 180 * M_PI, v3^v2);
    v3 = rot*v3;
    v3.normalize();
    v1 = library.length(cur_atom->element, atoms[connect4].element) * v3;
    frag4.Translate(connect4, cur_atom->coordinates + v1);
    // third bond (hydrogen)
    rot = rot.createRotation(19.5 / 180 * M_PI, delta^v2);
    delta = rot*delta;
    delta.normalize();
    v1 = library.length(cur_atom->element, hydrogen) * delta;
    toAdd.Residue = cur_atom->Residue;
    convertor.str(string());
    convertor << nextIndex;
    toAdd.PDBIndex = convertor.str();
    toAdd.IndexArray.clear();
    toAdd.IndexArray.insert(bond.second);
    toAdd.coordinates = cur_atom->coordinates + v1;
    side2.insert(pair<int, atom > (nextIndex, toAdd));
    cur_atom->IndexArray.insert(nextIndex);

    side1.insert(pair<int, atom > (bond.first, atoms[bond.first]));
    side2.insert(pair<int, atom > (bond.second, atoms[bond.second]));

    // rotate fragments to reflect best orientation
    v1 = atoms[connect1].coordinates - atoms[bond.first].coordinates;
    v2 = frag1.atoms[connect1].coordinates - atoms[bond.first].coordinates;
    normal = v2^v1;
    angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
    frag1.RotateLine(normal, connect1, angle);

    v1 = atoms[connect2].coordinates - atoms[bond.first].coordinates;
    v2 = frag2.atoms[connect2].coordinates - atoms[bond.first].coordinates;
    normal = v2^v1;
    angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
    frag2.RotateLine(normal, connect2, angle);

    v1 = atoms[connect3].coordinates - atoms[bond.second].coordinates;
    v2 = frag3.atoms[connect3].coordinates - atoms[bond.second].coordinates;
    normal = v2^v1;
    angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
    frag3.RotateLine(normal, connect3, angle);

    v1 = atoms[connect4].coordinates - atoms[bond.second].coordinates;
    v2 = frag4.atoms[connect4].coordinates - atoms[bond.second].coordinates;
    normal = v2^v1;
    angle = acos(v1 * v2 / sqrt(v1.length2() * v2.length2()));
    frag4.RotateLine(normal, connect4, angle);
    // recontruct the whole molecule in half
    for (map<int, atom>::iterator it = frag1.atoms.begin(); it != frag1.atoms.end(); ++it)
        side1.insert(*it);
    for (map<int, atom>::iterator it = frag2.atoms.begin(); it != frag2.atoms.end(); ++it)
        side1.insert(*it);
    for (map<int, atom>::iterator it = frag3.atoms.begin(); it != frag3.atoms.end(); ++it)
        side2.insert(*it);
    for (map<int, atom>::iterator it = frag4.atoms.begin(); it != frag4.atoms.end(); ++it)
        side2.insert(*it);

    // rotate to minimize steric hinderance
    angle = DihedralAngle(side1[connect1], side1[bond.first], side2[bond.second], side2[connect3]);
    frag1.atoms = side1;
    frag2.atoms = side2;
    frag2.RotateLine(bond_vector, bond.second, angle + M_PI / 3);
    nextIndex = 0;
    max_dist = 0;
    for (int i = 0; i != 3; ++i) {
        frag2.RotateLine(bond_vector, bond.second, 2 * M_PI / 3);
        dist = frag1.MolecularDistance(frag2);
        if (dist > max_dist) {
            nextIndex = i;
            max_dist = dist;
        }
    }
    atoms.clear();
    frag2.RotateLine(bond_vector, bond.second, 2 * M_PI / 3 * nextIndex + M_PI / 3);
    for (map<int, atom>::iterator it = frag1.atoms.begin(); it != frag1.atoms.end(); ++it)
        atoms.insert(*it);
    for (map<int, atom>::iterator it = frag2.atoms.begin(); it != frag2.atoms.end(); ++it)
        atoms.insert(*it);
    return 0;
}

// recursively trace all atoms in the branch

void Ligand::split_fragment(Ligand& target, int start, int end)
{
    atom a = atoms[start];
    if (!(target.atoms.insert(pair<int, atom > (start, a))).second)
        return;
    for (set<int>::iterator it = a.IndexArray.begin(); it != a.IndexArray.end(); ++it)
        if (*it != end)
            split_fragment(target, *it, end);
}

inline bool Ligand::checkHeteroAtom(list<int>& candidate)
{
    // ensure the whole thing is connected
    bool flag = true;
    list<int>::iterator list_it, list_it2;
    string carbon("C");
    atom* cur_atom;
    int notcarbon(0), count;
    for (list_it = candidate.begin(); list_it != candidate.end(); ++list_it) {
        count = 0;
        cur_atom = &atoms[*list_it];
        // a heterocyclic ring will not contains too many heteroatom
        if (cur_atom->element != carbon)
            ++notcarbon;
        cur_atom = &atoms[*list_it];
        for (list_it2 = candidate.begin(); list_it2 != candidate.end(); ++list_it2)
            if (cur_atom->IndexArray.find(*list_it2) != cur_atom->IndexArray.end())
                ++count;
        // an atom is connected to exactly 2 others in this ring
        if (count != 2)
            flag = false;
    }
    //if (notcarbon > 2)
    //	flag = false;
    return flag;
}

// determine whether both molecules consist ring structure that can be joined

int Ligand::JoinRing(Ligand& ref)
{
    int ligandRingCount, fragmentRingCount, count(0), step(0);
    list<int> ring;
    ligandRingCount = DetectRing(ring);
    // no ring structure in either molecule, pointless to execute this method
    if (ligandRingCount == 0)
        return -1;
    // find out a segment of ring that has no other connection, consecutive atoms
    // first: index of atom, second: frequency
    map<int, int> ligandPairs, fragmentPairs;
    map<int, int>::iterator map_it;
    atom cur_atom;
    // scan through all rings of the ligand, find out those occurred once
    for (int i = 0; i != ligandRingCount; ++i) {
        DetectRing(ring, i);
        for (list<int>::iterator it = ring.begin(); it != ring.end(); ++it) {
            map_it = ligandPairs.find(*it);
            if (map_it == ligandPairs.end())
                ligandPairs.insert(pair<int, int>(*it, 1));
            else
                map_it->second++;
        }
    }
    fragmentRingCount = ref.DetectRing(ring);
    if (fragmentRingCount == 0)
        return -1;
    // scan through all rings of the fragment, find out those occurred once
    for (int i = 0; i != fragmentRingCount; ++i) {
        ref.DetectRing(ring, i);
        for (list<int>::iterator it = ring.begin(); it != ring.end(); ++it) {
            map_it = fragmentPairs.find(*it);
            if (map_it == fragmentPairs.end())
                fragmentPairs.insert(pair<int, int>(*it, 1));
            else
                map_it->second++;
        }
    }
    // remove those indice that have occurred more than once
    map_it = ligandPairs.begin();
    while (map_it != ligandPairs.end())
        if (map_it->second != 1)
            ligandPairs.erase(map_it++);
        else
            ++map_it;
    map_it = fragmentPairs.begin();
    while (map_it != fragmentPairs.end())
        if (map_it->second != 1)
            fragmentPairs.erase(map_it++);
        else
            ++map_it;
    if (ligandPairs.empty() || fragmentPairs.empty()) return -1;
    // randomly pick an edge form ligand and fragment
    pair<int, int> edgeLigand, edgeFragment;
    bool ready = false;
    int ligandRemoveIndex1(-1), ligandRemoveIndex2(-1), fragmentRemoveIndex1(-1), fragmentRemoveIndex2(-1), ligandRefIndex, fragmentRefIndex;
    string hydrogen("H"), ligandElement1, ligandElement2, fragmentElement1, fragmentElement2;

    while (count != 100 && !ready) {
        // pick a pair from the ligand
        map_it = ligandPairs.begin();
//        step = generator_int() % ligandPairs.size();
	step = 0 % ligandPairs.size();
        advance(map_it, step);
        edgeLigand.first = map_it->first;
        cur_atom = atoms[edgeLigand.first];
        // from connected atoms, find one that would be an isolated edge
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it) {
            if (ligandPairs.find(*it) != ligandPairs.end()) {
                edgeLigand.second = *it;
                break;
            }
        }
        // pick a pair from the fragment
        map_it = fragmentPairs.begin();
//        step = generator_int() % fragmentPairs.size();
	step = 0 % fragmentPairs.size();
        advance(map_it, step);
        edgeFragment.first = map_it->first;
        cur_atom = ref.atoms[edgeFragment.first];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it) {
            if (fragmentPairs.find(*it) != fragmentPairs.end()) {
                edgeFragment.second = *it;
                break;
            }
        }
        // sometimes the selected atom does not have a free edge
        if (edgeLigand.second == 0 || edgeFragment.second == 0) {
            ++count;
            continue;
        }
        ligandElement1 = atoms[edgeLigand.first].element;
        ligandElement2 = atoms[edgeLigand.second].element;
        fragmentElement1 = ref.atoms[edgeFragment.first].element;
        fragmentElement2 = ref.atoms[edgeFragment.second].element;
        // the edges on ligand and fragment contains different element, cannot be aligned to merge
        if ((ligandElement1 != fragmentElement1 || ligandElement2 != fragmentElement2) &&
            (ligandElement2 != fragmentElement1 || ligandElement1 != fragmentElement2)) {
            ++count;
            continue;
        }
        // the difference in edge length is too much
        if (abs((atoms[edgeLigand.first].coordinates - atoms[edgeLigand.second].coordinates).length2()
                -(ref.atoms[edgeFragment.first].coordinates - ref.atoms[edgeFragment.second].coordinates).length2()) > 0.25) {
            ++count;
            continue;
        }
        // get the index of a hydrogen to be removed from each of the atom
        cur_atom = atoms[edgeLigand.first];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
            if (atoms[*it].element == hydrogen)
                ligandRemoveIndex1 = *it;
        cur_atom = atoms[edgeLigand.second];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
            if (atoms[*it].element == hydrogen)
                ligandRemoveIndex2 = *it;
        cur_atom = ref.atoms[edgeFragment.first];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
            if (ref.atoms[*it].element == hydrogen)
                fragmentRemoveIndex1 = *it;
        cur_atom = ref.atoms[edgeFragment.second];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
            if (ref.atoms[*it].element == hydrogen)
                fragmentRemoveIndex2 = *it;
        // some more advanced structure is being connected without a hydrogen found
        if (ligandRemoveIndex1 == -1 || ligandRemoveIndex2 == -1 || fragmentRemoveIndex1 == -1 || fragmentRemoveIndex2 == -1) {
            ++count;
            continue;
        }
        ready = true;
    }
    //if (atoms.begin()->first == 0) makeLog("invalid ligand in ring joining");
    // ring joining is not possible
    if (!ready) return -1;

    // 80% chance of not using this method
//    if (generator() > 0.8) return -1;

    // remove hydrogen atoms
    DeleteAtom(ligandRemoveIndex1);
    DeleteAtom(ligandRemoveIndex2);
    ref.DeleteAtom(fragmentRemoveIndex1);
    ref.DeleteAtom(fragmentRemoveIndex2);

    Vec3d v1, v2, normal;
    double angle, dist;
    // now use this flag to indicate how the atom elements corresponded
    ready = true;
    if (atoms[edgeLigand.first].element != ref.atoms[edgeFragment.first].element)
        ready = false;
    // move the fragment to the correct coordinates on the ligand
    v1 = atoms[edgeLigand.second].coordinates - atoms[edgeLigand.first].coordinates;
    if (ready) {
        ref.Translate(edgeFragment.first, atoms[edgeLigand.first].coordinates);
        v2 = ref.atoms[edgeFragment.second].coordinates - ref.atoms[edgeFragment.first].coordinates;
        // get the roration vector between the fragment and ligand
        normal = v2^v1;
        // compute angle between fragment and ligand
        angle = Angle(atoms[edgeLigand.second].coordinates, atoms[edgeLigand.first].coordinates, ref.atoms[edgeFragment.second].coordinates);
        // record the old distance
        dist = (ref.atoms[edgeFragment.second].coordinates - atoms[edgeLigand.second].coordinates).length();
        if (angle != 0)
            ref.RotateLine(normal, edgeFragment.first, angle);
        // rotate the other side if the distance has increased
        if ((ref.atoms[edgeFragment.second].coordinates - atoms[edgeLigand.second].coordinates).length() > dist)
            ref.RotateLine(normal, edgeFragment.first, -2 * angle);
        // find out another adjacent atom
        cur_atom = ref.atoms[edgeFragment.first];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
            if (ref.atoms[*it].element != hydrogen && *it != edgeFragment.second)
                fragmentRefIndex = *it;
    }
    else {
        ref.Translate(edgeFragment.second, atoms[edgeLigand.first].coordinates);
        v2 = ref.atoms[edgeFragment.first].coordinates - ref.atoms[edgeFragment.second].coordinates;
        normal = v2^v1;
        angle = Angle(atoms[edgeLigand.second].coordinates, atoms[edgeLigand.first].coordinates, ref.atoms[edgeFragment.first].coordinates);
        dist = (ref.atoms[edgeFragment.first].coordinates - atoms[edgeLigand.second].coordinates).length();
        if (angle != 0)
            ref.RotateLine(normal, edgeFragment.second, angle);
        if ((ref.atoms[edgeFragment.first].coordinates - atoms[edgeLigand.second].coordinates).length() > dist)
            ref.RotateLine(normal, edgeFragment.second, -2 * angle);
        cur_atom = ref.atoms[edgeFragment.second];
        for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
            if (ref.atoms[*it].element != hydrogen && *it != edgeFragment.first)
                fragmentRefIndex = *it;
    }
    // make a reference point on the ligand
    cur_atom = atoms[edgeLigand.second];
    for (set<int>::iterator it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
        if (atoms[*it].element != hydrogen && *it != edgeLigand.first)
            ligandRefIndex = *it;
    angle = DihedralAngle(atoms[ligandRefIndex], atoms[edgeLigand.second], atoms[edgeLigand.first], ref.atoms[fragmentRefIndex]);
    // try to planarise the rings
    // record max distance and its orientation in a map
    // try all 4 orientations
    map<double, double> best_angles;
    ref.RotateLine(v1, edgeFragment.first, -angle);
    dist = MolecularDistance(ref);
    best_angles.insert(pair<double, double>(dist, -angle));
    ref.RotateLine(v1, edgeFragment.first, M_PI);
    dist = MolecularDistance(ref);
    best_angles.insert(pair<double, double>(dist, M_PI - angle));
    ref.RotateLine(v1, edgeFragment.first, -M_PI + 2 * angle);
    dist = MolecularDistance(ref);
    best_angles.insert(pair<double, double>(dist, angle));
    ref.RotateLine(v1, edgeFragment.first, M_PI);
    dist = MolecularDistance(ref);
    best_angles.insert(pair<double, double>(dist, M_PI + angle));
    ref.RotateLine(v1, edgeFragment.first, -M_PI - angle);
    ref.RotateLine(v1, edgeFragment.first, best_angles.rbegin()->second);

    // redirect atom index to ligand
    int updateIndex, cascadeIndex(MaxIndex());
    set<int> tempIndice;
    ostringstream output;
    for (map<int, atom>::iterator it = ref.atoms.begin(); it != ref.atoms.end(); ++it) {
        // the edge on the fragment is being discarded, the originally connected atoms are now redirected onto the ligand
        if (it->first == edgeFragment.first || it->first == edgeFragment.second) continue;
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
        for (set<int>::iterator iter = toAdd.IndexArray.begin(); iter != toAdd.IndexArray.end(); ++iter) {
            if (*iter == edgeFragment.first) {
                // first-to-first
                if (ready) {
                    tempIndice.insert(edgeLigand.first);
                    // add the connect index to the ligand
                    atoms[edgeLigand.first].IndexArray.insert(updateIndex);
                }
                else {
                    tempIndice.insert(edgeLigand.second);
                    atoms[edgeLigand.second].IndexArray.insert(updateIndex);
                }
            }
            else if (*iter == edgeFragment.second) {
                if (ready) {
                    tempIndice.insert(edgeLigand.second);
                    atoms[edgeLigand.second].IndexArray.insert(updateIndex);
                }
                else {
                    tempIndice.insert(edgeLigand.first);
                    atoms[edgeLigand.first].IndexArray.insert(updateIndex);
                }
            }
            else {
                tempIndice.insert(*iter + cascadeIndex);
            }
        }
        // clear the original indice
        toAdd.IndexArray.clear();
        // move the indice to the original set
        for (set<int>::iterator iter = tempIndice.begin(); iter != tempIndice.end(); ++iter)
            toAdd.IndexArray.insert(*iter);
        // add atom to the current molecule
        atoms.insert(pair<int, atom > (updateIndex, toAdd));
    }
    return 0;
}
