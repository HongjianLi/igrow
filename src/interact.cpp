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

#include "interact.hpp"
#include <sstream>

namespace igrow
{
    using std::pair;
    using std::map;
    using std::set;
    using std::ostringstream;

    // overloaded method to accept file name

    ligand Interaction::mate(const path& file1, const path& file2)
    {
	ligand male;
	male.load(file1);
	ligand female;
	female.load(file2);
	return mate(male, female);
    }

    ligand Interaction::mate(ligand male, ligand female)
    {
	// since index starts with 1, just add to it
	int updateIndex, cascadeIndex = male.MaxIndex();
	set<int> tempIndice;
	ostringstream output;
	// append molecule to another by cascading the index
	for (map<int, atom>::iterator it = female.atoms.begin(); it != female.atoms.end(); ++it)
	{
	    atom toAdd = atom(it->second);
	    updateIndex = atoi(toAdd.PDBIndex.c_str());
	    updateIndex += cascadeIndex;
	    output.str(string());
	    output << updateIndex;
	    toAdd.PDBIndex = output.str();
	    tempIndice.clear();
	    for (set<int>::iterator iter = toAdd.IndexArray.begin(); iter != toAdd.IndexArray.end(); ++iter)
		tempIndice.insert(*iter + cascadeIndex);
	    toAdd.IndexArray.clear();
	    for (set<int>::iterator iter = tempIndice.begin(); iter != tempIndice.end(); ++iter)
		toAdd.IndexArray.insert(*iter);
	    male.atoms.insert(pair<int, atom > (updateIndex, toAdd));
	}

	// initialise sets
	set<int> toRemove;
	toRemove.clear();
	overlap.clear();
	scanned.clear();

	// obtain overlapping
	map<int, atom>::iterator it1, it2;
	atom *atom1, *atom2;
	// prevent comparison of atom in same molecule
	it1 = male.atoms.begin();
	// there are cascadeIndex number of atoms in first molecule
	while (it1 != male.atoms.end())
	{
	    // move iterator to appended molecule
	    it2 = it1;
	    ++it2;
	    // scan till the end, only cascadeIndex to the end number of atoms in second molecule
	    while (it2 != male.atoms.end())
	    {
		// skip index that can be found in the first molecule
		if (it2->first <= cascadeIndex)
		{
		    ++it2;
		    continue;
		}
		atom1 = &it1->second;
		atom2 = &it2->second;
		// there could only be 1 equivalent atom
		if (atom1->name == atom2->name && atom1->coordinates == atom2->coordinates)
		{
		    for (set<int>::iterator it = atom2->IndexArray.begin(); it != atom2->IndexArray.end(); ++it)
		    {
			atom1->IndexArray.insert(*it);
			male.atoms[*it].IndexArray.insert(it1->first);
		    }
		    for (set<int>::iterator it = atom1->IndexArray.begin(); it != atom1->IndexArray.end(); ++it)
		    {
			atom2->IndexArray.insert(*it);
			male.atoms[*it].IndexArray.insert(it2->first);
		    }
		    // register equivalent pair
		    overlap.insert(pair<int, int>(it1->first, it2->first));
		    overlap.insert(pair<int, int>(it2->first, it1->first));
		    // remove larger indexed one
		    toRemove.insert(it2->first);
		    break;
		}
		++it2;
	    }
	    ++it1;
	}

	// select some atoms to keep
	scan_recursive(male, male.MinIndex());

	// removal step
	map<int, atom>::reverse_iterator r_it = male.atoms.rbegin();
	while (r_it != male.atoms.rend())
	{
	    // atom not scanned
	    if (scanned.find(r_it->first) == scanned.end())
	    {
		r_it = male.DeleteAtom(r_it->first);
		// repeated atom
	    }
	    else if (toRemove.find(r_it->first) != toRemove.end())
	    {
		r_it = male.DeleteAtom(r_it->first);
	    }
	    else
	    {
		++r_it;
	    }
	}
	return male;
    }

    // traverse the molecule graph

    void Interaction::scan_recursive(ligand ref, int index)
    {
	scanned.insert(index);
	// make its counterpart scanned, prevent double counting
	if (overlap.find(index) != overlap.end())
	    scanned.insert(overlap[index]);

	set<int> IDs;
	set<int>::iterator it;
	IDs.clear();
	atom connect_atom, cur_atom = ref.atoms[index];
	int outgoingEdge(0);

	for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
	    // find not scanned neighbour
	    if (scanned.find(*it) == scanned.end())
		if (overlap.find(index) != overlap.end() && overlap.find(*it) == overlap.end())
		    ++outgoingEdge;

	if (outgoingEdge == 0)
	{
	    for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
		if (scanned.find(*it) == scanned.end())
		{
		    connect_atom = (*(ref.atoms.find(*it))).second;
		    if (connect_atom.ID == cur_atom.ID)
			scan_recursive(ref, *it);
		}
	    return;
	}

	// find possible IDs from connected atoms
	for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
	    if (scanned.find(*it) == scanned.end())
		IDs.insert(ref.atoms[*it].ID);

	if (IDs.empty()) return;

	it = IDs.begin();
	//	advance(it, generator_int()%IDs.size());
	advance(it, 0 % IDs.size());
	int toKeep = *it;

	// try all neighbours
	for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
	    if (scanned.find(*it) == scanned.end())
	    {
		connect_atom = (*(ref.atoms.find(*it))).second;
		// discard fragment of different ID
		if ((*(ref.atoms.find(*it))).second.ID == toKeep)
		    scan_recursive(ref, *it);
		else if (connect_atom.ID == cur_atom.ID)
		    // maintain core scaffold
		    scan_recursive(ref, *it);
	    }
    }

    ligand Interaction::merge(const path& file1, const path& file2)
    {
	ligand male;
	male.load(file1);
	ligand female;
	female.load(file2);
	return merge(male, female);
    }

    ligand Interaction::merge(ligand male, ligand female)
    {
	// since index starts with 1, just add to it
	int updateIndex, cascadeIndex = male.MaxIndex();
	set<int> tempIndice;
	ostringstream output;
	// append molecule to another by cascading the index
	for (map<int, atom>::iterator it = female.atoms.begin(); it != female.atoms.end(); ++it)
	{
	    atom toAdd = atom(it->second);
	    updateIndex = atoi(toAdd.PDBIndex.c_str());
	    updateIndex += cascadeIndex;
	    output.str(string());
	    output << updateIndex;
	    toAdd.PDBIndex = output.str();
	    tempIndice.clear();
	    for (set<int>::iterator iter = toAdd.IndexArray.begin(); iter != toAdd.IndexArray.end(); ++iter)
		tempIndice.insert(*iter + cascadeIndex);
	    toAdd.IndexArray.clear();
	    for (set<int>::iterator iter = tempIndice.begin(); iter != tempIndice.end(); ++iter)
		toAdd.IndexArray.insert(*iter);
	    male.atoms.insert(pair<int, atom > ((it->first) + cascadeIndex, toAdd));
	}

	// initialise sets
	set<int> toRemove;
	toRemove.clear();
	overlap.clear();
	scanned.clear();

	// obtain overlapping
	map<int, atom>::iterator it1, it2;
	atom *atom1, *atom2;
	// prevent comparison of atom in same molecule
	it1 = male.atoms.begin();
	// there are cascadeIndex number of atoms in first molecule
	while (it1 != male.atoms.end())
	{
	    // move iterator to appended molecule
	    it2 = it1;
	    ++it2;
	    // scan till the end, only cascadeIndex to the end number of atoms in second molecule
	    while (it2 != male.atoms.end())
	    {
		// skip index that can be found in the first molecule
		if (it2->first <= cascadeIndex)
		{
		    ++it2;
		    continue;
		}
		atom1 = &it1->second;
		atom2 = &it2->second;
		// there could only be 1 equivalent atom
		if (atom1->name == atom2->name && atom1->coordinates == atom2->coordinates)
		{
		    for (set<int>::iterator it = atom2->IndexArray.begin(); it != atom2->IndexArray.end(); ++it)
		    {
			atom1->IndexArray.insert(*it);
			male.atoms[*it].IndexArray.insert(it1->first);
		    }
		    for (set<int>::iterator it = atom1->IndexArray.begin(); it != atom1->IndexArray.end(); ++it)
		    {
			atom2->IndexArray.insert(*it);
			male.atoms[*it].IndexArray.insert(it2->first);
		    }
		    // register equivalent pair
		    overlap.insert(pair<int, int>(it1->first, it2->first));
		    overlap.insert(pair<int, int>(it2->first, it1->first));
		    // remove larger indexed one
		    toRemove.insert(it2->first);
		    break;
		}
		++it2;
	    }
	    ++it1;
	}

	// select some atoms to keep
	remove_invalid(male, male.MinIndex());

	// removal phase
	map<int, atom>::reverse_iterator r_it = male.atoms.rbegin();
	while (r_it != male.atoms.rend())
	{
	    if (scanned.find(r_it->first) == scanned.end())
	    {
		r_it = male.DeleteAtom(r_it->first);
		continue;
	    }
	    if (toRemove.find(r_it->first) != toRemove.end())
	    {
		male.DeleteAtom(r_it->first);
		continue;
	    }
	    ++r_it;
	}
	return male;
    }

    void Interaction::remove_invalid(ligand ref, int index)
    {
	scanned.insert(index);
	// also get its counterpart
	if (overlap.find(index) != overlap.end())
	    scanned.insert(overlap[index]);

	set<int> IDs;
	set<int>::iterator it;
	atom cur_atom = ref.atoms[index];
	for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
	    if (scanned.find(*it) == scanned.end())
		IDs.insert(ref.atoms[*it].ID);

	if (IDs.empty()) return;

	int toKeep(cur_atom.ID), outgoingEdge(0);

	for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
	    // find not scanned neighbour
	    if (scanned.find(*it) == scanned.end())
		if (overlap.find(index) != overlap.end() && overlap.find(*it) == overlap.end())
		    ++outgoingEdge;

	// possible collision at the same bond, choose one to stay
	if (IDs.size() > 1 && outgoingEdge > 1)
	{
	    it = IDs.begin();
	    //		advance(it, generator_int()%IDs.size());
	    advance(it, 0 % IDs.size());
	    toKeep = *it;
	}

	for (it = cur_atom.IndexArray.begin(); it != cur_atom.IndexArray.end(); ++it)
	    if (scanned.find(*it) == scanned.end())
	    {
		// go to the one which won
		if ((*(ref.atoms.find(*it))).second.ID == toKeep)
		    remove_invalid(ref, *it);
		// if there is only one fragment, always adopt it
		if (outgoingEdge == 1)
		    remove_invalid(ref, *it);
		// maintain core scaffold
		if ((*(ref.atoms.find(*it))).second.ID == cur_atom.ID)
		    remove_invalid(ref, *it);
	    }
    }

}
