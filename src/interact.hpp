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

#ifndef IGROW_INTERACT_HPP
#define IGROW_INTERACT_HPP

#include "common.hpp"
#include "ligand.hpp"
#include <map>

namespace igrow
{

    class Interaction
    {
        // set of atoms that have been visited
        std::set<int> scanned;
        // pairs of atoms that are considered to be equal in terms of type and position
        std::map<int, int> overlap;
    public:


        Ligand* mate(std::string Filename1, std::string Filename2);
        // selectively combine two molecules to produce a new one
        Ligand* mate(Ligand* male, Ligand* female);
        Ligand* merge(std::string Filename1, std::string Filename2);
        // maximally merge two molecules to produce a new one
        Ligand* merge(Ligand* male, Ligand* female);
    protected:
        // traverse the molecule and decide which atom to select
        void scan_recursive(Ligand* ref, int index);
        // traverse the molecule to find out which atom to be removed
        void remove_invalid(Ligand* ref, int index);
    };

}

#endif
