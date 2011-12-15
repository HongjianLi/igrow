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

#include "ligand.hpp"

namespace igrow
{
    using boost::filesystem::path;

    class Interaction
    {
    public:

        ligand mate(ligand male, ligand female);
        // maximally merge two molecules to produce a new one
        ligand merge(ligand male, ligand female);

    protected:
        // traverse the molecule and decide which atom to select
        void scan_recursive(ligand ref, int index);
        // traverse the molecule to find out which atom to be removed
        void remove_invalid(ligand ref, int index);
        
    private:
        // set of atoms that have been visited
        std::set<int> scanned;
        
        // pairs of atoms that are considered to be equal in terms of type and position
        std::map<int, int> overlap;
    };

}

#endif
