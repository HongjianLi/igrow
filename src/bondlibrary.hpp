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

#ifndef IGROW_BONDLIBRARY_HPP
#define IGROW_BONDLIBRARY_HPP

#include "common.hpp"

namespace igrow
{

    class bond_library
    {
    public:

        static const fl BADBOND_THRESHOLD;
        
        // enumeration of bond types

        enum BOND_TYPE
        {
            SINGLE_BOND = 1,
            DOUBLE_BOND = 2,
            TRIPLE_BOND = 3,
        };

        // return the bond length between 2 given element
        double length(string element1, string element2);

        // an unused function
        bool badBond(string element1, string element2, double test_length);

        // determine type of bond using atomic distance
        BOND_TYPE type(string element1, string element2, double test_length);
        
        // return the molecular weight of an element
        double weight(string element);
        
        // todo: add a method by taking 2 atoms
        
    private:
        // private methods in determining bond length
        inline double checkHydrogen(string element);
        inline double checkCarbon(string element);
        inline double checkNitrogen(string element);
        inline double checkOxygen(string element);
        inline double checkFlorine(string element);
        inline double checkSulphur(string element);
        inline double checkPhosphorus(string element);
        inline double checkChlorine(string element);
        inline double checkBromine(string element);
        inline double checkIodine(string element);
    };

}

#endif
