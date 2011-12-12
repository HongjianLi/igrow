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

#include "bondlibrary.hpp"

using namespace std;

namespace igrow
{

    double bond_library::length(string element1, string element2)
    {
	bool multiwords = false;
	if (element1.length() > 1 && isalpha(element1[1]))
	    multiwords = true;
	if (!multiwords && element1[0] == 'C')
	    return checkCarbon(element2);
	if (element1[0] == 'N')
	    return checkNitrogen(element2);
	if (element1[0] == 'O')
	    return checkOxygen(element2);
	if (element1[0] == 'H')
	    return checkHydrogen(element2);
	if (element1[0] == 'P')
	    return checkPhosphorus(element2);
	if (element1[0] == 'S')
	    return checkSulphur(element2);
	if (element1[0] == 'F')
	    return checkFlorine(element2);
	if (multiwords && element1[0] == 'C' && element1[0] == 'L')
	    return checkChlorine(element2);
	if (multiwords && element1[0] == 'B' && element1[0] == 'R')
	    return checkBromine(element2);
	if (element1[0] == 'I')
	    return checkIodine(element2);
	return 0;
    }

    bool bond_library::badBond(string element1, string element2, double test_length)
    {
	double length_ = length(element1, element2);
	if (test_length - length_ > BADBOND_THRESHOLD || length_ - test_length > BADBOND_THRESHOLD)
	    return false;
	return true;
    }

    bond_library::BOND_TYPE bond_library::type(string element1, string element2, double test_length)
    {
	string carbon("C"), nitrogen("N"), oxygen("O"), phosphorus("P");
	// linear separation of atomic distance to determine bond type
	if (element1 == carbon && element2 == carbon)
	{
	    if (test_length > 1.540)
		return SINGLE_BOND;
	    if (test_length < 1.200)
		return TRIPLE_BOND;
	    if (test_length > 1.340)
	    {
		if (1.540 - test_length < test_length - 1.340)
		    return SINGLE_BOND;
		else
		    return DOUBLE_BOND;
	    }
	    else
	    {
		if (1.340 - test_length < test_length - 1.200)
		    return DOUBLE_BOND;
		else
		    return TRIPLE_BOND;
	    }
	}
	if ((element1 == carbon && element2 == nitrogen) || (element1 == nitrogen && element2 == carbon))
	{
	    if (test_length > 1.465)
		return SINGLE_BOND;
	    if (test_length < 1.136)
		return TRIPLE_BOND;
	    if (test_length > 1.340)
	    {
		if (1.465 - test_length < test_length - 1.279)
		    return SINGLE_BOND;
		else
		    return DOUBLE_BOND;
	    }
	    else
	    {
		if (1.279 - test_length < test_length - 1.136)
		    return DOUBLE_BOND;
		else
		    return TRIPLE_BOND;
	    }
	}
	if ((element1 == carbon && element2 == oxygen) || (element1 == oxygen && element2 == carbon))
	{
	    /*if (test_length > 1.413)
		    return SINGLE_BOND;
	    if (test_length < 1.230)
		    return DOUBLE_BOND;*/
	    // partial double could probably classified as single bond
	    if (test_length > 1.322)
		return SINGLE_BOND;
	    else
		return DOUBLE_BOND;
	}
	if ((element1 == carbon && element2 == phosphorus) || (element1 == phosphorus && element2 == carbon))
	{
	    if (test_length > 1.621)
		return SINGLE_BOND;
	    else
		return DOUBLE_BOND;
	}
	return SINGLE_BOND;
    }

    double bond_library::weight(string element)
    {
	bool multiwords = false;
	if (element.length() > 1 && isalpha(element[1]))
	    multiwords = true;
	if (multiwords)
	{
	    if (element[0] == 'C' && element[1] == 'L')
		return 35.453;
	    if (element[0] == 'B' && element[1] == 'R')
		return 79.904;
	}
	else
	{
	    switch (element[0])
	    {
		case 'H':
		    return 1.0079;
		case 'C':
		    return 12.011;
		case 'N':
		    return 14.007;
		case 'O':
		    return 15.999;
		case 'F':
		    return 18.998;
		case 'P':
		    return 30.974;
		case 'S':
		    return 32.066;
		case 'I':
		    return 126.90;
		default:
		    return 1;
	    }
	}
	return 1;
    }

    // statements are arranged in accordance to occurrence

    inline double bond_library::checkHydrogen(string element)
    {
	if (element[0] == 'C')
	    return 1.059;
	if (element[0] == 'N')
	    return 1.009;
	if (element[0] == 'O')
	    return 0.967;
	if (element[0] == 'S')
	    return 1.013;
	return 0;
    }

    inline double bond_library::checkCarbon(string element)
    {
	bool multiwords = false;
	if (element.length() > 1 && isalpha(element[1]))
	    multiwords = true;
	if (!multiwords && element[0] == 'C')
	    return 1.530;
	if (element[0] == 'H')
	    return 1.059;
	if (element[0] == 'O')
	    return 1.413;
	if (element[0] == 'N')
	    return 1.469;
	if (element[0] == 'S')
	    return 1.819;
	if (element[0] == 'F')
	    return 1.399;
	if (element[0] == 'P')
	    return 1.85;
	if (multiwords && element[0] == 'C' && element[1] == 'L')
	    return 1.790;
	if (multiwords && element[0] == 'B' && element[1] == 'R')
	    return 1.910;
	if (element[0] == 'I')
	    return 2.162;
	return 0;
    }

    inline double bond_library::checkNitrogen(string element)
    {
	if (element[0] == 'C')
	    return 1.469;
	if (element[0] == 'N')
	    return 1.425;
	if (element[0] == 'H')
	    return 1.009;
	if (element[0] == 'O')
	    return 1.463;
	return 0;
    }

    inline double bond_library::checkOxygen(string element)
    {
	if (element[0] == 'C')
	    return 1.413;
	if (element[0] == 'H')
	    return 0.967;
	if (element[0] == 'N')
	    return 1.463;
	if (element[0] == 'O')
	    return 1.469;
	if (element[0] == 'S')
	    return 1.577;
	if (element[0] == 'P')
	    return 1.540;
	return 0;
    }

    inline double bond_library::checkFlorine(string element)
    {
	if (element[0] == 'C')
	    return 1.399;
	if (element[0] == 'P')
	    return 1.560;
	return 0;
    }

    inline double bond_library::checkSulphur(string element)
    {
	if (element[0] == 'C')
	    return 1.819;
	if (element[0] == 'O')
	    return 1.577;
	if (element[0] == 'S')
	    return 2.048;
	if (element[0] == 'H')
	    return 1.013;
	if (element[0] == 'P')
	    return 2.090;
	return 0;
    }

    inline double bond_library::checkPhosphorus(string element)
    {
	bool multiwords = false;
	if (element.length() > 1 && isalpha(element[1]))
	    multiwords = true;
	if (element[0] == 'F')
	    return 1.560;
	if (multiwords && element[0] == 'C' && element[1] == 'L')
	    return 2.040;
	if (multiwords && element[0] == 'B' && element[1] == 'R')
	    return 2.220;
	if (element[0] == 'I')
	    return 2.430;
	if (element[0] == 'S')
	    return 2.090;
	if (element[0] == 'P')
	    return 2.235;
	if (element[0] == 'O')
	    return 1.540;
	if (element[0] == 'C')
	    return 1.85;
	// double bond has length 1.480
	return 0;
    }

    inline double bond_library::checkChlorine(string element)
    {
	if (element[0] == 'C')
	    return 1.790;
	if (element[0] == 'P')
	    return 2.040;
	return 0;
    }

    inline double bond_library::checkBromine(string element)
    {
	if (element[0] == 'C')
	    return 1.910;
	if (element[0] == 'P')
	    return 2.220;
	return 0;
    }

    inline double bond_library::checkIodine(string element)
    {
	if (element[0] == 'C')
	    return 2.162;
	if (element[0] == 'P')
	    return 2.430;
	return 0;
    }

}
