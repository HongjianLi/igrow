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

#ifndef IGROW_COMMON_HPP
#define IGROW_COMMON_HPP

#include <vector>
#include <string>
#include <boost/assert.hpp>

namespace igrow
{
	// These classes are widely used across the entire program.
	using std::runtime_error;
	using std::vector;
	using std::string;

	/// igrow uses double precision floating point computation by default.
	/// This could possible be demoted to single precision for better performance.
	typedef double fl;
}

#endif
