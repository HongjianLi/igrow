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

#pragma once
#ifndef IGROW_MUTATION_TASK_HPP
#define IGROW_MUTATION_TASK_HPP

#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/path.hpp>
#include "ligand.hpp"

namespace igrow
{
	using std::vector;
	using boost::ptr_vector;
	using boost::filesystem::path;
	
	/// Task for creating a child ligand from two parent ligands by mutation.
	int mutation_task(ptr_vector<ligand>& ligands, const size_t index, const path& filename, const size_t num_elitists, const vector<path>& fragments, const validator& v, const size_t seed, const size_t max_failures, size_t num_failures);
}

#endif
