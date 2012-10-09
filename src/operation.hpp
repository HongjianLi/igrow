/*

   Copyright (c) 2012, The Chinese University of Hong Kong

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
#ifndef IGROW_OPERATION_HPP
#define IGROW_OPERATION_HPP

#include <vector>
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/atomic.hpp>
#include "ligand.hpp"

namespace igrow
{
	using std::vector;
	using boost::ptr_vector;
	using boost::filesystem::path;

	/// Represents a GA operation, be it either mutation or crossover.
	class operation
	{
	public:
		/// Constructs a GA operation.
		explicit operation(ptr_vector<ligand>& ligands, const size_t num_elitists, const vector<path>& fragments, const validator& v, const size_t max_failures, boost::atomic<size_t>& num_failures) : ligands(ligands), num_elitists(num_elitists), fragments(fragments), num_fragments(fragments.size()), v(v), max_failures(max_failures), num_failures(num_failures) {}

		/// Task for creating a child ligand from two parent ligands by mutation.
		/// @exception maximum_failures_reached_error Thrown when the number of failures reaches the user specified maximum number of failures.
		void mutation_task(const size_t index, const path& p, const size_t seed);

		/// Task for creating a child ligand from one parent ligands by splitting.
		/// @exception maximum_failures_reached_error Thrown when the number of failures reaches the user specified maximum number of failures.
		void split_task(const size_t index, const path& p, const size_t seed);

		/// Task for creating a child ligand from two parent ligands by crossover.
		/// @exception maximum_failures_reached_error Thrown when the number of failures reaches the user specified maximum number of failures.
		void crossover_task(const size_t index, const path& p, const size_t seed);

	protected:
		ptr_vector<ligand>& ligands;
		const size_t num_elitists;
		const vector<path>& fragments;
		const size_t num_fragments;
		const validator& v;
		const size_t max_failures;
		boost::atomic<size_t>& num_failures;

//		/// Represent a GA operation for removal of ligand duplicates.
//		class operation_code
//		{
//		public:
//			size_t operation; ///< 0 for mutation and 1 for crossover.
//			path parent1; ///< Parent ligand 1.
//			path parent2; ///< Parent ligand 2.
//			size_t connector1; ///< Serial number of the connecting atom of parent 1.
//			size_t connector2; ///< Serial number of the connecting atom of parent 2.
//
//			/// Constructs an operation code, which is uniquely identified by the type of GA operation, the two parent ligands and the two connector atoms. Note that operation is required because both mutation and crossover can occur on the same parent ligands and connector atoms but result in different child ligands.
//			operation_code(const size_t operation, const path& parent1, const path& parent2, const size_t connector1, const size_t connector2) : operation(operation), parent1(parent1), parent2(parent2), connector1(connector1), connector2(connector2) {}
//		};
	};

	class maximum_failures_reached_error : public std::domain_error
	{
	public:
		maximum_failures_reached_error(const size_t max_failures) : std::domain_error("The number of failures has reached " + lexical_cast<string>(max_failures)) {}
	};
}

#endif
