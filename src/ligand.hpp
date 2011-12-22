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

#ifndef IGROW_LIGAND_HPP
#define IGROW_LIGAND_HPP

#include <boost/filesystem/path.hpp>
#include <boost/flyweight.hpp>
#include <boost/flyweight/key_value.hpp>
#include <boost/flyweight/no_tracking.hpp>
#include "atom.hpp"

namespace igrow
{
	using boost::filesystem::path;

	// Represents a ligand.

	class ligand
	{
	public:
		static const fl pi;

		const path p;
		size_t num_hb_donors; ///< Number of hydrogen bond donors.
		size_t num_hb_acceptors; ///< Number of hydrogen bond acceptors.
		fl mw; ///< Molecular weight.		
		fl logp; ///< LogP.
		fl free_energy; ///< Predicted free energy obtained by external docking.

		ligand() {}
		ligand(const path& p);
		
		void save(const path& file) const;
		ligand* mutate(const ligand& lig) const;
/*
		// add a fragment to the molecule by replacing a hydrogen in the original structure
		void mutate(ligand fragment);
		// a distance between two molecules given by the sum of minimum distance of each atom
		double MolecularDistance(ligand& other);
		// obtain the index of one hydrogen of this molecule
		int IndexOfRandomHydrogen();
		// move atom of index to origin along with its connected atoms
		void Translate(int index, Vec3d origin);
		// rotate along the line given by v1-v2 using indexed atom as pivot in terms of radian
		void RotateLine(Vec3d v1, Vec3d v2, int index, double radian);
		// rotate along the line normal using indexed atom as pivot in terms of radian
		void RotateLine(Vec3d normal, int index, double radian);
		// the dihedral angle among the four atoms
		double DihedralAngle(const Vec3d& a1, const Vec3d& a2, const Vec3d& a3, const Vec3d& a4);
*/
	};

	/// For sorting ptr_vector<ligand>.
	inline bool operator<(const ligand& a, const ligand& b)
	{
		return a.free_energy < b.free_energy;
	}

	/// For extracting path out of a ligand.
	class ligand_path_extractor
	{
	public:
		const path& operator()(const ligand& lig) const
		{
			return lig.p;
		}
	};

	/// Define flyweight type for ligand.
	using namespace boost::flyweights;
	typedef	flyweight<key_value<path, ligand, ligand_path_extractor>, no_tracking> ligand_flyweight;
}

#endif
