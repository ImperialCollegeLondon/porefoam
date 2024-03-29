/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
	This file is part of foam-extend.

	foam-extend is free software: you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published by the
	Free Software Foundation, either version 3 of the License, or (at your
	option) any later version.

	foam-extend is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
	Foam::pointMeshMapper

Description
	Class holds all the necessary information for mapping fields associated
	with pointMesh.

SourceFiles
	pointMeshMapper.C

\*---------------------------------------------------------------------------*/

#ifndef pointMeshMapper_H
#define pointMeshMapper_H

#include "mapPolyMesh.H"
#include "pointMapper.H"
#include "pointBoundaryMeshMapper.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class pointMesh;
class mapPolyMesh;


class pointMeshMapper
{
	// Private data

		//- Reference to mesh
		const pointMesh& mesh_;

		//- Point mapper
		pointMapper pointMap_;

		//- Boundary mapper
		pointBoundaryMeshMapper boundaryMap_;

		//- Reset patch flag.  Set true for the rebuilt patch (no mapping)
		boolList resetPatchFlag_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		pointMeshMapper(const pointMeshMapper&);

		//- Disallow default bitwise assignment
		void operator=(const pointMeshMapper&);


public:

	// Constructors

		//- Construct from pointMesh
		pointMeshMapper(const pointMesh& mesh, const mapPolyMesh& mpm)
		:
			mesh_(mesh),
			pointMap_(mpm),
			boundaryMap_(mesh, pointMap_, mpm),
			resetPatchFlag_(mpm.resetPatchFlag())
		{}


	// Member Functions


		//- Return reference to mesh fields belong to
		const pointMesh& mesh() const
		{
			return mesh_;
		}

		//- Return reference to objectRegistry storing fields. Can be
		//  removed once fields stored on pointMesh.
		const objectRegistry& thisDb() const
		{
			return mesh_();
		}

		//- Return point mapper
		const morphFieldMapper& pointMap() const
		{
			return pointMap_;
		}

		//- Return boundary mapper
		const pointBoundaryMeshMapper& boundaryMap() const
		{
			return boundaryMap_;
		}

		//- Return reset patch flag
		const boolList& resetPatchFlag() const
		{
			return resetPatchFlag_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
