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
	Foam::pointPatch

Description
	Basic pointPatch represents a set of points from the mesh.

SourceFiles
	pointPatch.C

\*---------------------------------------------------------------------------*/

#ifndef pointPatch_H
#define pointPatch_H

#include "labelList.H"
#include "vectorField.H"
#include "triFaceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

class pointBoundaryMesh;
class pointConstraint;


class pointPatch
{
	// Private data

		//- Reference to boundary mesh
		const pointBoundaryMesh& boundaryMesh_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		pointPatch(const pointPatch&);

		//- Disallow default bitwise assignment
		void operator=(const pointPatch&);


protected:

	// Protected Member Functions

		// The geometry initialisation is called by pointBoundaryMesh
		friend class pointBoundaryMesh;

		//- Initialise the calculation of the patch geometry
		virtual void initGeometry()
		{}

		//- Calculate the patch geometry
		virtual void calcGeometry()
		{}

		//- Correct patches after moving points
		virtual void movePoints()
		{}

		//- Initialise update of the patch topology
		virtual void initUpdateMesh()
		{}

		//- Update of the patch topology
		virtual void updateMesh()
		{}


public:

	//- Runtime type information
	TypeName("basePatch");


	// Constructor

		pointPatch
		(
			const pointBoundaryMesh& bm
		)
		:
			boundaryMesh_(bm)
		{}


	// Destructor

		virtual ~pointPatch()
		{}


	// Member Functions

		//- Return name
		virtual const word& name() const = 0;

		//- Return size
		virtual label size() const = 0;

		//- Return the index of this patch in the pointBoundaryMesh
		virtual label index() const = 0;

		//- Return boundaryMesh reference
		const pointBoundaryMesh& boundaryMesh() const
		{
			return boundaryMesh_;
		}

		//- Return true if this patch field is coupled
		virtual bool coupled() const
		{
			return false;
		}

		//- Return mesh points
		virtual const labelList& meshPoints() const = 0;

		//- Return mesh points
		virtual const vectorField& localPoints() const = 0;

		//- Return  point normals
		virtual const vectorField& pointNormals() const = 0;

		//- Accumulate the effect of constraint direction of this patch
		virtual void applyConstraint
		(
			const label pointi,
			pointConstraint&
		) const
		{}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
