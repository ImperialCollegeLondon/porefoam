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
	Foam::fvBoundaryMesh

Description
	Foam::fvBoundaryMesh

SourceFiles
	fvBoundaryMesh.C

\*---------------------------------------------------------------------------*/

#ifndef fvBoundaryMesh_H
#define fvBoundaryMesh_H

#include "fvPatchList.H"
#include "lduInterfacePtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;
class polyBoundaryMesh;


class fvBoundaryMesh
:
	public fvPatchList
{
private:

	// Private data

		//- Reference to mesh
		const fvMesh& mesh_;


	// Private Member Functions

		//- Disable default copy construct
		fvBoundaryMesh(const fvBoundaryMesh&);

		//- Disallow assignment
		void operator=(const fvBoundaryMesh&);

		//- Add fvPatches corresponding to the polyBoundaryMesh
		void addFvPatches();

		//- Reset fvPatches corresponding to the polyBoundaryMesh
		//  When the mesh is resized, the existing patches may be preserved
		//  or reset depending on the resetFlag
		void resetFvPatches(const boolList& resetFvPatchFlag);


protected:

		//- Update boundary based on new polyBoundaryMesh
		void readUpdate();


public:

	friend class fvMesh;


	// Constructors

		//- Construct from polyBoundaryMesh
		explicit fvBoundaryMesh(const fvMesh&);


	// Destructor - default


	// Member functions

		// Access

			//- Return the mesh reference
			const fvMesh& mesh() const
			{
				return mesh_;
			}

			//- Return a list of pointers for each patch
			//  with only those pointing to interfaces being set
			lduInterfacePtrsList interfaces() const;


		//- Correct patches after moving points
		void movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
