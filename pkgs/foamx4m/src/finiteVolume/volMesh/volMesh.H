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
	Foam::volMesh

Description
	Mesh data needed to do the Finite Volume discretisation.

\*---------------------------------------------------------------------------*/

#ifndef volMesh_H
#define volMesh_H

#include "GeoMesh.H"
#include "fvMesh.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class volMesh
:
	public GeoMesh<fvMesh>
{

public:

	// Constructors

		//- Construct from fvMesh
		explicit volMesh(const fvMesh& mesh)
		:
			GeoMesh<fvMesh>(mesh)
		{}


	// Member Functions

		//- Return size
		label size() const
		{
			return size(mesh_);
		}

		//- Return size
		static label size(const Mesh& mesh)
		{
			return mesh.nCells();
		}

		//- Return cell centres
		const volVectorField& C()
		{
			return mesh_.C();
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
