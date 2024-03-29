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
	Foam::upwindCPCCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef upwindCPCCellToFaceStencilObject_H
#define upwindCPCCellToFaceStencilObject_H

#include "extendedUpwindCellToFaceStencil.H"
#include "CPCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class upwindCPCCellToFaceStencilObject
:
	public MeshObject<fvMesh, upwindCPCCellToFaceStencilObject>,
	public extendedUpwindCellToFaceStencil
{

public:

	TypeName("upwindCFCCellToFaceStencil");

	// Constructors

		//- Construct from uncompacted face stencil
		explicit upwindCPCCellToFaceStencilObject
		(
			const fvMesh& mesh,
			const bool pureUpwind,
			const scalar minOpposedness
		)
		:
			MeshObject<fvMesh, upwindCPCCellToFaceStencilObject>(mesh),
			extendedUpwindCellToFaceStencil
			(
				CPCCellToFaceStencil(mesh),
				pureUpwind,
				minOpposedness
			)
		{
			if (extendedCellToFaceStencil::debug)
			{
				Info<< "Generated off-centred stencil " << type()
					<< nl << endl;
				writeStencilStats(Info, ownStencil(), ownMap());
			}
		}


	// Destructor

		virtual ~upwindCPCCellToFaceStencilObject()
		{}


	// Member functions

		//- Delete the stencil when the mesh moves
		virtual bool movePoints() const
		{
			deleteObject();
			return true;
		}

		//- Delete the stencil when mesh updates
		virtual bool updateMesh(const mapPolyMesh&) const
		{
			deleteObject();
			return true;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
