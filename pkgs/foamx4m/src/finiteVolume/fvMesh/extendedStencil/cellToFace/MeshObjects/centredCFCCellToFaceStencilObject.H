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
	Foam::centredCFCCellToFaceStencilObject

Description

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef centredCFCCellToFaceStencilObject_H
#define centredCFCCellToFaceStencilObject_H

#include "extendedCentredCellToFaceStencil.H"
#include "CFCCellToFaceStencil.H"
#include "MeshObject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class centredCFCCellToFaceStencilObject
:
	public MeshObject<fvMesh, centredCFCCellToFaceStencilObject>,
	public extendedCentredCellToFaceStencil
{

	// Private Member Functions

		//- Disallow default bitwise copy construct
		centredCFCCellToFaceStencilObject
		(
			const centredCFCCellToFaceStencilObject&
		);

		//- Disallow default bitwise assignment
		void operator=(const centredCFCCellToFaceStencilObject&);


public:

	TypeName("centredCFCCellToFaceStencil");

	// Constructors

		//- Construct from uncompacted face stencil
		explicit centredCFCCellToFaceStencilObject
		(
			const fvMesh& mesh
		)
		:
			MeshObject<fvMesh, centredCFCCellToFaceStencilObject>(mesh),
			extendedCentredCellToFaceStencil(CFCCellToFaceStencil(mesh))
		{
			if (extendedCellToFaceStencil::debug)
			{
				Info<< "Generated centred stencil " << type()
					<< nl << endl;
				writeStencilStats(Info, stencil(), map());
			}
		}


	//- Destructor
	virtual ~centredCFCCellToFaceStencilObject()
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
