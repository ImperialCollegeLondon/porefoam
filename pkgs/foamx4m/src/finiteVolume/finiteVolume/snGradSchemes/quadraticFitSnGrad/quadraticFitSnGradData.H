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
	quadraticFitSnGradData

Description
	Data for the quadratic fit correction snGrad scheme

SourceFiles
	quadraticFitSnGradData.C

\*---------------------------------------------------------------------------*/

#ifndef quadraticFitSnGradData_H
#define quadraticFitSnGradData_H

#include "MeshObject.H"
#include "fvMesh.H"
#include "extendedCellToFaceStencil.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class quadraticFitSnGradData
:
	public MeshObject<fvMesh, quadraticFitSnGradData>
{
	// Private data

		const scalar centralWeight_;
		const label dim_;

		//- minimum stencil size
		const label minSize_;

		//- Extended stencil addressing
		extendedCellToFaceStencil stencil_;

		//- For each cell in the mesh store the values which multiply the
		//  values of the stencil to obtain the gradient for each direction
		List<scalarList> fit_;


	// Private member functions

		//- Find the normal direction and i, j and k directions for face faci
		static void findFaceDirs
		(
			vector& idir,        // value changed in return
			vector& jdir,        // value changed in return
			vector& kdir,        // value changed in return
			const fvMesh& mesh,
			const label faci
		);

		label calcFit(const List<point>&, const label faci);


public:

	TypeName("quadraticFitSnGradData");


	// Constructors

		explicit quadraticFitSnGradData
		(
			const fvMesh& mesh,
			const scalar cWeight
		);


	// Destructor

		virtual ~quadraticFitSnGradData()
		{};


	// Member functions

		//- Return reference to the stencil
		const extendedCellToFaceStencil& stencil() const
		{
			return stencil_;
		}

		//- Return reference to fit coefficients
		const List<scalarList>& fit() const { return fit_; }

		//- Delete the data when the mesh moves not implemented
		virtual bool movePoints();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
