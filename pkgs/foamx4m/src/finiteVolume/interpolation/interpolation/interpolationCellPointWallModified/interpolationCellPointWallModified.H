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
	Foam::interpolationCellPoint

Description
	Same as interpolationCellPoint, but if interpolating a wall face, uses
	cell centre value instead

\*---------------------------------------------------------------------------*/

#ifndef interpolationCellPointWallModified_H
#define interpolationCellPointWallModified_H

#include "interpolationCellPoint.H"
#include "cellPointWeightWallModified.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class interpolationCellPointWallModified
:
	public interpolationCellPoint<Type>
{
public:

	//- Runtime type information
	TypeName("cellPointWallModified");


	// Constructors

		//- Construct from components
		interpolationCellPointWallModified
		(
			const GeometricField<Type, fvPatchField, volMesh>& psi
		);


	// Member Functions

		//- Interpolate field for the given cellPointWeight
		inline Type interpolate(const cellPointWeightWallModified& cpw) const;

		//- Interpolate field to the given point in the given cell
		inline Type interpolate
		(
			const vector& position,
			const label nCell,
			const label facei = -1
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "interpolationCellPointWallModifiedI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "interpolationCellPointWallModified.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
