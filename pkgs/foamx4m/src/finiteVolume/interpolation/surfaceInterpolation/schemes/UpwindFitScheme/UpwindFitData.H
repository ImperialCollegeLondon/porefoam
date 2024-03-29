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
	Foam::UpwindFitData

Description
	Data for the quadratic fit correction interpolation scheme to be used with
	upwind biased stencil.
	- linearCorrection = true : fit calculated for corrected linear scheme
	- linearCorrection = false : fit calculated for corrected upwind scheme

SourceFiles
	UpwindFitData.C

\*---------------------------------------------------------------------------*/

#ifndef UpwindFitData_H
#define UpwindFitData_H

#include "FitData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class extendedUpwindCellToFaceStencil;


template<class Polynomial>
class UpwindFitData
:
	public FitData
	<
		UpwindFitData<Polynomial>,
		extendedUpwindCellToFaceStencil,
		Polynomial
	>
{
	// Private data

		//- For each face of the mesh store the coefficients to multiply the
		//  stencil cell values by if the flow is from the owner
		List<scalarList> owncoeffs_;

		//- For each face of the mesh store the coefficients to multiply the
		//  stencil cell values by if the flow is from the neighbour
		List<scalarList> neicoeffs_;


	// Private member functions

		//- Calculate the fit for the all the mesh faces
		//  and set the coefficients
		void calcFit();


public:

	TypeName("UpwindFitData");


	// Constructors

		//- Construct from components
		UpwindFitData
		(
			const fvMesh& mesh,
			const extendedUpwindCellToFaceStencil& stencil,
			const bool linearCorrection,
			const scalar linearLimitFactor,
			const scalar centralWeight
		);


	//- Destructor
	virtual ~UpwindFitData()
	{}


	// Member functions

		//- Return reference to owner fit coefficients
		const List<scalarList>& owncoeffs() const
		{
			return owncoeffs_;
		}

		//- Return reference to neighbour fit coefficients
		const List<scalarList>& neicoeffs() const
		{
			return neicoeffs_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "UpwindFitData.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
