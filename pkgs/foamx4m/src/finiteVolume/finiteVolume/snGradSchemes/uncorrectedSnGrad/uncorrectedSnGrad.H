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
	Foam::fv::uncorrectedSnGrad

Description
	Simple central-difference snGrad scheme without non-orthogonal correction.

SourceFiles
	uncorrectedSnGrad.C

\*---------------------------------------------------------------------------*/

#ifndef uncorrectedSnGrad_H
#define uncorrectedSnGrad_H

#include "snGradScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class uncorrectedSnGrad
:
	public snGradScheme<Type>
{
	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const uncorrectedSnGrad&);


public:

	//- Runtime type information
	TypeName("uncorrected");


	// Constructors

		//- Construct from mesh
		uncorrectedSnGrad(const fvMesh& mesh)
		:
			snGradScheme<Type>(mesh)
		{}


		//- Construct from mesh and data stream
		uncorrectedSnGrad(const fvMesh& mesh, Istream&)
		:
			snGradScheme<Type>(mesh)
		{}


	// Destructor

		virtual ~uncorrectedSnGrad();


	// Member Functions

		//- Return the interpolation weighting factors for the given field
		virtual tmp<surfaceScalarField> deltaCoeffs
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			return this->mesh().deltaCoeffs();
		}

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return false;
		}

		//- Return the explicit correction to the uncorrectedSnGrad
		//  for the given field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction(const GeometricField<Type, fvPatchField, volMesh>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "uncorrectedSnGrad.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
