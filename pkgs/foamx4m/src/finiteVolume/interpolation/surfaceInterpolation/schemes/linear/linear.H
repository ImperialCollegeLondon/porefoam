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
	Foam::linear

Description
	Central-differencing interpolation scheme class

SourceFiles
	linear.C

\*---------------------------------------------------------------------------*/

#ifndef linear_H
#define linear_H

#include "surfaceInterpolationScheme.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class linear
:
	public surfaceInterpolationScheme<Type>
{
	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const linear&);


public:

	//- Runtime type information
	TypeName("linear");


	// Constructors

		//- Construct from mesh
		linear(const fvMesh& mesh)
		:
			surfaceInterpolationScheme<Type>(mesh)
		{}

		//- Construct from Istream
		linear(const fvMesh& mesh, Istream&)
		:
			surfaceInterpolationScheme<Type>(mesh)
		{}

		//- Construct from faceFlux and Istream
		linear
		(
			const fvMesh& mesh,
			const surfaceScalarField&,
			Istream&
		)
		:
			surfaceInterpolationScheme<Type>(mesh)
		{}


	// Member Functions

		//- Return the interpolation weighting factors
		tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			return this->mesh().surfaceInterpolation::weights();
		}
};


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
linearInterpolate(const GeometricField<Type, fvPatchField, volMesh>& vf)
{
	return surfaceInterpolationScheme<Type>::interpolate
	(
		vf,
		vf.mesh().surfaceInterpolation::weights()
	);
}


template<class Type>
tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
linearInterpolate(const tmp<GeometricField<Type, fvPatchField, volMesh> >& tvf)
{
	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tinterp =
		linearInterpolate(tvf());
	tvf.clear();
	return tinterp;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
