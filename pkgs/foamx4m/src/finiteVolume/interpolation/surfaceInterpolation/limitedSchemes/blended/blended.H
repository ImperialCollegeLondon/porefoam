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
	Foam::blended

Description
	linear/upwind blended differencing scheme.

SourceFiles
	blended.C

\*---------------------------------------------------------------------------*/

#ifndef blended_H
#define blended_H

#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class blended
:
	public limitedSurfaceInterpolationScheme<Type>
{
	// Private data

		const scalar blendingFactor_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		blended(const blended&);

		//- Disallow default bitwise assignment
		void operator=(const blended&);


public:

	//- Runtime type information
	TypeName("blended");


	// Constructors

		//- Construct from mesh, faceFlux and blendingFactor
		blended
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			const scalar blendingFactor
		)
		:
			limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
			blendingFactor_(blendingFactor)
		{}

		//- Construct from mesh and Istream.
		//  The name of the flux field is read from the Istream and looked-up
		//  from the mesh objectRegistry
		blended
		(
			const fvMesh& mesh,
			Istream& is
		)
		:
			limitedSurfaceInterpolationScheme<Type>(mesh, is),
			blendingFactor_(readScalar(is))
		{}

		//- Construct from mesh, faceFlux and Istream
		blended
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& is
		)
		:
			limitedSurfaceInterpolationScheme<Type>(mesh, faceFlux),
			blendingFactor_(readScalar(is))
		{}


	// Member Functions

		//- Return the interpolation limiter
		virtual tmp<surfaceScalarField> limiter
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const
		{
			return tmp<surfaceScalarField>
			(
				new surfaceScalarField
				(
					IOobject
					(
					    "blendedLimiter",
					    this->mesh().time().timeName(),
					    this->mesh()
					),
					this->mesh(),
					dimensionedScalar
					(
					    "blendedLimiter",
					    dimless,
					    1 - blendingFactor_
					)
				)
			);
		}

		//- Return the interpolation weighting factors
		virtual tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			return
				blendingFactor_*this->mesh().surfaceInterpolation::weights()
			  + (1 - blendingFactor_)*pos(this->faceFlux_);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
