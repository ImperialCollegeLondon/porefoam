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
	Foam::limitWith

Description
	limitWith differencing scheme limits the specified scheme with the
	specified limiter.

SourceFiles
	limitWith.C

\*---------------------------------------------------------------------------*/

#ifndef limitWith_H
#define limitWith_H

#include "surfaceInterpolationScheme.H"
#include "limitedSurfaceInterpolationScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class limitWith
:
	public surfaceInterpolationScheme<Type>
{
	// Private Member Functions

		//- Interpolation scheme
		tmp<surfaceInterpolationScheme<Type> > tInterp_;

		//- Limiter
		tmp<limitedSurfaceInterpolationScheme<Type> > tLimiter_;


		//- Disallow default bitwise copy construct
		limitWith(const limitWith&);

		//- Disallow default bitwise assignment
		void operator=(const limitWith&);


public:

	//- Runtime type information
	TypeName("limitWith");


	// Constructors

		//- Construct from mesh and Istream.
		//  The name of the flux field is read from the Istream and looked-up
		//  from the mesh objectRegistry
		limitWith
		(
			const fvMesh& mesh,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh),
			tInterp_
			(
				surfaceInterpolationScheme<Type>::New(mesh, is)
			),
			tLimiter_
			(
				limitedSurfaceInterpolationScheme<Type>::New(mesh, is)
			)
		{}

		//- Construct from mesh, faceFlux and Istream
		limitWith
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& is
		)
		:
			surfaceInterpolationScheme<Type>(mesh),
			tInterp_
			(
				surfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
			),
			tLimiter_
			(
				limitedSurfaceInterpolationScheme<Type>::New(mesh, faceFlux, is)
			)
		{}


	// Member Functions

		//- Return the interpolation weighting factors
		virtual tmp<surfaceScalarField> weights
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			return tLimiter_().weights
			(
				vf,
				tInterp_().weights(vf),
				tLimiter_().limiter(vf)
			);
		}

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return tInterp_().corrected();
		}

		//- Return the explicit correction to the face-interpolate
		//  for the given field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction(const GeometricField<Type, fvPatchField, volMesh>& vf) const
		{
			return tLimiter_().limiter(vf)*tInterp_().correction(vf);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
