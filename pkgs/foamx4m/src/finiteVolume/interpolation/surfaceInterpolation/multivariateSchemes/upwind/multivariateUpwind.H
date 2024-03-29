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
	Foam::multivariateUpwind

Description
	Upwind-difference form of the multivariate surfaceInterpolationScheme.

SourceFiles
	multivariateUpwindmake.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateUpwind_H
#define multivariateUpwind_H

#include "multivariateSurfaceInterpolationScheme.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class multivariateUpwind
:
	public multivariateSurfaceInterpolationScheme<Type>
{
	// Private data

		const surfaceScalarField& faceFlux_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		multivariateUpwind(const multivariateUpwind&);

		//- Disallow default bitwise assignment
		void operator=(const multivariateUpwind&);


public:

	//- Runtime type information
	TypeName("upwind");


	// Constructors

		//- Construct for field, faceFlux and Istream
		multivariateUpwind
		(
			const fvMesh& mesh,
			const typename multivariateSurfaceInterpolationScheme<Type>::
				fieldTable& fields,
			const surfaceScalarField& faceFlux,
			Istream& schemeData
		)
		:
			multivariateSurfaceInterpolationScheme<Type>
			(
				mesh,
				fields,
				faceFlux,
				schemeData
			),
			faceFlux_(faceFlux)
		{}


	// Member Operators

		//- surfaceInterpolationScheme sub-class returned by operator(field)
		//  which is used as the interpolation scheme for the field
		class fieldScheme
		:
			public multivariateSurfaceInterpolationScheme<Type>::fieldScheme
		{
			// Private data

				const surfaceScalarField& faceFlux_;

		public:

			// Constructors

				fieldScheme
				(
					const GeometricField<Type, fvPatchField, volMesh>& field,
					const surfaceScalarField& faceFlux
				)
				:
					multivariateSurfaceInterpolationScheme<Type>::
					    fieldScheme(field),
					faceFlux_(faceFlux)
				{}


			// Member Functions

				//- Return the interpolation weighting factors
				tmp<surfaceScalarField> weights
				(
					const GeometricField<Type, fvPatchField, volMesh>&
				) const
				{
					return pos(faceFlux_);
				}
		};

		tmp<surfaceInterpolationScheme<Type> > operator()
		(
			const GeometricField<Type, fvPatchField, volMesh>& field
		) const
		{
			return tmp<surfaceInterpolationScheme<Type> >
			(
				new fieldScheme(field, faceFlux_)
			);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
