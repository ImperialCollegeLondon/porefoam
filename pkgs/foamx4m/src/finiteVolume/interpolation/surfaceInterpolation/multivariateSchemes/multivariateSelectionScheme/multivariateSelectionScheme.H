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
	Foam::multivariateSelectionScheme

Description
	Generic multi-variate discretisation scheme class for which any of the
	NVD, CNVD or NVDV schemes may be selected for each variable.

SourceFiles
	multivariateSelectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateSelectionScheme_H
#define multivariateSelectionScheme_H

#include "multivariateSurfaceInterpolationScheme.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class multivariateSelectionScheme
:
	public multivariateSurfaceInterpolationScheme<Type>
{
	// Private data

		//- Schemes dictionary
		dictionary schemes_;

		//- Face flux
		const surfaceScalarField& faceFlux_;

		//- Combined weights
		surfaceScalarField weights_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		multivariateSelectionScheme(const multivariateSelectionScheme&);

		//- Disallow default bitwise assignment
		void operator=(const multivariateSelectionScheme&);


public:

	//- Runtime type information
	TypeName("multivariateSelection");


	// Constructors

		//- Construct for field, faceFlux and Istream
		multivariateSelectionScheme
		(
			const fvMesh& mesh,
			const typename multivariateSurfaceInterpolationScheme<Type>::
				fieldTable& fields,
			const surfaceScalarField& faceFlux,
			Istream& schemeData
		);


	// Member Operators

		//- surfaceInterpolationScheme sub-class returned by operator(field)
		//  which is used as the interpolation scheme for the field
		class fieldScheme
		:
			public multivariateSurfaceInterpolationScheme<Type>::fieldScheme
		{
			// Private data

				const surfaceScalarField& weights_;

		public:

			// Constructors

				fieldScheme
				(
					const GeometricField<Type, fvPatchField, volMesh>& field,
					const surfaceScalarField& weights
				)
				:
					multivariateSurfaceInterpolationScheme<Type>::
					    fieldScheme(field),
					weights_(weights)
				{}


			// Member Functions

				//- Return the interpolation weighting factors
				tmp<surfaceScalarField> weights
				(
					const GeometricField<Type, fvPatchField, volMesh>&
				) const
				{
					return weights_;
				}
		};

		tmp<surfaceInterpolationScheme<Type> > operator()
		(
			const GeometricField<Type, fvPatchField, volMesh>& field
		) const
		{
			return tmp<surfaceInterpolationScheme<Type> >
			(
				new fieldScheme(field, weights_)
			);
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "multivariateSelectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
