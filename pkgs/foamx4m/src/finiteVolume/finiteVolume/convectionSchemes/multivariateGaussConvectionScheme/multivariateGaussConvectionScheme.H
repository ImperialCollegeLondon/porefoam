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
	Foam::fv::multivariateGaussConvectionScheme

Description
	Basic second-order convection using face-gradients and Gauss' theorem.

SourceFiles
	multivariateGaussConvectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef multivariateGaussConvectionScheme_H
#define multivariateGaussConvectionScheme_H

#include "convectionScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class multivariateGaussConvectionScheme
:
	public fv::convectionScheme<Type>
{
	// Private data

		tmp<multivariateSurfaceInterpolationScheme<Type> > tinterpScheme_;


public:

	//- Runtime type information
	TypeName("Gauss");


	// Constructors

		//- Construct from flux and Istream
		multivariateGaussConvectionScheme
		(
			const fvMesh& mesh,
			const typename multivariateSurfaceInterpolationScheme<Type>::
				fieldTable& fields,
			const surfaceScalarField& faceFlux,
			Istream& is
		)
		:
			convectionScheme<Type>(mesh, faceFlux),
			tinterpScheme_
			(
				multivariateSurfaceInterpolationScheme<Type>::New
				(
					mesh, fields, faceFlux, is
				)
			)
		{}


	// Member Functions

		tmp<multivariateSurfaceInterpolationScheme<Type> >
		interpolationScheme() const
		{
			return tinterpScheme_;
		}

		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > interpolate
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

		tmp<fvMatrix<Type> > fvmDiv
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "multivariateGaussConvectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
