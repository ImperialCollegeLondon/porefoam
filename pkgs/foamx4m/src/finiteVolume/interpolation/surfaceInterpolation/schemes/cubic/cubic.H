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
	Foam::cubic

Description
	Cubic interpolation scheme class derived from linear and returns
	linear weighting factors but also applies an explicit correction.

SourceFiles
	cubic.C

\*---------------------------------------------------------------------------*/

#ifndef cubic_H
#define cubic_H

#include "linear.H"
#include "gaussGrad.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class cubic
:
	public linear<Type>
{
	// Private Member Functions

		//- Disallow default bitwise copy construct
		cubic(const cubic&);

		//- Disallow default bitwise assignment
		void operator=(const cubic&);


public:

	//- Runtime type information
	TypeName("cubic");


	// Constructors

		//- Construct from mesh
		cubic(const fvMesh& mesh)
		:
			linear<Type>(mesh)
		{}

		//- Construct from mesh and Istream
		cubic
		(
			const fvMesh& mesh,
			Istream&
		)
		:
			linear<Type>(mesh)
		{}

		//- Construct from mesh, faceFlux and Istream
		cubic
		(
			const fvMesh& mesh,
			const surfaceScalarField&,
			Istream&
		)
		:
			linear<Type>(mesh)
		{}


	// Member Functions

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return true;
		}

		//- Return the explicit correction to the face-interpolate
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction
		(
			const GeometricField<Type, fvPatchField, volMesh>& vf
		) const
		{
			const fvMesh& mesh = this->mesh();

			// calculate the appropriate interpolation factors
			const surfaceScalarField& lambda = mesh.weights();

			surfaceScalarField kSc =
				lambda*(scalar(1) - lambda*(scalar(3) - scalar(2)*lambda));

			surfaceScalarField kVecP = sqr(scalar(1) - lambda)*lambda;
			surfaceScalarField kVecN = sqr(lambda)*(lambda - scalar(1));

			tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
			(
				new GeometricField<Type, fvsPatchField, surfaceMesh>
				(
					IOobject
					(
					    vf.name(),
					    mesh.time().timeName(),
					    mesh
					),
					surfaceInterpolationScheme<Type>::interpolate(vf, kSc, -kSc)
				)
			);

			GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

			for (direction cmpt=0; cmpt<pTraits<Type>::nComponents; cmpt++)
			{
				sfCorr.replace
				(
					cmpt,
					sfCorr.component(cmpt)
				  + (
					    surfaceInterpolationScheme
					    <
					        typename outerProduct
					        <
					            vector,
					            typename pTraits<Type>::cmptType
					        >::type
					    >::interpolate
					    (
					        fv::gaussGrad
					        <typename pTraits<Type>::cmptType>(mesh)
					       .grad(vf.component(cmpt)),
					        kVecP,
					        kVecN
					    ) & mesh.Sf()
					)/mesh.magSf()/mesh.surfaceInterpolation::deltaCoeffs()
				);
			}

			forAll (sfCorr.boundaryField(), pi)
			{
				if (!sfCorr.boundaryField()[pi].coupled())
				{
					sfCorr.boundaryField()[pi] = pTraits<Type>::zero;
				}
			}

			return tsfCorr;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
