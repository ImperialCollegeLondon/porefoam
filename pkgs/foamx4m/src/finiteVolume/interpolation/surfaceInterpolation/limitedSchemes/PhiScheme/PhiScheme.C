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

\*---------------------------------------------------------------------------*/

#include "volFields.H"
#include "surfaceFields.H"
#include "fvcGrad.H"
#include "coupledFvPatchFields.H"
#include "surfaceInterpolate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type, class PhiLimiter>
tmp<surfaceScalarField> PhiScheme<Type, PhiLimiter>::limiter
(
	const GeometricField<Type, fvPatchField, volMesh>& phi
) const
{
	const fvMesh& mesh = this->mesh();

	tmp<surfaceScalarField> tLimiter
	(
		new surfaceScalarField
		(
			IOobject
			(
				"PhiLimiter",
				mesh.time().timeName(),
				mesh
			),
			mesh,
			dimless
		)
	);
	surfaceScalarField& Limiter = tLimiter();

	const surfaceScalarField& CDweights = mesh.surfaceInterpolation::weights();

	const surfaceVectorField& Sf = mesh.Sf();
	const surfaceScalarField& magSf = mesh.magSf();

	const unallocLabelList& owner = mesh.owner();
	const unallocLabelList& neighbour = mesh.neighbour();

	tmp<surfaceScalarField> tUflux = this->faceFlux_;

	if (this->faceFlux_.dimensions() == dimDensity*dimVelocity*dimArea)
	{
		const volScalarField& rho =
			phi.db().objectRegistry::template
			lookupObject<volScalarField>("rho");

		tUflux = this->faceFlux_/fvc::interpolate(rho);
	}
	else if (this->faceFlux_.dimensions() != dimVelocity*dimArea)
	{
		FatalErrorIn
		(
			"PhiScheme<PhiLimiter>::limiter"
			"(const GeometricField<Type, fvPatchField, volMesh>& phi)"
		)   << "dimensions of faceFlux are not correct"
			<< exit(FatalError);
	}

	const surfaceScalarField& Uflux = tUflux();

	scalarField& pLimiter = Limiter.internalField();

	forAll(pLimiter, face)
	{
		pLimiter[face] = PhiLimiter::limiter
		(
			CDweights[face],
			Uflux[face],
			phi[owner[face]],
			phi[neighbour[face]],
			Sf[face],
			magSf[face]
		);
	}


	surfaceScalarField::Boundary& bLimiter =
		Limiter.boundaryField();

	forAll(bLimiter, patchI)
	{
		scalarField& pLimiter = bLimiter[patchI];

		if (bLimiter[patchI].coupled())
		{
			const scalarField& pCDweights = CDweights.boundaryField()[patchI];
			const vectorField& pSf = Sf.boundaryField()[patchI];
			const scalarField& pmagSf = magSf.boundaryField()[patchI];
			const scalarField& pFaceFlux = Uflux.boundaryField()[patchI];
			Field<Type> pphiP =
				phi.boundaryField()[patchI].patchInternalField();
			Field<Type> pphiN =
				phi.boundaryField()[patchI].patchNeighbourField();

			forAll(pLimiter, face)
			{
				pLimiter[face] = PhiLimiter::limiter
				(
					pCDweights[face],
					pFaceFlux[face],
					pphiP[face],
					pphiN[face],
					pSf[face],
					pmagSf[face]
				);
			}
		}
		else
		{
			pLimiter = 1.0;
		}
	}

	return tLimiter;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
