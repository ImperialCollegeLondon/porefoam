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

#include "linearUpwindV.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
Foam::tmp<Foam::GeometricField<Type, Foam::fvsPatchField, Foam::surfaceMesh> >
Foam::linearUpwindV<Type>::correction
(
	const GeometricField<Type, fvPatchField, volMesh>& vf
) const
{
	const fvMesh& mesh = this->mesh();

	tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tsfCorr
	(
		new GeometricField<Type, fvsPatchField, surfaceMesh>
		(
			IOobject
			(
				"linearUpwindCorrection(" + vf.name() + ')',
				mesh.time().timeName(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE,
				false
			),
			mesh,
			dimensioned<Type>(vf.name(), vf.dimensions(), pTraits<Type>::zero)
		)
	);

	GeometricField<Type, fvsPatchField, surfaceMesh>& sfCorr = tsfCorr();

	const surfaceScalarField& faceFlux = this->faceFlux_;
	const surfaceScalarField& w = mesh.weights();

	const labelList& own = mesh.owner();
	const labelList& nei = mesh.neighbour();

	const volVectorField& C = mesh.C();
	const surfaceVectorField& Cf = mesh.Cf();

	// Due to gradient cacheing, must take a tmp field
	// HJ, 22/Apr/2016
	tmp
	<
		GeometricField
		<
			typename outerProduct<vector, Type>::type,
			fvPatchField,
			volMesh
		>
	> tgradVf = gradScheme_().grad(vf, gradSchemeName_);

	const GeometricField
	<
		typename outerProduct<vector, Type>::type,
		fvPatchField,
		volMesh
	>& gradVf = tgradVf();

	// Note: in order for the patchNeighbourField to be correct on coupled
	// boundaries, correctBoundaryConditions needs to be called.
	// The call shall be moved into the library fvc operators
	// Moved to cached gradScheme.  HJ, 22/Apr/2016
//     gradVf.correctBoundaryConditions();

	forAll(faceFlux, facei)
	{
		vector maxCorr;

		if (faceFlux[facei] > 0)
		{
			maxCorr =
				(1 - w[facei])*(vf[nei[facei]] - vf[own[facei]]);

			sfCorr[facei] =
				(Cf[facei] - C[own[facei]]) & gradVf[own[facei]];
		}
		else
		{
			maxCorr =
				w[facei]*(vf[own[facei]] - vf[nei[facei]]);

			sfCorr[facei] =
				(Cf[facei] - C[nei[facei]]) & gradVf[nei[facei]];
		}

		scalar sfCorrs = magSqr(sfCorr[facei]);
		scalar maxCorrs = sfCorr[facei] & maxCorr;

		if (sfCorrs > 0)
		{
			if (maxCorrs < 0)
			{
				sfCorr[facei] = vector::zero;
			}
			else if (sfCorrs > maxCorrs)
			{
				sfCorr[facei] *= maxCorrs/(sfCorrs + VSMALL);
			}
		}
		else if (sfCorrs < 0)
		{
			if (maxCorrs > 0)
			{
				sfCorr[facei] = vector::zero;
			}
			else if (sfCorrs < maxCorrs)
			{
				sfCorr[facei] *= maxCorrs/(sfCorrs - VSMALL);
			}
		}
	}

	// Added missing treatment of coupled boundaries.  HJ, 27/Jul/2011

	typename GeometricField<Type, fvsPatchField, surfaceMesh>::
		Boundary& bSfCorr = sfCorr.boundaryField();

	forAll(bSfCorr, patchi)
	{
		fvsPatchField<Type>& pSfCorr = bSfCorr[patchi];

		if (pSfCorr.coupled())
		{
			const fvPatch& p = mesh.boundary()[patchi];

			const unallocLabelList& pOwner = p.faceCells();

			const vectorField& pCf = Cf.boundaryField()[patchi];

			const scalarField& pFaceFlux = faceFlux.boundaryField()[patchi];

			const scalarField& pW = w.boundaryField()[patchi];

			Field<Type> vfOwn =
				vf.boundaryField()[patchi].patchInternalField();

			Field<Type> vfNei =
				vf.boundaryField()[patchi].patchNeighbourField();


			Field<typename outerProduct<vector, Type>::type> pGradVfNei =
				gradVf.boundaryField()[patchi].patchNeighbourField();

			// Build the d-vectors
			// Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
			vectorField pd = p.delta();

			forAll(pOwner, facei)
			{
				vector maxCorr;

				label own = pOwner[facei];

				if (pFaceFlux[facei] > 0)
				{
					maxCorr =
					    (1.0 - pW[facei])*(vfNei[facei] - vfOwn[facei]);

					pSfCorr[facei] = (pCf[facei] - C[own]) & gradVf[own];
				}
				else
				{
					maxCorr =
					    pW[facei]*(vfOwn[facei] - vfNei[facei]);

					pSfCorr[facei] =
					    (pCf[facei] - pd[facei] - C[own]) & pGradVfNei[facei];
				}

				scalar pSfCorrs = magSqr(pSfCorr[facei]);
				scalar maxCorrs = pSfCorr[facei] & maxCorr;

				if (pSfCorrs > 0)
				{
					if (maxCorrs < 0)
					{
					    pSfCorr[facei] = vector::zero;
					}
					else if (pSfCorrs > maxCorrs)
					{
					    pSfCorr[facei] *= maxCorrs/(pSfCorrs + VSMALL);
					}
				}
				else if (pSfCorrs < 0)
				{
					if (maxCorrs > 0)
					{
					    pSfCorr[facei] = vector::zero;
					}
					else if (pSfCorrs < maxCorrs)
					{
					    pSfCorr[facei] *= maxCorrs/(pSfCorrs - VSMALL);
					}
				}
			}
		}
	}

	return tsfCorr;
}


namespace Foam
{
	makelimitedSurfaceInterpolationTypeScheme(linearUpwindV, vector)
}

// ************************************************************************* //
