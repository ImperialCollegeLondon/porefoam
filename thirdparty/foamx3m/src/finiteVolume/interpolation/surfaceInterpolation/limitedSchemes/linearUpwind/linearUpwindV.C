/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
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

	GeometricField
		<typename outerProduct<vector, Type>::type, fvPatchField, volMesh>
		gradVf = gradScheme_().grad(vf);

	// Note: in order for the patchNeighbourField to be correct on coupled
	// boundaries, correctBoundaryConditions needs to be called.
	// The call shall be moved into the library fvc operators
	gradVf.correctBoundaryConditions();

	forAll(faceFlux, faceI)
	{
		vector maxCorr;

		if (faceFlux[faceI] > 0)
		{
			maxCorr =
				(1 - w[faceI])*(vf[nei[faceI]] - vf[own[faceI]]);

			sfCorr[faceI] =
				(Cf[faceI] - C[own[faceI]]) & gradVf[own[faceI]];
		}
		else
		{
			maxCorr =
				w[faceI]*(vf[own[faceI]] - vf[nei[faceI]]);

			sfCorr[faceI] =
				(Cf[faceI] - C[nei[faceI]]) & gradVf[nei[faceI]];
		}

		scalar sfCorrs = magSqr(sfCorr[faceI]);
		scalar maxCorrs = sfCorr[faceI] & maxCorr;

		if (sfCorrs > 0)
		{
			if (maxCorrs < 0)
			{
				sfCorr[faceI] = vector::zero;
			}
			else if (sfCorrs > maxCorrs)
			{
				sfCorr[faceI] *= maxCorrs/(sfCorrs + VSMALL);
			}
		}
		else if (sfCorrs < 0)
		{
			if (maxCorrs > 0)
			{
				sfCorr[faceI] = vector::zero;
			}
			else if (sfCorrs < maxCorrs)
			{
				sfCorr[faceI] *= maxCorrs/(sfCorrs - VSMALL);
			}
		}
	}

	// Added missing treatment of coupled boundaries.  HJ, 27/Jul/2011

	typename GeometricField<Type, fvsPatchField, surfaceMesh>::Boundary& bSfCorr = sfCorr.boundaryField();

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

			forAll(pOwner, faceI)
			{
				vector maxCorr;

				label own = pOwner[faceI];

				if (pFaceFlux[faceI] > 0)
				{
					maxCorr =
						(1.0 - pW[faceI])*(vfNei[faceI] - vfOwn[faceI]);

					pSfCorr[faceI] = (pCf[faceI] - C[own]) & gradVf[own];
				}
				else
				{
					maxCorr =
						pW[faceI]*(vfOwn[faceI] - vfNei[faceI]);

					pSfCorr[faceI] =
						(pCf[faceI] - pd[faceI] - C[own]) & pGradVfNei[faceI];
				}

				scalar pSfCorrs = magSqr(pSfCorr[faceI]);
				scalar maxCorrs = pSfCorr[faceI] & maxCorr;

				if (pSfCorrs > 0)
				{
					if (maxCorrs < 0)
					{
						pSfCorr[faceI] = vector::zero;
					}
					else if (pSfCorrs > maxCorrs)
					{
						pSfCorr[faceI] *= maxCorrs/(pSfCorrs + VSMALL);
					}
				}
				else if (pSfCorrs < 0)
				{
					if (maxCorrs > 0)
					{
						pSfCorr[faceI] = vector::zero;
					}
					else if (pSfCorrs < maxCorrs)
					{
						pSfCorr[faceI] *= maxCorrs/(pSfCorrs - VSMALL);
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
