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

Description
	Specialisation of leastSquaresDivScheme for vectors.
	Needed for implicit fvmDiv operator for block coupled systems.

\*---------------------------------------------------------------------------*/

#include "vectorLeastSquaresDivScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<BlockLduSystem<vector, scalar> > leastSquaresDivScheme<vector>::fvmUDiv
(
	const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
	// Get weights
	// Note: weights are only used for boundary contributions in coefficients
	// HJ, 30/Jul/2015
	surfaceScalarField weights = this->tinterpScheme_().weights(vf);

	const fvMesh& mesh = vf.mesh();

	const unallocLabelList& own = mesh.owner();
	const unallocLabelList& nei = mesh.neighbour();

	volScalarField cellV
	(
		IOobject
		(
			"cellV",
			vf.instance(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimVolume, 0),
		zeroGradientFvPatchScalarField::typeName
	);
	cellV.internalField() = mesh.V();
	cellV.correctBoundaryConditions();
	const scalarField& cellVIn = cellV.internalField();

	tmp<BlockLduSystem<vector, scalar> > tbs
	(
		new BlockLduSystem<vector, scalar>(mesh)
	);
	BlockLduSystem<vector, scalar>& bs = tbs();
	scalarField& source = bs.source();

	// Grab ldu parts of block matrix as linear always
	CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
	CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
	CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

	// Get references to least square vectors
	const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);

	const surfaceVectorField& ownLs = lsv.pVectors();
	const surfaceVectorField& neiLs = lsv.nVectors();

	// Internal field
	const vectorField& ownLsIn = ownLs.internalField();
	const vectorField& neiLsIn = neiLs.internalField();
//     const scalarField& wIn = weights.internalField();

	label owner, neighbour;

	forAll (nei, faceI)
	{
		owner = own[faceI];
		neighbour = nei[faceI];

		u[faceI] = cellVIn[owner]*ownLsIn[faceI];
		l[faceI] = cellVIn[neighbour]*neiLsIn[faceI];

		// Caution - this is NOT negSumDiag(). VV, 17/July/2014.
		d[owner] -= u[faceI];
		d[neighbour] -= l[faceI];
	}

	// Boundary contributions
	forAll(vf.boundaryField(), patchI)
	{
		const fvPatchVectorField& pf = vf.boundaryField()[patchI];
		const fvPatch& patch = pf.patch();
		const vectorField& pownLs = ownLs.boundaryField()[patchI];
		const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
		const unallocLabelList& fc = patch.faceCells();

		// Part of diagonal contribution irrespective of the patch type
		forAll (pf, faceI)
		{
			const label cellI = fc[faceI];
			d[cellI] -= cellVIn[cellI]*pownLs[faceI];
		}

		if (patch.coupled())
		{
			const vectorField& pneiLs = neiLs.boundaryField()[patchI];
			const scalarField cellVInNei =
				cellV.boundaryField()[patchI].patchNeighbourField();

			CoeffField<vector>::linearTypeField& pcoupleUpper =
				bs.coupleUpper()[patchI].asLinear();
			CoeffField<vector>::linearTypeField& pcoupleLower =
				bs.coupleLower()[patchI].asLinear();

			// Coupling  and diagonal contributions
			forAll (pf, faceI)
			{
				pcoupleUpper[faceI] -= cellVIn[fc[faceI]]*pownLs[faceI];
				pcoupleLower[faceI] -= cellVInNei[faceI]*pneiLs[faceI];
			}
		}
		else
		{
			const vectorField internalCoeffs(pf.valueInternalCoeffs(pw));
			const vectorField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

			// Diagonal and source contributions depending on the patch type
			forAll (pf, faceI)
			{
				const label cellI = fc[faceI];

				d[cellI] += cellVIn[cellI]*
					cmptMultiply(pownLs[faceI], internalCoeffs[faceI]);

				source[cellI] -= cellVIn[cellI]*
					(pownLs[faceI] & boundaryCoeffs[faceI]);
			}
		}
	}

	return tbs;
}


template<>
tmp<BlockLduSystem<vector, scalar> > leastSquaresDivScheme<vector>::fvmUDiv
(
	const surfaceScalarField& flux,
	const GeometricField<vector, fvPatchField, volMesh>& vf
) const
{
	// Get weights
	// Note: weights are only used for boundary contributions in coefficients
	// HJ, 30/Jul/2015
	surfaceScalarField weights = this->tinterpScheme_().weights(vf);

	const fvMesh& mesh = vf.mesh();

	const unallocLabelList& own = mesh.owner();
	const unallocLabelList& nei = mesh.neighbour();

	volScalarField cellV
	(
		IOobject
		(
			"cellV",
			vf.instance(),
			mesh,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		mesh,
		dimensionedScalar("zero", dimVolume, 0),
		zeroGradientFvPatchScalarField::typeName
	);
	cellV.internalField() = mesh.V();
	cellV.correctBoundaryConditions();
	const scalarField& cellVIn = cellV.internalField();

	tmp<BlockLduSystem<vector, scalar> > tbs
	(
		new BlockLduSystem<vector, scalar>(mesh)
	);
	BlockLduSystem<vector, scalar>& bs = tbs();
	scalarField& source = bs.source();

	// Grab ldu parts of block matrix as linear always
	CoeffField<vector>::linearTypeField& d = bs.diag().asLinear();
	CoeffField<vector>::linearTypeField& u = bs.upper().asLinear();
	CoeffField<vector>::linearTypeField& l = bs.lower().asLinear();

	// Get references to least square vectors
	const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);

	const surfaceVectorField& ownLs = lsv.pVectors();
	const surfaceVectorField& neiLs = lsv.nVectors();

	// Internal field
	const vectorField& ownLsIn = ownLs.internalField();
	const vectorField& neiLsIn = neiLs.internalField();
//     const scalarField& wIn = weights.internalField();
	const scalarField& fluxIn = flux.internalField();

	forAll (nei, faceI)
	{
		label owner = own[faceI];
		label neighbour = nei[faceI];

		u[faceI] = cellVIn[owner]*ownLsIn[faceI]*fluxIn[faceI];
		l[faceI] = cellVIn[neighbour]*neiLsIn[faceI]*fluxIn[faceI];

		// Caution - this is NOT negSumDiag(). VV, 17/July/2014.
		d[owner] -= u[faceI];
		d[neighbour] -= l[faceI];
	}

	// Boundary contributions
	forAll(vf.boundaryField(), patchI)
	{
		const fvPatchVectorField& pf = vf.boundaryField()[patchI];
		const fvPatch& patch = pf.patch();
		const vectorField& pownLs = ownLs.boundaryField()[patchI];
		const fvsPatchScalarField& pw = weights.boundaryField()[patchI];
		const unallocLabelList& fc = patch.faceCells();

		const scalarField& pFlux = flux.boundaryField()[patchI];

		// Part of diagonal contribution irrespective of the patch type
		forAll (pf, faceI)
		{
			const label cellI = fc[faceI];
			d[cellI] -= cellVIn[cellI]*pownLs[faceI]*pFlux[faceI];
		}

		if (patch.coupled())
		{
			const vectorField& pneiLs = neiLs.boundaryField()[patchI];
			const scalarField cellVInNei =
				cellV.boundaryField()[patchI].patchNeighbourField();

			CoeffField<vector>::linearTypeField& pcoupleUpper =
				bs.coupleUpper()[patchI].asLinear();
			CoeffField<vector>::linearTypeField& pcoupleLower =
				bs.coupleLower()[patchI].asLinear();

			// Coupling  and diagonal contributions
			forAll (pf, faceI)
			{
				pcoupleUpper[faceI] -=
					cellVIn[fc[faceI]]*pownLs[faceI]*pFlux[faceI];
				pcoupleLower[faceI] -=
					cellVInNei[faceI]*pneiLs[faceI]*pFlux[faceI];
			}
		}
		else
		{
			const vectorField internalCoeffs(pf.valueInternalCoeffs(pw));
			const vectorField boundaryCoeffs(pf.valueBoundaryCoeffs(pw));

			// Diagonal and source contributions depending on the patch type
			forAll (pf, faceI)
			{
				const label cellI = fc[faceI];
				d[cellI] += cellVIn[cellI]*pFlux[faceI]*
					cmptMultiply(pownLs[faceI], internalCoeffs[faceI]);
				source[cellI] -= cellVIn[cellI]*pFlux[faceI]*
					(pownLs[faceI] & boundaryCoeffs[faceI]);
			}
		}
	}

	return tbs;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
