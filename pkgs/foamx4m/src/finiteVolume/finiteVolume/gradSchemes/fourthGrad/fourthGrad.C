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

#include "fourthGrad.H"
#include "leastSquaresGrad.H"
#include "gaussGrad.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "surfaceMesh.H"
#include "GeometricField.H"
#include "zeroGradientFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class Type>
tmp
<
	GeometricField
	<
		typename outerProduct<vector, Type>::type, fvPatchField, volMesh
	>
>
fourthGrad<Type>::calcGrad
(
	const GeometricField<Type, fvPatchField, volMesh>& vsf,
	const word& name
) const
{
	// The fourth-order gradient is calculated in two passes.  First,
	// the standard least-square gradient is assembled.  Then, the
	// fourth-order correction is added to the second-order accurate
	// gradient to complete the accuracy.

	typedef typename outerProduct<vector, Type>::type GradType;

	const fvMesh& mesh = vsf.mesh();

	// Assemble the second-order least-square gradient
	// Calculate the second-order least-square gradient
	tmp<GeometricField<GradType, fvPatchField, volMesh> > tsecondfGrad
		= leastSquaresGrad<Type>(mesh).grad(vsf);
	const GeometricField<GradType, fvPatchField, volMesh>& secondfGrad =
		tsecondfGrad();

	tmp<GeometricField<GradType, fvPatchField, volMesh> > tfGrad
	(
		new GeometricField<GradType, fvPatchField, volMesh>
		(
			IOobject
			(
				name,
				vsf.instance(),
				mesh,
				IOobject::NO_READ,
				IOobject::NO_WRITE
			),
			secondfGrad
		)
	);
	GeometricField<GradType, fvPatchField, volMesh>& fGrad = tfGrad();

	const vectorField& C = mesh.C();

	const surfaceScalarField& lambda = mesh.weights();

	// Get reference to least square vectors
	const leastSquaresVectors& lsv = leastSquaresVectors::New(mesh);
	const surfaceVectorField& ownLs = lsv.pVectors();
	const surfaceVectorField& neiLs = lsv.nVectors();

	// owner/neighbour addressing
	const unallocLabelList& own = mesh.owner();
	const unallocLabelList& nei = mesh.neighbour();

	// Assemble the fourth-order gradient

	// Internal faces
	forAll(own, facei)
	{
		Type dDotGradDelta = 0.5*
		(
		   (C[nei[facei]] - C[own[facei]])
		 & (secondfGrad[nei[facei]] - secondfGrad[own[facei]])
		);

		fGrad[own[facei]] -= lambda[facei]*ownLs[facei]*dDotGradDelta;
		fGrad[nei[facei]] -= (1.0 - lambda[facei])*neiLs[facei]*dDotGradDelta;
	}

	// Boundary faces
	forAll(vsf.boundaryField(), patchi)
	{
		if (secondfGrad.boundaryField()[patchi].coupled())
		{
			const fvPatch& p = fGrad.boundaryField()[patchi].patch();

			const fvsPatchVectorField& patchOwnLs =
				ownLs.boundaryField()[patchi];

			const scalarField& lambdap = lambda.boundaryField()[patchi];

			// Build the d-vectors
			// Better version of d-vectors: Zeljko Tukovic, 25/Apr/2010
			vectorField pd = p.delta();

			const unallocLabelList& faceCells = p.faceCells();

			Field<GradType> neighbourSecondfGrad =
				secondfGrad.boundaryField()[patchi].patchNeighbourField();

			forAll(faceCells, patchFaceI)
			{
				fGrad[faceCells[patchFaceI]] -=
					0.5*lambdap[patchFaceI]*patchOwnLs[patchFaceI]
				   *(
					    pd[patchFaceI]
					  & (
					        neighbourSecondfGrad[patchFaceI]
					      - secondfGrad[faceCells[patchFaceI]]
					    )
					);
			}
		}
	}

	fGrad.correctBoundaryConditions();
	gaussGrad<Type>::correctBoundaryConditions(vsf, fGrad);

	return tfGrad;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
