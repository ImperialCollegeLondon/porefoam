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

#include "GaussSeidelSmoother.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(GaussSeidelSmoother, 0);

	lduSmoother::addsymMatrixConstructorToTable<GaussSeidelSmoother>
		addGaussSeidelSmootherSymMatrixConstructorToTable_;

	lduSmoother::addasymMatrixConstructorToTable<GaussSeidelSmoother>
		addGaussSeidelSmootherAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::GaussSeidelSmoother::GaussSeidelSmoother
(
	const lduMatrix& matrix,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const FieldField<Field, scalar>& coupleIntCoeffs,
	const lduInterfaceFieldPtrsList& interfaces
)
:
	lduSmoother
	(
		matrix,
		coupleBouCoeffs,
		coupleIntCoeffs,
		interfaces
	)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::GaussSeidelSmoother::smooth
(
	scalarField& x,
	const lduMatrix& matrix,
	const scalarField& b,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const lduInterfaceFieldPtrsList& interfaces,
	const direction cmpt,
	const label nSweeps
)
{
	scalar* __restrict__ xPtr = x.begin();

	const label nCells = x.size();

	scalarField bPrime(nCells);
	scalar* __restrict__ bPrimePtr = bPrime.begin();

	const scalar* const __restrict__ diagPtr = matrix.diag().begin();
	const scalar* const __restrict__ upperPtr =
		matrix.upper().begin();
	const scalar* const __restrict__ lowerPtr =
		matrix.lower().begin();

	const label* const __restrict__ uPtr =
		matrix.lduAddr().upperAddr().begin();

	const label* const __restrict__ ownStartPtr =
		matrix.lduAddr().ownerStartAddr().begin();


	// Coupled boundary initialisation.  The coupled boundary is treated
	// as an effective jacobi interface in the boundary.
	// Note: there is a change of sign in the coupled
	// interface update.  The reason for this is that the
	// internal coefficients are all located at the l.h.s. of
	// the matrix whereas the "implicit" coefficients on the
	// coupled boundaries are all created as if the
	// coefficient contribution is of a b-kind (i.e. they
	// have a sign as if they are on the r.h.s. of the matrix.
	// To compensate for this, it is necessary to turn the
	// sign of the contribution.

	// Handled by LHS switch on initMatrixInterfaces and updateMatrixInterfaces
	// HJ, 22/May/2013

	for (label sweep = 0; sweep < nSweeps; sweep++)
	{
		bPrime = b;

		// Update from lhs
		matrix.initMatrixInterfaces
		(
			coupleBouCoeffs,
			interfaces,
			x,
			bPrime,
			cmpt,
			true         // switch to lhs
		);

		// Update from lhs
		matrix.updateMatrixInterfaces
		(
			coupleBouCoeffs,
			interfaces,
			x,
			bPrime,
			cmpt,
			true         // switch to lhs
		);

		scalar curX;
		label fStart;
		label fEnd = ownStartPtr[0];

		for (label cellI = 0; cellI < nCells; cellI++)
		{
			// Start and end of this row
			fStart = fEnd;
			fEnd = ownStartPtr[cellI + 1];

			// Get the accumulated neighbour side
			curX = bPrimePtr[cellI];

			// Accumulate the owner product side
			for (label curFace = fStart; curFace < fEnd; curFace++)
			{
				curX -= upperPtr[curFace]*xPtr[uPtr[curFace]];
			}

			// Finish current x
			curX /= diagPtr[cellI];

			// Distribute the neighbour side using current x
			for (label curFace = fStart; curFace < fEnd; curFace++)
			{
				bPrimePtr[uPtr[curFace]] -= lowerPtr[curFace]*curX;
			}

			xPtr[cellI] = curX;
		}
	}
}


void Foam::GaussSeidelSmoother::smooth
(
	scalarField& x,
	const scalarField& b,
	const direction cmpt,
	const label nSweeps
) const
{
	smooth
	(
		x,
		matrix_,
		b,
		coupleBouCoeffs_,
		interfaces_,
		cmpt,
		nSweeps
	);
}


// ************************************************************************* //
