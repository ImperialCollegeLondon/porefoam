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
	Multiply a given vector (second argument) by the matrix or its transpose
	and return the result in the first argument.

\*---------------------------------------------------------------------------*/

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::lduMatrix::Amul
(
	scalarField& Ax,
	const scalarField& x,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const lduInterfaceFieldPtrsList& interfaces,
	const direction cmpt
) const
{
	// Reset multiplication result to zero
	// HJ, 5/Nov/2007
	Ax = 0;

	// Initialise the update of coupled interfaces
	initMatrixInterfaces
	(
		coupleBouCoeffs,
		interfaces,
		x,
		Ax,
		cmpt
	);

	// AmulCore must be additive to account for initialisation step
	// in ldu interfaces.  HJ, 6/Nov/2007
	AmulCore(Ax, x);

	// Update coupled interfaces
	updateMatrixInterfaces
	(
		coupleBouCoeffs,
		interfaces,
		x,
		Ax,
		cmpt
	);
}


void Foam::lduMatrix::AmulCore
(
	scalarField& Ax,
	const scalarField& x
) const
{
	scalar* __restrict__ AxPtr = Ax.begin();

	const scalar* const __restrict__ xPtr = x.begin();

	// Protection for multiplication of incomplete matrices
	// HJ, 19/Sep/2008
	if (hasDiag())
	{
		const scalar* const __restrict__ diagPtr = diag().begin();

		 const label nCells = diag().size();
		for ( label cell=0; cell<nCells; cell++)
		{
			// AmulCore must be additive to account for initialisation step
			// in ldu interfaces.  HJ, 6/Nov/2007
			AxPtr[cell] += diagPtr[cell]*xPtr[cell];
		}
	}

	// Protection for multiplication of incomplete matrices
	// HJ, 19/Sep/2008
	if (hasUpper() || hasLower())
	{
		const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
		const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

		const scalar* const __restrict__ upperPtr = upper().begin();
		const scalar* const __restrict__ lowerPtr = lower().begin();

		 const label nFaces = upper().size();

		for ( label face=0; face<nFaces; face++)
		{
			AxPtr[uPtr[face]] += lowerPtr[face]*xPtr[lPtr[face]];
			AxPtr[lPtr[face]] += upperPtr[face]*xPtr[uPtr[face]];
		}
	}
}


void Foam::lduMatrix::Tmul
(
	scalarField& Tx,
	const scalarField& x,
	const FieldField<Field, scalar>& coupleIntCoeffs,
	const lduInterfaceFieldPtrsList& interfaces,
	const direction cmpt
) const
{
	// Reset multiplication result to zero
	// HJ, 5/Nov/2007
	Tx = 0;

	// Initialise the update of coupled interfaces
	initMatrixInterfaces
	(
		coupleIntCoeffs,
		interfaces,
		x,
		Tx,
		cmpt
	);

	// TmulCore must be additive to account for initialisation step
	// in ldu interfaces.  HJ, 6/Nov/2007
	TmulCore(Tx, x);

	// Update coupled interfaces
	updateMatrixInterfaces
	(
		coupleIntCoeffs,
		interfaces,
		x,
		Tx,
		cmpt
	);
}


void Foam::lduMatrix::TmulCore
(
	scalarField& Tx,
	const scalarField& x
) const
{
	scalar* __restrict__ TxPtr = Tx.begin();

	const scalar* const __restrict__ xPtr = x.begin();

	// Protection for multiplication of incomplete matrices
	// HJ, 19/Sep/2008
	if (hasDiag())
	{
		const scalar* const __restrict__ diagPtr = diag().begin();

		 const label nCells = diag().size();
		for ( label cell=0; cell<nCells; cell++)
		{
			// TmulCore must be additive to account for initialisation step
			// in ldu interfaces.  HJ, 6/Nov/2007
			TxPtr[cell] += diagPtr[cell]*xPtr[cell];
		}
	}

	// Protection for multiplication of incomplete matrices
	// HJ, 19/Sep/2008
	if (hasUpper() || hasLower())
	{
		const label* const __restrict__ uPtr = lduAddr().upperAddr().begin();
		const label* const __restrict__ lPtr = lduAddr().lowerAddr().begin();

		const scalar* const __restrict__ lowerPtr = lower().begin();
		const scalar* const __restrict__ upperPtr = upper().begin();

		 const label nFaces = upper().size();
		for ( label face=0; face<nFaces; face++)
		{
			TxPtr[uPtr[face]] += upperPtr[face]*xPtr[lPtr[face]];
			TxPtr[lPtr[face]] += lowerPtr[face]*xPtr[uPtr[face]];
		}
	}
}


void Foam::lduMatrix::sumA
(
	scalarField& sumA,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const lduInterfaceFieldPtrsList& interfaces
) const
{
	scalar* __restrict__ sumAPtr = sumA.begin();

	const scalar* __restrict__ diagPtr = diag().begin();

	const label* __restrict__ uPtr = lduAddr().upperAddr().begin();
	const label* __restrict__ lPtr = lduAddr().lowerAddr().begin();

	const scalar* __restrict__ lowerPtr = lower().begin();
	const scalar* __restrict__ upperPtr = upper().begin();

	 const label nCells = diag().size();
	 const label nFaces = upper().size();

	for ( label cell=0; cell<nCells; cell++)
	{
		sumAPtr[cell] = diagPtr[cell];
	}

	for ( label face=0; face<nFaces; face++)
	{
		sumAPtr[uPtr[face]] += lowerPtr[face];
		sumAPtr[lPtr[face]] += upperPtr[face];
	}

	// Add the interface internal coefficients to diagonal
	// and the interface boundary coefficients to the sum-off-diagonal
	forAll(interfaces, patchI)
	{
		if (interfaces.set(patchI))
		{
			const unallocLabelList& pa = lduAddr().patchAddr(patchI);
			const scalarField& pCoeffs = coupleBouCoeffs[patchI];

			forAll(pa, face)
			{
				sumAPtr[pa[face]] -= pCoeffs[face];
			}
		}
	}
}


void Foam::lduMatrix::residual
(
	scalarField& rA,
	const scalarField& x,
	const scalarField& b,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const lduInterfaceFieldPtrsList& interfaces,
	const direction cmpt
) const
{
	// Reset multiplication result to zero
	// HJ, 5/Nov/2007
	rA = 0;

	// Standard implementation
	Amul(rA, x, coupleBouCoeffs, interfaces, cmpt);

	const scalar* const __restrict__ bPtr = b.begin();
	scalar* __restrict__ rAPtr = rA.begin();

	forAll (rA, cell)
	{
		rAPtr[cell] = bPtr[cell] - rAPtr[cell];
	}
}


Foam::tmp<Foam::scalarField> Foam::lduMatrix::residual
(
	const scalarField& x,
	const scalarField& b,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const lduInterfaceFieldPtrsList& interfaces,
	const direction cmpt
) const
{
	tmp<scalarField> trA(new scalarField(x.size()));
	residual(trA(), x, b, coupleBouCoeffs, interfaces, cmpt);
	return trA;
}


// ************************************************************************* //
