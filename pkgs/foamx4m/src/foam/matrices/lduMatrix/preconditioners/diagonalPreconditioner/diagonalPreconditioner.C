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

#include "diagonalPreconditioner.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(diagonalPreconditioner, 0);

	lduPreconditioner::
		addsymMatrixConstructorToTable<diagonalPreconditioner>
		adddiagonalPreconditionerSymMatrixConstructorToTable_;

	lduPreconditioner::
		addasymMatrixConstructorToTable<diagonalPreconditioner>
		adddiagonalPreconditionerAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::diagonalPreconditioner::diagonalPreconditioner
(
	const lduMatrix& matrix,
	const FieldField<Field, scalar>& coupleBouCoeffs,
	const FieldField<Field, scalar>& coupleIntCoeffs,
	const lduInterfaceFieldPtrsList& interfaces,
	const dictionary&
)
:
	lduPreconditioner
	(
		matrix,
		coupleBouCoeffs,
		coupleIntCoeffs,
		interfaces
	),
	rD(matrix.diag().size())
{
	scalar* __restrict__ rDPtr = rD.begin();
	const scalar* __restrict__ DPtr = matrix_.diag().begin();

	 label nCells = rD.size();

	// Generate reciprocal diagonal
	for ( label cell=0; cell<nCells; cell++)
	{
		rDPtr[cell] = 1.0/DPtr[cell];
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::diagonalPreconditioner::precondition
(
	scalarField& wA,
	const scalarField& rA,
	const direction
) const
{
	scalar* __restrict__ wAPtr = wA.begin();
	const scalar* __restrict__ rAPtr = rA.begin();
	const scalar* __restrict__ rDPtr = rD.begin();

	 label nCells = wA.size();

	for ( label cell=0; cell<nCells; cell++)
	{
		wAPtr[cell] = rDPtr[cell]*rAPtr[cell];
	}
}


// ************************************************************************* //
