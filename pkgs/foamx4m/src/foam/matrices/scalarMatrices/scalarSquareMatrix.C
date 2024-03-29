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

#include "scalarSquareMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::scalarSquareMatrix::scalarSquareMatrix()
{}


Foam::scalarSquareMatrix::scalarSquareMatrix(const label mSize)
:
	SquareMatrix<scalar>(mSize, 0.0)
{}


Foam::scalarSquareMatrix::scalarSquareMatrix
(
	const label mSize,
	const scalar v
)
:
	SquareMatrix<scalar>(mSize, v)
{}


// Foam::scalarSquareMatrix::scalarSquareMatrix(const scalarSquareMatrix& matrix)
// :
//     SquareMatrix<scalar>(matrix)
// {}


Foam::scalarSquareMatrix::scalarSquareMatrix(const SquareMatrix<scalar>& matrix)
:
	SquareMatrix<scalar>(matrix)
{}


Foam::scalarSquareMatrix::scalarSquareMatrix(Istream& is)
:
	SquareMatrix<scalar>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::scalarSquareMatrix::LUDecompose
(
	scalarSquareMatrix& matrix,
	labelList& pivotIndices
)
{
	label n = matrix.n();
	scalar vv[n];

	for ( label i=0; i<n; i++)
	{
		scalar largestCoeff = 0.0;
		scalar temp;
		const scalar* __restrict__ matrixi = matrix[i];

		for ( label j=0; j<n; j++)
		{
			if ((temp = mag(matrixi[j])) > largestCoeff)
			{
				largestCoeff = temp;
			}
		}

		if (largestCoeff == 0.0)
		{
			FatalErrorIn
			(
				"scalarSquareMatrix::LUdecompose"
				"(scalarSquareMatrix& matrix, labelList& rowIndices)"
			)   << "Singular matrix" << exit(FatalError);
		}

		vv[i] = 1.0/largestCoeff;
	}

	for ( label j=0; j<n; j++)
	{
		scalar* __restrict__ matrixj = matrix[j];

		for ( label i=0; i<j; i++)
		{
			scalar* __restrict__ matrixi = matrix[i];

			scalar sum = matrixi[j];
			for ( label k=0; k<i; k++)
			{
				sum -= matrixi[k]*matrix[k][j];
			}
			matrixi[j] = sum;
		}

		label iMax = 0;

		scalar largestCoeff = 0.0;
		for ( label i=j; i<n; i++)
		{
			scalar* __restrict__ matrixi = matrix[i];
			scalar sum = matrixi[j];

			for ( label k=0; k<j; k++)
			{
				sum -= matrixi[k]*matrix[k][j];
			}

			matrixi[j] = sum;

			scalar temp;
			if ((temp = vv[i]*mag(sum)) >= largestCoeff)
			{
				largestCoeff = temp;
				iMax = i;
			}
		}

		pivotIndices[j] = iMax;

		if (j != iMax)
		{
			scalar* __restrict__ matrixiMax = matrix[iMax];

			for ( label k=0; k<n; k++)
			{
				Swap(matrixj[k], matrixiMax[k]);
			}

			vv[iMax] = vv[j];
		}

		if (matrixj[j] == 0.0)
		{
			matrixj[j] = SMALL;
		}

		if (j != n-1)
		{
			scalar rDiag = 1.0/matrixj[j];

			for ( label i=j+1; i<n; i++)
			{
				matrix[i][j] *= rDiag;
			}
		}
	}
}


Foam::scalarSquareMatrix Foam::scalarSquareMatrix::LUinvert() const
{
	scalarSquareMatrix luMatrix = *this;

	scalarSquareMatrix luInvert(luMatrix.n());
	scalarField column(luMatrix.n());

	labelList pivotIndices(luMatrix.n());

	LUDecompose(luMatrix, pivotIndices);

	for (label j = 0; j < luMatrix.n(); j++)
	{
		for (label i = 0; i < luMatrix.n(); i++)
		{
			column[i] = 0.0;
		}

		column[j] = 1.0;

		LUBacksubstitute(luMatrix, pivotIndices, column);

		for (label i = 0; i < luMatrix.n(); i++)
		{
			luInvert[i][j] = column[i];
		}
	}

	return luInvert;
}


// ************************************************************************* //
