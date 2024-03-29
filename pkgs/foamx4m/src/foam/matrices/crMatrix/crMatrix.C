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
	crMatrix

Description
	Sparse matrix in compressed row format

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved

\*----------------------------------------------------------------------------*/

#include "crMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::crMatrix::crMatrix()
:
	refCount(),
	crAddr_(),
	coeffs_()
{}


// Construct from addressing
Foam::crMatrix::crMatrix(const crAddressing& addr)
:
	refCount(),
	crAddr_(addr),
	coeffs_(addr.nEntries())
{}


// Construct given size.  Column and coefficients set later
Foam::crMatrix::crMatrix
(
	const label nRows,
	const label nCols,
	const labelList& count
)
:
	refCount(),
	crAddr_(nRows, nCols, count),
	coeffs_(crAddr_.nEntries())
{}


// Construct from components of addressing
Foam::crMatrix::crMatrix
(
	const label nRows,
	const label nCols,
	const labelList& row,
	const labelList& col
)
:
	refCount(),
	crAddr_(nRows, nCols, row, col),
	coeffs_(crAddr_.nEntries())
{}

// Construct as copy
Foam::crMatrix::crMatrix(const crMatrix& m)
:
	refCount(),
	crAddr_(m.crAddr_),
	coeffs_(m.coeffs_)
{}


// Construct as copy of tmp<crMatrix> deleting argument
Foam::crMatrix::crMatrix(const tmp<crMatrix>& tm)
:
	refCount(),
	crAddr_(tm().crAddr()),
	coeffs_(tm().coeffs())
{
	tm.clear();
}


// Construct from Istream
Foam::crMatrix::crMatrix(Istream& is)
:
	refCount(),
	crAddr_(is),
	coeffs_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return transpose
Foam::tmp<Foam::crMatrix> Foam::crMatrix::T() const
{
	// My addressing
	const label myNRows = crAddr().nRows();
	const labelList& myRow = crAddr().rowStart();
	const labelList& myCol = crAddr().column();
	const scalarField& myCoeffs = coeffs();

	// Create transpose
	tmp<crMatrix> ttranspose(new crMatrix(crAddr().T()));
	crMatrix& transpose = ttranspose();

	scalarField& tCoeffs = transpose.coeffs();

	// Reset coefficients to zero to treat degenerate matrices
	tCoeffs = 0;

	// Transpose addressing
	const labelList& tRow = transpose.crAddr().rowStart();
	const labelList& tCol = transpose.crAddr().column();

	// Set transpose coefficients

	label ja;

	for (label ia = 0; ia < myNRows; ia++)
	{
		for (label jpa = myRow[ia]; jpa < myRow[ia + 1]; jpa++)
		{
			ja = myCol[jpa];

			for (label jpb = tRow[ja]; jpb < tRow[ja + 1]; jpb++)
			{
				if (tCol[jpb] == ia)
				{
					tCoeffs[jpb] = myCoeffs[jpa];
					break;
				}
			}
		}
	}

	return ttranspose;
}

// Calculate b += A*x
void Foam::crMatrix::dotPlus(scalarField& b, const scalarField& x) const
{
	const labelList& row = crAddr_.rowStart();
	const labelList& col = crAddr_.column();

	forAll (b, i)
	{
		scalar& bi = b[i];

		for (label ip = row[i]; ip < row[i + 1];  ip++)
		{
			bi += coeffs_[ip]*x[col[ip]];
		}
	}
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::crMatrix::operator=(const crMatrix& rhs)
{
	// Check for assignment to self
	if (this == &rhs)
	{
		FatalErrorIn("Foam::crMatrix::operator=(const Foam::crMatrix&)")
			<< "Attempted assignment to self"
			<< abort(FatalError);
	}

	crAddr_ = rhs.crAddr_;
	coeffs_ = rhs.coeffs_;
}


void Foam::crMatrix::operator=(const tmp<crMatrix>& trhs)
{
	operator=(trhs());
	trhs.clear();
}


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const crMatrix& m)
{
	os << m.crAddr_ << m.coeffs_ << endl;

	return os;
}


// ************************************************************************* //
