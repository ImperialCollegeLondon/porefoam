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

#include "simpleMatrix.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::simpleMatrix<Type>::simpleMatrix(const label mSize)
:
	scalarSquareMatrix(mSize),
	source_(mSize)
{}


template<class Type>
Foam::simpleMatrix<Type>::simpleMatrix
(
	const label mSize,
	const scalar coeffVal,
	const Type& sourceVal
)
:
	scalarSquareMatrix(mSize, coeffVal),
	source_(mSize, sourceVal)
{}


template<class Type>
Foam::simpleMatrix<Type>::simpleMatrix
(
	const scalarSquareMatrix& matrix,
	const Field<Type>& source
)
:
	scalarSquareMatrix(matrix),
	source_(source)
{}


template<class Type>
Foam::simpleMatrix<Type>::simpleMatrix(Istream& is)
:
	scalarSquareMatrix(is),
	source_(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Foam::Field<Type> Foam::simpleMatrix<Type>::solve() const
{
	// Since matrix and source are trashed during solution,
	// a copy is made.  HJ, 23/Dec/2008
	scalarSquareMatrix tmpMatrix = *this;
	Field<Type> sourceSol = source_;

	scalarSquareMatrix::solve(tmpMatrix, sourceSol);

	return sourceSol;
}


template<class Type>
Foam::Field<Type> Foam::simpleMatrix<Type>::LUsolve() const
{
	// Since matrix and source are trashed during solution,
	// a copy is made.  HJ, 23/Dec/2008
	scalarSquareMatrix luMatrix = *this;
	Field<Type> sourceSol = source_;

	scalarSquareMatrix::LUsolve(luMatrix, sourceSol);

	return sourceSol;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void Foam::simpleMatrix<Type>::operator=(const simpleMatrix<Type>& m)
{
	if (this == &m)
	{
		FatalErrorIn("simpleMatrix<Type>::operator=(const simpleMatrix<Type>&)")
			<< "Attempted assignment to self"
			<< abort(FatalError);
	}

	if (n() != m.n())
	{
		FatalErrorIn("simpleMatrix<Type>::operator=(const simpleMatrix<Type>&)")
			<< "Different size matrices"
			<< abort(FatalError);
	}

	if (source_.size() != m.source_.size())
	{
		FatalErrorIn("simpleMatrix<Type>::operator=(const simpleMatrix<Type>&)")
			<< "Different size source vectors"
			<< abort(FatalError);
	}

	scalarSquareMatrix::operator=(m);
	source_ = m.source_;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Type>
Foam::simpleMatrix<Type> Foam::operator+
(
	const simpleMatrix<Type>& m1,
	const simpleMatrix<Type>& m2
)
{
	return simpleMatrix<Type>
	(
		static_cast<const scalarSquareMatrix&>(m1)
	  + static_cast<const scalarSquareMatrix&>(m2),
		m1.source_ + m2.source_
	);
}


template<class Type>
Foam::simpleMatrix<Type> Foam::operator-
(
	const simpleMatrix<Type>& m1,
	const simpleMatrix<Type>& m2
)
{
	return simpleMatrix<Type>
	(
		static_cast<const scalarSquareMatrix&>(m1)
	  - static_cast<const scalarSquareMatrix&>(m2),
		m1.source_ - m2.source_
	);
}


template<class Type>
Foam::simpleMatrix<Type> Foam::operator*
(
	const scalar s,
	const simpleMatrix<Type>& m
)
{
	return simpleMatrix<Type>(s*m.matrix_, s*m.source_);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Foam::Ostream& Foam::operator<<(Ostream& os, const simpleMatrix<Type>& m)
{
	os << static_cast<const scalarSquareMatrix&>(m) << nl << m.source_;
	return os;
}


// ************************************************************************* //
