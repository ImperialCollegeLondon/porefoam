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
	A storage mechanism which allows setting of the fixed value and
	consequently recovering the equation for a single row of the matrix as
	well as the source. The equation is taken out of the matrix using a
	variant of compact matrix storage mechanism.

\*---------------------------------------------------------------------------*/

#include "constraint.H"
#include "demandDrivenData.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
scalar constraint<Type>::componentOfValue
(
	const Type& v,
	const direction d
) const
{
	return v.component(d);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
template<class Type>
constraint<Type>::constraint
(
	const label row,
	const Type value,
	const Type& fixedCmpts
)
:
	rowID_(row),
	value_(value),
	fixedComponents_(fixedCmpts),
	matrixCoeffsSet_(false),
	upperCoeffsOwnerPtr_(nullptr),
	upperCoeffsNeighbourPtr_(nullptr),
	lowerCoeffsOwnerPtr_(nullptr),
	lowerCoeffsNeighbourPtr_(nullptr)
{}


// Construct from matrix and components
template<class Type>
template<template<class> class Matrix>
constraint<Type>::constraint
(
	const Matrix<Type>& matrix,
	const label row,
	const Type value,
	const Type& fixedCmpts
)
:
	rowID_(row),
	value_(value),
	fixedComponents_(fixedCmpts),
	matrixCoeffsSet_(false),
	upperCoeffsOwnerPtr_(nullptr),
	upperCoeffsNeighbourPtr_(nullptr),
	lowerCoeffsOwnerPtr_(nullptr),
	lowerCoeffsNeighbourPtr_(nullptr)
{
	setMatrix(matrix);
}


// Construct as copy
template<class Type>
constraint<Type>::constraint(const constraint& e)
:
	rowID_(e.rowID_),
	value_(e.value_),
	fixedComponents_(e.fixedComponents_),
	matrixCoeffsSet_(false),
	upperCoeffsOwnerPtr_(nullptr),
	upperCoeffsNeighbourPtr_(nullptr),
	lowerCoeffsOwnerPtr_(nullptr),
	lowerCoeffsNeighbourPtr_(nullptr)
{}


// Construct from Istream
template<class Type>
constraint<Type>::constraint(Istream& is)
:
	rowID_(readLabel(is)),
	value_(is),
	fixedComponents_(is),
	matrixCoeffsSet_(false),
	upperCoeffsOwnerPtr_(nullptr),
	upperCoeffsNeighbourPtr_(nullptr),
	lowerCoeffsOwnerPtr_(nullptr),
	lowerCoeffsNeighbourPtr_(nullptr)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
constraint<Type>::~constraint()
{
	deleteDemandDrivenData(upperCoeffsOwnerPtr_);
	deleteDemandDrivenData(upperCoeffsNeighbourPtr_);

	deleteDemandDrivenData(lowerCoeffsOwnerPtr_);
	deleteDemandDrivenData(lowerCoeffsNeighbourPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
scalar constraint<Type>::diagCoeff() const
{
	if (matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"scalar constraint<Type>::diagCoeff() const"
		)   << "matrix coefficients not set"
			<< abort(FatalError);
	}

	return diagCoeff_;
}


template<class Type>
Type constraint<Type>::source() const
{
	if (matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"Type constraint<Type>::source() const"
		)   << "matrix coefficients not set"
			<< abort(FatalError);
	}

	return source_;
}


template<class Type>
const scalarField& constraint<Type>::upperCoeffsOwner() const
{
	if (!upperCoeffsOwnerPtr_ || !matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"const scalarField& constraint<Type>::"
			"upperCoeffsOwner() const"
		)   << "upper matrix coefficients not set"
			<< abort(FatalError);
	}

	return *upperCoeffsOwnerPtr_;
}


template<class Type>
const scalarField& constraint<Type>::upperCoeffsNeighbour() const
{
	if (!upperCoeffsNeighbourPtr_ || !matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"const scalarField& constraint<Type>::"
			"upperCoeffsNeighbour() const"
		)   << "upper matrix coefficients not set"
			<< abort(FatalError);
	}

	return *upperCoeffsNeighbourPtr_;
}


template<class Type>
const scalarField& constraint<Type>::lowerCoeffsOwner() const
{
	if (!lowerCoeffsOwnerPtr_ || !matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"const scalarField& constraint<Type>::"
			"lowerCoeffsOwner() const"
		)   << "lower matrix coefficients not set"
			<< abort(FatalError);
	}

	return *lowerCoeffsOwnerPtr_;
}


template<class Type>
const scalarField& constraint<Type>::lowerCoeffsNeighbour() const
{
	if (!lowerCoeffsNeighbourPtr_ || !matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"const scalarField& constraint<Type>::"
			"lowerCoeffsNeighbour() const"
		)   << "lower matrix coefficients not set"
			<< abort(FatalError);
	}

	return *lowerCoeffsNeighbourPtr_;
}


template<class Type>
void constraint<Type>::combine
(
	const constraint<Type>& e
)
{
	for
	(
		direction cmptI = 0;
		cmptI < pTraits<Type>::nComponents;
		cmptI++
	)
	{
		if
		(
			e.fixedComponents_.component(cmptI)
		  > fixedComponents_.component(cmptI)
		)
		{
			fixedComponents_.component(cmptI) =
				e.fixedComponents_.component(cmptI);

			value_.replace(cmptI, e.value_.component(cmptI));
		}
	}
}


template<class Type>
void constraint<Type>::clearMatrix()
{
	matrixCoeffsSet_ = false;

	diagCoeff_ = 0;

	source_ = pTraits<Type>::zero;

	deleteDemandDrivenData(upperCoeffsOwnerPtr_);
	deleteDemandDrivenData(upperCoeffsNeighbourPtr_);

	deleteDemandDrivenData(lowerCoeffsOwnerPtr_);
	deleteDemandDrivenData(lowerCoeffsNeighbourPtr_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Type>
void constraint<Type>::operator=
(
	const constraint<Type>& rhs
)
{
	// Check for assignment to self
	if (this == &rhs)
	{
		FatalErrorIn
		(
			"constraint::operator=(const constraint&)"
		)   << "attempted assignment to self"
			<< abort(FatalError);
	}

	rowID_ = rhs.rowID_;

	value_ = rhs.value_;

	fixedComponents_ = rhs.fixedComponents_;

	clearMatrix();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class Type>
Ostream& operator<<(Ostream& os, const constraint<Type>& e)
{
	os  << e.rowID_ << nl
		<< e.value_ << nl
		<< e.fixedComponents_ << nl;

	os.check("Ostream& operator<<(Ostream&, constraint<Type>&");

	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
