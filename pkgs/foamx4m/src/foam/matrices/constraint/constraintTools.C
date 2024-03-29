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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
template<template<class> class Matrix>
void constraint<Type>::setMatrix
(
	const Matrix<Type>& matrix
)
{
	if (matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"const scalarField& constraint<Type>::setMatrix"
		)   << "(const Matrix<Type>& matrix)"
			<< "matrix coefficients already set"
			<< abort(FatalError);
	}

	matrixCoeffsSet_ = true;

	if (matrix.hasDiag())
	{
		diagCoeff_ = matrix.diag()[rowID_];
	}

	source_ = matrix.source()[rowID_];

	const label startFaceOwn =
		matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID_];
	const label endFaceOwn =
		matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID_ + 1];

	const label startFaceNbr =
		matrix.psi().mesh().lduAddr().losortStartAddr()[rowID_];
	const label endFaceNbr =
		matrix.psi().mesh().lduAddr().losortStartAddr()[rowID_ + 1];

	const unallocLabelList& losort = matrix.psi().mesh().lduAddr().losortAddr();

	if (matrix.hasUpper())
	{
		// get the upper coefficients

		const scalarField& matrixUpper = matrix.upper();

		// owner side
		upperCoeffsOwnerPtr_ = new scalarField(endFaceOwn - startFaceOwn);

		scalarField& uOwn = *upperCoeffsOwnerPtr_;

		label faceIndex = startFaceOwn;

		forAll(uOwn, uOwnI)
		{
			uOwn[uOwnI] = matrixUpper[faceIndex];

			faceIndex++;
		}

		// neighbour side
		upperCoeffsNeighbourPtr_ = new scalarField(endFaceNbr - startFaceNbr);

		scalarField& uNbr = *upperCoeffsNeighbourPtr_;

		faceIndex = startFaceNbr;

		forAll (uNbr, uNbrI)
		{
			uNbr[uNbrI] = matrixUpper[losort[faceIndex]];

			faceIndex++;
		}
	}

	if (matrix.hasLower())
	{
		// get the lower coefficients

		const scalarField& matrixLower = matrix.lower();

		// owner side
		lowerCoeffsOwnerPtr_ = new scalarField(endFaceOwn - startFaceOwn);

		scalarField& lOwn = *lowerCoeffsOwnerPtr_;

		label faceIndex = startFaceOwn;

		forAll (lOwn, lOwnI)
		{
			lOwn[lOwnI] = matrixLower[faceIndex];

			faceIndex++;
		}

		// neighbour side
		lowerCoeffsNeighbourPtr_ = new scalarField(endFaceNbr - startFaceNbr);

		scalarField& lNbr = *lowerCoeffsNeighbourPtr_;

		faceIndex = startFaceNbr;

		forAll (lNbr, lNbrI)
		{
			lNbr[lNbrI] = matrixLower[losort[faceIndex]];

			faceIndex++;
		}
	}
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::eliminateEquation
(
	Matrix<Type>& matrix,
	const label rowID,
	const Type& value
)
{
	// Record equation as eliminated
	matrix.eliminatedEqns().insert(rowID);

	Field<Type>& source = matrix.source();

	const label startFaceOwn =
		matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID];
	const label endFaceOwn =
		matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID + 1];

	const label startFaceNbr =
		matrix.psi().mesh().lduAddr().losortStartAddr()[rowID];
	const label endFaceNbr =
		matrix.psi().mesh().lduAddr().losortStartAddr()[rowID + 1];

	const unallocLabelList& owner = matrix.psi().mesh().lduAddr().lowerAddr();
	const unallocLabelList& neighbour =
		matrix.psi().mesh().lduAddr().upperAddr();
	const unallocLabelList& losort =
		matrix.psi().mesh().lduAddr().losortAddr();

	// My index =  rowID
	if (matrix.symmetric())
	{
		// get the coefficients
		scalarField* coeffs = nullptr;

		if (matrix.hasUpper())
		{
			coeffs = &matrix.upper();
		}
		else if (matrix.hasLower())
		{
			coeffs = &matrix.lower();
		}

		scalarField& matrixCoeffs = *coeffs;

		// owner side
		label faceIndex = startFaceOwn;

		while (faceIndex < endFaceOwn)
		{
			// add contribution to source of the neighbour (I am the owner)
			source[neighbour[faceIndex]] -= matrixCoeffs[faceIndex]*value;

			// I am the owner =  owner[faceIndex]
			// set off-diagonal coefficient to zero
			matrixCoeffs[faceIndex] = 0.0;

			faceIndex++;
		}

		// neighbour side
		faceIndex = startFaceNbr;

		while (faceIndex < endFaceNbr)
		{
			// add contribution to source of owner (I am the neighbour)
			source[owner[losort[faceIndex]]] -=
				matrixCoeffs[losort[faceIndex]]*value;

			// I am the neighbour = neighbour[losort[faceIndex]]
			// set off-diagonal coefficient to zero
			matrixCoeffs[losort[faceIndex]] = 0.0;

			faceIndex++;
		}
	}
	else if (matrix.asymmetric())
	{
		scalarField& matrixUpper = matrix.upper();

		scalarField& matrixLower = matrix.lower();

		// owner side
		label faceIndex = startFaceOwn;

		while (faceIndex < endFaceOwn)
		{
			// add contribution to source of the neighbour (I am the owner)
			source[neighbour[faceIndex]] -= matrixLower[faceIndex]*value;

			// set off-diagonal coefficient to zero
			matrixLower[faceIndex] = 0.0;

			faceIndex++;
		}

		// neighbour side
		faceIndex = startFaceNbr;

		while (faceIndex < endFaceNbr)
		{
			// add contribution to source of owner (I am the neighbour)
			source[owner[losort[faceIndex]]] -=
				matrixUpper[losort[faceIndex]]*value;

			// set off-diagonal coefficient to zero
			matrixUpper[losort[faceIndex]] = 0.0;

			faceIndex++;
		}
	}
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::eliminateEquation(Matrix<Type>& matrix) const
{
	eliminateEquation(matrix, rowID_, value_);
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::eliminateEquation
(
	Matrix<Type>& matrix,
	const direction d,
	scalarField& sourceCmpt
) const
{
	// Record equation as eliminated
	matrix.eliminatedEqns().insert(rowID_);

	const Type& fc = fixedComponents();

	const scalar fcOfD = componentOfValue(fc, d);

	if (fcOfD > SMALL)
	{
		const label startFaceOwn =
			matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID_];
		const label endFaceOwn =
			matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID_ + 1];

		const label startFaceNbr =
			matrix.psi().mesh().lduAddr().losortStartAddr()[rowID_];
		const label endFaceNbr =
			matrix.psi().mesh().lduAddr().losortStartAddr()[rowID_ + 1];

		const unallocLabelList& owner =
			matrix.psi().mesh().lduAddr().lowerAddr();
		const unallocLabelList& neighbour =
			matrix.psi().mesh().lduAddr().upperAddr();
		const unallocLabelList& losort =
			matrix.psi().mesh().lduAddr().losortAddr();

		// My index =  rowID_
		if (matrix.symmetric())
		{
			// get the coefficients
			scalarField* coeffs = nullptr;

			if (matrix.hasUpper())
			{
				coeffs = &matrix.upper();
			}
			else if (matrix.hasLower())
			{
				coeffs = &matrix.lower();
			}

			scalarField& matrixCoeffs = *coeffs;

			// owner side
			label faceIndex = startFaceOwn;

			while (faceIndex < endFaceOwn)
			{
				// add contribution to source of the neighbour (I am the owner)
				sourceCmpt[neighbour[faceIndex]] -=
					fcOfD*matrixCoeffs[faceIndex]*componentOfValue(value(), d);

				// I am the owner =  owner[faceIndex]
				// set off-diagonal coefficient to zero
				matrixCoeffs[faceIndex] *= 1.0 - fcOfD;

				faceIndex++;
			}

			// neighbour side
			faceIndex = startFaceNbr;

			while (faceIndex < endFaceNbr)
			{
				// add contribution to source of owner (I am the neighbour)
				sourceCmpt[owner[losort[faceIndex]]] -=
					fcOfD*matrixCoeffs[losort[faceIndex]]
					*componentOfValue(value(), d);

				// I am the neighbour = neighbour[losort[faceIndex]]
				// set off-diagonal coefficient to zero
				matrixCoeffs[losort[faceIndex]] *= 1.0 - fcOfD;

				faceIndex++;
			}
		}
		else if (matrix.asymmetric())
		{
			scalarField& matrixUpper = matrix.upper();

			scalarField& matrixLower = matrix.lower();

			// owner side
			label faceIndex = startFaceOwn;

			while (faceIndex < endFaceOwn)
			{
				// add contribution to source of the neighbour (I am the owner)
				sourceCmpt[neighbour[faceIndex]] -=
					fcOfD*matrixLower[faceIndex]*componentOfValue(value(), d);

				// set off-diagonal coefficient to zero
				matrixLower[faceIndex] *= 1.0 - fcOfD;

				faceIndex++;
			}

			// neighbour side
			faceIndex = startFaceNbr;

			while (faceIndex < endFaceNbr)
			{
				// add contribution to source of owner (I am the neighbour)
				sourceCmpt[owner[losort[faceIndex]]] -=
					fcOfD*matrixUpper[losort[faceIndex]]
					*componentOfValue(value(), d);

				// set off-diagonal coefficient to zero
				matrixUpper[losort[faceIndex]] *= 1.0 - fcOfD;

				faceIndex++;
			}
		}
	}
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::setSource
(
	Matrix<Type>& matrix,
	const label rowID,
	const Type& value
)
{
	matrix.source()[rowID] = matrix.diag()[rowID]*value;
	matrix.psi()[rowID] = value;
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::setSource(Matrix<Type>& matrix) const
{
	setSource(matrix, rowID_, value_);
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::setSourceDiag
(
	Matrix<Type>& matrix,
	const direction d,
	scalarField& psiCmpt,
	scalarField& sourceCmpt
) const
{
	const Type& fc = fixedComponents();

	if (componentOfValue(fc, d) > SMALL)
	{
		sourceCmpt[rowID()] = componentOfValue(fc, d)*matrix.diag()[rowID()]
			*componentOfValue(value(), d);

		// set the solution to the right value as well
		psiCmpt[rowID()] =
			componentOfValue(fc, d)*componentOfValue(value(), d);
	}
}


template<class Type>
template<template<class> class Matrix>
void constraint<Type>::reconstructMatrix
(
	Matrix<Type>& matrix
) const
{
	if (!matrixCoeffsSet_)
	{
		FatalErrorIn
		(
			"void constraint<Type>::reconstructMatrix("
			"Matrix<Type>& matrix)"
		)   << "matrix coefficients not set"
			<< abort(FatalError);
	}

	if (matrix.hasDiag())
	{
		matrix.diag()[rowID_] = diagCoeff_;
	}

	const label startFaceOwn =
		matrix.psi().mesh().lduAddr().ownerStartAddr()[rowID_];

	const label startFaceNbr =
		matrix.psi().mesh().lduAddr().losortStartAddr()[rowID_];

	const unallocLabelList& losort = matrix.psi().mesh().lduAddr().losortAddr();

	if (matrix.hasUpper())
	{
		// get the upper coefficients

		scalarField& matrixUpper = matrix.upper();

		// owner side
		const scalarField& uOwn = upperCoeffsOwner();

		label faceIndex = startFaceOwn;

		forAll(uOwn, uOwnI)
		{
			matrixUpper[faceIndex] = uOwn[uOwnI];

			faceIndex++;
		}

		// neighbour side
		const scalarField& uNbr = upperCoeffsNeighbour();

		faceIndex = startFaceNbr;

		forAll (uNbr, uNbrI)
		{
			matrixUpper[losort[faceIndex]] = uNbr[uNbrI];

			faceIndex++;
		}
	}

	if (matrix.hasLower())
	{
		// get the lower coefficients

		scalarField& matrixLower = matrix.lower();

		// owner side
		const scalarField& lOwn = lowerCoeffsOwner();

		label faceIndex = startFaceOwn;

		forAll (lOwn, lOwnI)
		{
			matrixLower[faceIndex] = lOwn[lOwnI];

			faceIndex++;
		}

		// neighbour side
		const scalarField& lNbr = lowerCoeffsNeighbour();

		faceIndex = startFaceNbr;

		forAll (lNbr, lNbrI)
		{
			matrixLower[losort[faceIndex]] = lNbr[lNbrI];

			faceIndex++;
		}
	}
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
