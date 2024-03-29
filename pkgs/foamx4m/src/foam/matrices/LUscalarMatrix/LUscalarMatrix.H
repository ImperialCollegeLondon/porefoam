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
	Foam::LUscalarMatrix

Description
	Foam::LUscalarMatrix

SourceFiles
	LUscalarMatrix.C

\*---------------------------------------------------------------------------*/

#ifndef LUscalarMatrix_H
#define LUscalarMatrix_H

#include "scalarMatrices.H"
#include "labelList.H"
#include "FieldField.H"
#include "lduInterfaceFieldPtrsList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class lduMatrix;
class procLduMatrix;


class LUscalarMatrix
:
	public scalarSquareMatrix
{
	// Private data

		//- Processor matrix offsets
		labelList procOffsets_;

		//- The pivot indices used in the LU decomposition
		labelList pivotIndices_;

		//- Convert the given lduMatrix into this LUscalarMatrix
		void convert
		(
			const lduMatrix& ldum,
			const FieldField<Field, scalar>& interfaceCoeffs,
			const lduInterfaceFieldPtrsList& interfaces
		);

		//- Convert the given list of procLduMatrix into this LUscalarMatrix
		//  on the master processor
		void convert(const PtrList<procLduMatrix>& lduMatrices);


		//- Print the ratio of the mag-sum of the off-diagonal coefficients
		//  to the mag-diagonal
		void printDiagonalDominance() const;


public:

	// Constructors

		//- Construct from scalarSquareMatrix and perform LU decomposition
		LUscalarMatrix(const scalarSquareMatrix&);

		//- Construct from lduMatrix and perform LU decomposition
		LUscalarMatrix
		(
			const lduMatrix&,
			const FieldField<Field, scalar>& interfaceCoeffs,
			const lduInterfaceFieldPtrsList& interfaces
		);


	// Member Functions

		//- Solve the matrix using the LU decomposition with pivoting
		//  returning the solution in the source
		template<class T>
		void solve(Field<T>& source) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "LUscalarMatrixTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
