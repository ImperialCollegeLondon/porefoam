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
	Foam::SVD

Description
	Singular value decomposition of a rectangular matrix.

SourceFiles
	SVDI.H
	SVD.C

\*---------------------------------------------------------------------------*/

#ifndef SVD_H
#define SVD_H

#include "scalarMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



class SVD
{
	// Private data

		//- Rectangular matrix with the same dimensions as the input
		scalarRectangularMatrix U_;

		//- square matrix V
		scalarRectangularMatrix V_;

		//- The singular values
		DiagonalMatrix<scalar> S_;

		//- The matrix product V S^(-1) U^T
		scalarRectangularMatrix VSinvUt_;

		//- The number of zero singular values
		label nZeros_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		SVD(const SVD&);

		//- Disallow default bitwise assignment
		void operator=(const SVD&);

		template<class T>
		inline const T sign(const T& a, const T& b);


public:

	// Constructors

		//- Construct from a rectangular Matrix
		SVD(const scalarRectangularMatrix& A, const scalar minCondition = 0);


	// Access functions

		//- Return U
		inline const scalarRectangularMatrix& U() const;

		//- Return the square matrix V
		inline const scalarRectangularMatrix& V() const;

		//- Return the singular values
		inline const scalarDiagonalMatrix& S() const;

		//- Return VSinvUt (the pseudo inverse)
		inline const scalarRectangularMatrix& VSinvUt() const;

		//- Return the number of zero singular values
		inline label nZeros() const;

		//- Return the minimum non-zero singular value
		inline scalar minNonZeroS() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SVDI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
