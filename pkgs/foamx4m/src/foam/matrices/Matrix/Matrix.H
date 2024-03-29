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
	Foam::Matrix

Description
	A templated 2D matrix of objects of \<T\>, where the n x m matrix
	dimensions are known and used for subscript bounds checking, etc.

SourceFiles
	Matrix.C
	MatrixI.H
	MatrixIO.C

\*---------------------------------------------------------------------------*/

#ifndef Matrix_H
#define Matrix_H

#include "List.H"
#include "label.H"
#include "bool.H"
#include "autoPtr.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Form, class Type> class Matrix;

template<class Form, class Type> Istream& operator>>
(
	Istream&,
	Matrix<Form, Type>&
);

template<class Form, class Type> Ostream& operator<<
(
	Ostream&,
	const Matrix<Form, Type>&
);



template<class Form, class Type>
class Matrix
{
	// Private data

		//- Row pointers
		Type** __restrict__ v_;

		//- Number of rows and columns in Matrix.
		label n_, m_;

		//- Allocate the storage for the row-pointers and the data
		//  and set the row pointers
		void allocate();


public:

	// Static data members

		//- Empty matrix
		static const Matrix<Form, Type> zero;


	// Static Member Functions

		//- Return a null Matrix
		inline static const Matrix<Form, Type>& null();


	// Constructors

		//- Null constructor.
		inline Matrix();

		//- Construct given number of rows and columns.
		Matrix(const label n, const label m);

		//- Construct with given number of rows and columns
		//  and value for all elements.
		Matrix(const label n, const label m, const Type&);

		//- Copy constructor.
	Matrix(const Matrix<Form, Type>&);

		//- Construct from Istream.
		Matrix(Istream&);

		//- Clone
		inline autoPtr<Matrix<Form, Type> > clone() const;


	// Destructor

		~Matrix();


	// Member functions

		// Access

			//- Return the number of rows
			inline label n() const;

			//- Return the number of columns
			inline label m() const;

			//- Return the number of elements in matrix (n*m)
			inline label size() const;


		// Check

			//- Check index i is within valid range (0 ... n-1).
			inline void checki(const label i) const;

			//- Check index j is within valid range (0 ... m-1).
			inline void checkj(const label j) const;


		// Edit

			//- Clear the Matrix, i.e. set sizes to zero.
			void clear();

			//- Transfer the contents of the argument Matrix into this Matrix
			//  and annull the argument Matrix.
			void transfer(Matrix<Form, Type>&);


		//- Return the transpose of the matrix
		Form T() const;


	// Member operators

		//- Return subscript-checked element of Matrix.
		inline Type* operator[](const label);

		//- Return subscript-checked element of constant Matrix.
		inline const Type* operator[](const label) const;

		//- Assignment operator. Takes linear time.
		void operator=(const Matrix<Form, Type>&);

		//- Assignment of all entries to the given value
		void operator=(const Type&);


	// IOstream operators

		//- Read Matrix from Istream, discarding contents of existing Matrix.
		friend Istream& operator>> <Form, Type>
		(
			Istream&,
			Matrix<Form, Type>&
		);

		// Write Matrix to Ostream.
		friend Ostream& operator<< <Form, Type>
		(
			Ostream&,
			const Matrix<Form, Type>&
		);
};


// Global functions and operators

template<class Form, class Type> const Type& max(const Matrix<Form, Type>&);
template<class Form, class Type> const Type& min(const Matrix<Form, Type>&);

template<class Form, class Type> Form operator-(const Matrix<Form, Type>&);

template<class Form, class Type> Form operator+
(
	const Matrix<Form, Type>&,
	const Matrix<Form, Type>&
);

template<class Form, class Type> Form operator-
(
	const Matrix<Form, Type>&,
	const Matrix<Form, Type>&
);

template<class Form, class Type> Form operator*
(
	const scalar,
	const Matrix<Form, Type>&
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#	include "MatrixI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "Matrix.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
