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
	Foam::dimensioned

Description
	Generic dimensioned Type class

SourceFiles
	dimensionedType.C

\*---------------------------------------------------------------------------*/

#ifndef dimensionedType_H
#define dimensionedType_H

#include "word.H"
#include "direction.H"
#include "dimensionSet.H"
#include "VectorSpace.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class dictionary;

// Forward declaration of friend functions and operators
template<class Type> class dimensioned;

template<class Type>
Istream& operator>>(Istream&, dimensioned<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const dimensioned<Type>&);



template <class Type>
class dimensioned
{
	// private data

		//- Variable name
		word name_;

		//- The dimension set
		dimensionSet dimensions_;

		//- The data value
		Type value_;


public:

	//- Component type
	typedef typename pTraits<Type>::cmptType cmptType;


	// Constructors

		//- Construct given a name, a value and its dimensionSet.
		dimensioned(const word&, const dimensionSet&, const Type);

		//- Construct from a dimensioned<Type> changing the name.
		dimensioned(const word&, const dimensioned<Type>&);

		//- Construct given a value (creates dimensionless value).
		dimensioned(const Type& t)
		:
			name_(::Foam::name(t)),
			dimensions_(dimless),
			value_(t)
		{}

		//- Construct from Istream.
		explicit dimensioned(Istream&);

		//- Construct from an Istream with a given name
		dimensioned(const word&, Istream&);

		//- Construct from an Istream with a given name and dimensions
		dimensioned(const word&, const dimensionSet&, Istream&);

		//- Construct from dictionary, with default value.
		static dimensioned<Type> lookupOrDefault
		(
			const word&,
			const dictionary&,
			const Type& defaultValue = pTraits<Type>::zero,
			const dimensionSet& dims = dimless
		);

		//- Construct from dictionary, with default value.
		//  If the value is not found, it is added into the dictionary.
		static dimensioned<Type> lookupOrAddToDict
		(
			const word&,
			dictionary&,
			const Type& defaultValue = pTraits<Type>::zero,
			const dimensionSet& dims = dimless
		);


	// Member functions

		//- Return const reference to name.
		inline const word& name() const
		{
			return name_;
		}

		//- Return non-const reference to name.
		inline word& name()
		{
			return name_;
		}


		//- Return const reference to dimensions.
		const dimensionSet& dimensions() const
		{
			return dimensions_;
		}

		//- Return non-const reference to dimensions.
		inline dimensionSet& dimensions()
		{
			return dimensions_;
		}

		//- Return const reference to value.
		inline const Type& value() const
		{
			return value_;
		}

		//- Return non-const reference to value.
		inline Type& value()
		{
			return value_;
		}

		//- Return a component as a dimensioned<cmptType>
		dimensioned<cmptType> component(const direction) const;

		//- Return a component with a dimensioned<cmptType>
		void replace(const direction, const dimensioned<cmptType>&);

		//- Return transpose.
		dimensioned<Type> T() const;

		//- Update the value of dimensioned<Type> if found in the dictionary.
		bool readIfPresent(const dictionary&);


	// Member operators

		//- Return a component as a dimensioned<cmptType>
		dimensioned<cmptType> operator[](const direction) const;

		void operator+=(const dimensioned<Type>&);
		void operator-=(const dimensioned<Type>&);
		void operator*=(const scalar);
		void operator/=(const scalar);


	// IOstream operators

		friend Istream& operator>> <Type>
		(Istream&, dimensioned<Type>&);

		friend Ostream& operator<< <Type>
		(Ostream&, const dimensioned<Type>&);
};


// * * * * * * * * * * * * * * * Global Operators  * * * * * * * * * * * * * //

template<class Type, int r>
dimensioned<typename powProduct<Type, r>::type>
pow
(
	const dimensioned<Type>&,
	typename powProduct<Type, r>::type
  = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type>
dimensioned<typename outerProduct<Type, Type>::type>
sqr(const dimensioned<Type>&);

template<class Type>
dimensioned<scalar> magSqr(const dimensioned<Type>&);

template<class Type>
dimensioned<scalar> mag(const dimensioned<Type>&);

template<class Type>
dimensioned<Type> cmptMultiply
(
	const dimensioned<Type>&,
	const dimensioned<Type>&
);

template<class Type>
dimensioned<Type> cmptDivide
(
	const dimensioned<Type>&,
	const dimensioned<Type>&
);

template<class Type>
dimensioned<Type> max(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> min(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
bool operator>(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
bool operator<(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator+(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator-(const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator-(const dimensioned<Type>&, const dimensioned<Type>&);

template<class Type>
dimensioned<Type> operator*
(
	const dimensioned<scalar>&,
	const dimensioned<Type>&
);

template<class Type>
dimensioned<Type> operator/
(
	const dimensioned<Type>&,
	const dimensioned<scalar>&
);


// Products
// ~~~~~~~~

#define PRODUCT_OPERATOR(product, op, opFunc)                                 \
					                                                          \
template<class Type1, class Type2>                                            \
dimensioned<typename product<Type1, Type2>::type>                             \
operator op(const dimensioned<Type1>&, const dimensioned<Type2>&);            \
					                                                          \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
dimensioned<typename product<Type, Form>::type>                               \
operator op                                                                   \
(                                                                             \
	const dimensioned<Type>&,                                                 \
	const VectorSpace<Form,Cmpt,nCmpt>&                                       \
);                                                                            \
					                                                          \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
dimensioned<typename product<Form, Type>::type>                               \
operator op                                                                   \
(                                                                             \
	const VectorSpace<Form,Cmpt,nCmpt>&,                                      \
	const dimensioned<Type>&                                                  \
);

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "dimensionedType.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
