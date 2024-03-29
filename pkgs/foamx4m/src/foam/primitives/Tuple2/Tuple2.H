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
	Foam::Tuple2

Description
	A 2-tuple: storing two objects of different types.

SeeAlso
	Foam::Pair for storing two objects of identical types.

\*---------------------------------------------------------------------------*/

#ifndef Tuple2_H
#define Tuple2_H

#include "Istream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Type1, class Type2>
class Tuple2;

template<class Type1, class Type2>
inline bool operator==
(
	const Tuple2<Type1, Type2>&,
	const Tuple2<Type1, Type2>&
);

template<class Type1, class Type2>
inline bool operator!=
(
	const Tuple2<Type1, Type2>&,
	const Tuple2<Type1, Type2>&
);

template<class Type1, class Type2>
inline Istream& operator>>(Istream&, Tuple2<Type1, Type2>&);

template<class Type1, class Type2>
inline Ostream& operator<<(Ostream&, const Tuple2<Type1, Type2>&);



template<class Type1, class Type2>
class Tuple2
{
	// Private data

		Type1 f_;
		Type2 s_;


public:

	// Constructors

		//- Null constructor for lists
		inline Tuple2()
		{}

		//- Construct from components
		inline Tuple2(const Type1& f, const Type2& s)
		:
			f_(f),
			s_(s)
		{}

		//- Construct from Istream
		inline Tuple2(Istream& is)
		{
			is >> *this;
		}


	// Member Functions

		//- Return first
		inline const Type1& first() const
		{
			return f_;
		}

		//- Return first
		inline Type1& first()
		{
			return f_;
		}

		//- Return second
		inline const Type2& second() const
		{
			return s_;
		}

		//- Return second
		inline Type2& second()
		{
			return s_;
		}

		//- Return reverse pair
		inline Tuple2<Type2, Type1> reverseTuple2() const
		{
			return Tuple2<Type2, Type1>(second(), first());
		}


	// Friend Operators

		friend bool operator== <Type1, Type2>
		(
			const Tuple2<Type1, Type2>& a,
			const Tuple2<Type1, Type2>& b
		);

		friend bool operator!= <Type1, Type2>
		(
			const Tuple2<Type1, Type2>& a,
			const Tuple2<Type1, Type2>& b
		);


	// IOstream operators

		//- Read Tuple2 from Istream, discarding contents of existing Tuple2.
		friend Istream& operator>> <Type1, Type2>
		(
			Istream& is,
			Tuple2<Type1, Type2>& t2
		);

		// Write Tuple2 to Ostream.
		friend Ostream& operator<< <Type1, Type2>
		(
			Ostream& os,
			const Tuple2<Type1, Type2>& t2
		);
};


template<class Type1, class Type2>
inline bool operator==
(
	const Tuple2<Type1, Type2>& a,
	const Tuple2<Type1, Type2>& b
)
{
	return
	(
		(a.first() == b.first()) && (a.second() == b.second())
	);
}


template<class Type1, class Type2>
inline bool operator!=
(
	const Tuple2<Type1, Type2>& a,
	const Tuple2<Type1, Type2>& b
)
{
	return !(a == b);
}


template<class Type1, class Type2>
inline Istream& operator>>(Istream& is, Tuple2<Type1, Type2>& t2)
{
	is.readBegin("Tuple2");
	is >> t2.f_ >> t2.s_;
	is.readEnd("Tuple2");

	// Check state of Istream
	is.check("operator>>(Istream&, Tuple2<Type1, Type2>&)");

	return is;
}


template<class Type1, class Type2>
inline Ostream& operator<<(Ostream& os, const Tuple2<Type1, Type2>& t2)
{
	os  << token::BEGIN_LIST
		<< t2.f_ << token::SPACE << t2.s_
		<< token::END_LIST;

	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
