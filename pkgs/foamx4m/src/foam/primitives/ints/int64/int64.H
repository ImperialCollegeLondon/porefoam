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

Primitive
	int64

Description
	64bit integer

SourceFiles
	int64.C
	int64IO.C

\*---------------------------------------------------------------------------*/

#ifndef int64_H
#define int64_H

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <climits>
#include <cstdlib>

#include "word.H"
#include "pTraits.H"
#include "direction.H"

#ifndef UINT64_MIN
#define UINT64_MIN 0
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a word representation of an int64
word name(const int64_t);

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

int64_t readInt64(Istream&);
bool read(const char*, int64_t&);
Istream& operator>>(Istream&, int64_t&);
Ostream& operator<<(Ostream&, const int64_t);

#if WM_ARCH_OPTION == 64 && darwin && __clang__
	Istream& operator>>(Istream&, long&);
	Ostream& operator<<(Ostream&, const long);
#endif

#if defined(mingw)
	Istream& operator>>(Istream&, off_t&);
	Ostream& operator<<(Ostream&, const off_t);
#endif

//- Template specialization for pTraits<int64_t>
template<>
class pTraits<int64_t>
{
	int64_t p_;

public:

	//- Component type
	typedef int64_t cmptType;

	// Member constants

		enum
		{
			dim = 3,         //!< Dimensionality of space
			rank = 0,        //!< Rank of int64_t is 0
			nComponents = 1  //!< Number of components in int64_t is 1
		};


	// Static data members

		static const char* const typeName;
		static const char* componentNames[];
		static const int64_t zero;
		static const int64_t one;
		static const int64_t min;
		static const int64_t max;
		static const int64_t rootMax;
		static const int64_t rootMin;


	// Constructors

		//- Construct from primitive
		explicit pTraits(const int64_t&);

		//- Construct from Istream
		pTraits(Istream&);


	// Member Functions

		//- Access to the int64_t value
		operator int64_t() const
		{
			return p_;
		}

		//- Access to the int value
		operator int64_t&()
		{
			return p_;
		}
};


inline int64_t mag(const int64_t l)
{
	return ::labs(l);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
