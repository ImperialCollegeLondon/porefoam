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
	uint64

Description
	64bit uinteger

SourceFiles
	uint64.C
	uint64IO.C

\*---------------------------------------------------------------------------*/

#ifndef uint64_H
#define uint64_H

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

//- Return a word representation of an uint64
word name(const uint64_t);

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint64_t readUint64(Istream&);
bool read(const char*, uint64_t&);
Istream& operator>>(Istream&, uint64_t&);
Ostream& operator<<(Ostream&, const uint64_t);

#if WM_ARCH_OPTION == 64 && darwin && __clang__
Ostream& operator<<(Ostream&, const uintptr_t);
#endif

//- Template specialization for pTraits<uint64_t>
template<>
class pTraits<uint64_t>
{
	uint64_t p_;

public:

	//- Component type
	typedef uint64_t cmptType;

	// Member constants

		enum
		{
			dim = 3,         //!< Dimensionality of space
			rank = 0,        //!< Rank of uint64_t is 0
			nComponents = 1  //!< Number of components in uint64_t is 1
		};


	// Static data members

		static const char* const typeName;
		static const char* componentNames[];
		static const uint64_t zero;
		static const uint64_t one;
		static const uint64_t min;
		static const uint64_t max;
		static const uint64_t rootMax;
		static const uint64_t rootMin;


	// Constructors

		//- Construct from primitive
		explicit pTraits(const uint64_t&);

		//- Construct from Istream
		pTraits(Istream&);


	// Member Functions

		//- Access to the uint64_t value
		operator uint64_t() const
		{
			return p_;
		}

		//- Access to the uint64_t value
		operator uint64_t&()
		{
			return p_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
