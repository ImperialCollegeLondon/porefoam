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
	uint32

Description
	32bit uinteger

SourceFiles
	uint32.C
	uint32IO.C

\*---------------------------------------------------------------------------*/

#ifndef uint32_H
#define uint32_H

#define __STDC_LIMIT_MACROS
#include <stdint.h>
#include <climits>
#include <cstdlib>

#include "word.H"
#include "pTraits.H"
#include "direction.H"

#ifndef UINT32_MIN
#define UINT32_MIN 0
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class Istream;
class Ostream;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Return a word representation of an uint32
word name(const uint32_t);

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

uint32_t readUint32(Istream&);
bool read(const char*, uint32_t&);
Istream& operator>>(Istream&, uint32_t&);
Ostream& operator<<(Ostream&, const uint32_t);

//- Template specialization for pTraits<uint32_t>
template<>
class pTraits<uint32_t>
{
	uint32_t p_;

public:

	//- Component type
	typedef uint32_t cmptType;

	// Member constants

		enum
		{
			dim = 3,         //!< Dimensionality of space
			rank = 0,        //!< Rank of uint32_t is 0
			nComponents = 1  //!< Number of components in uint32_t is 1
		};


	// Static data members

		static const char* const typeName;
		static const char* componentNames[];
		static const uint32_t zero;
		static const uint32_t one;
		static const uint32_t min;
		static const uint32_t max;
		static const uint32_t rootMax;
		static const uint32_t rootMin;


	// Constructors

		//- Construct from primitive
		explicit pTraits(const uint32_t&);

		//- Construct from Istream
		pTraits(Istream&);


	// Member Functions

		//- Access to the uint32_t value
		operator uint32_t() const
		{
			return p_;
		}

		//- Access to the uint32_t value
		operator uint32_t&()
		{
			return p_;
		}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
