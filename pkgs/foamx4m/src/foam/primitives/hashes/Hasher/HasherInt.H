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

InNamespace
	Foam

Description
	Optimized hashing functions.

SeeAlso
	Hasher.H

\*---------------------------------------------------------------------------*/

#ifndef HasherInt_H
#define HasherInt_H

#include "Hasher.H"

#include <cstddef>
#include <stdint.h>    // C++0x uses <cstdint>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	//- An optimized version of Hasher
	//  @param[in] data - an array of uint32_t values
	//  @param[in] len  - the number of values (not bytes)
	//  @param[in] seed - the previous hash, or an arbitrary value
	unsigned HasherInt(const uint32_t*, size_t length, unsigned seed = 0);

	//- An optimized version of Hasher, returning dual hash values
	//  @param[in] data - an array of uint32_t values
	//  @param[in] len  - the number of values (not bytes)
	//  @param[in] hash1 - the previous hash, or an arbitrary value
	//    on output, the primary hash value
	//  @param[in] hash1 - the previous hash, or an arbitrary value
	//    on output, the secondary hash value
	unsigned HasherDual
	(
		const uint32_t*,
		size_t len,
		unsigned& hash1,
		unsigned& hash2
	);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
