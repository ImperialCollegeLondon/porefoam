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

Typedef
	Foam::uLabel

Description
	A uLabel is an uint32_t or uint64_t as specified by the pre-processor macro
	WM_LABEL_SIZE.

	A readLabel function is defined so that uLabel can be constructed from
	Istream.

\*---------------------------------------------------------------------------*/

#ifndef uLabel_H
#define uLabel_H

#include "uint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define UINT_ADD_SIZE(x,s,y) x ## s ## y
#define UINT_ADD_DEF_SIZE(x,s,y) UINT_ADD_SIZE(x,s,y)
#define UINT_SIZE(x,y) UINT_ADD_DEF_SIZE(x,WM_LABEL_SIZE,y)

#if WM_LABEL_SIZE != 32 && WM_LABEL_SIZE != 64
	#error "uLabel.H: WM_LABEL_SIZE must be set to either 32 or 64"
#endif


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

typedef UINT_SIZE(uint, _t) uLabel;

static const uLabel uLabelMin = UINT_SIZE(UINT, _MIN);
static const uLabel uLabelMax = UINT_SIZE(UINT, _MAX);

inline uLabel readULabel(Istream& is)
{
	return UINT_SIZE(readUint,) (is);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//- Raise one uLabel to the power of another
uLabel pow(uLabel a, uLabel b);

//- Evaluate n! : 0 < n <= 12
uLabel factorial(uLabel n);


inline uLabel& setComponent(uLabel& l, const direction)
{
	return l;
}

inline uLabel component(const uLabel l, const direction)
{
	return l;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "uLabelSpecific.H"

#undef UINT_ADD_SIZE
#undef UINT_ADD_DEF_SIZE
#undef UINT_SIZE

#endif

// ************************************************************************* //
