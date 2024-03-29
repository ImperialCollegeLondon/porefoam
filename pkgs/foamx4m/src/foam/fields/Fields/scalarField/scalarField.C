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

Description
	Specialisation of Field\<T\> for scalar.

\*---------------------------------------------------------------------------*/

#include "scalarField.H"

#define TEMPLATE
#include "FieldFunctionsM.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<>
tmp<scalarField> scalarField::component(const direction) const
{
	return *this;
}

void component(scalarField& sf, const UList<scalar>& f, const direction)
{
	sf = f;
}

template<>
void scalarField::replace(const direction, const UList<scalar>& sf)
{
	*this = sf;
}

template<>
void scalarField::replace(const direction, const scalar& s)
{
	*this = s;
}


void stabilise(scalarField& res, const UList<scalar>& sf, const scalar s)
{
	TFOR_ALL_F_OP_FUNC_S_F
	(
		scalar, res, =, ::Foam::stabilise, scalar, s, scalar, sf
	)
}

tmp<scalarField> stabilise(const UList<scalar>& sf, const scalar s)
{
	tmp<scalarField> tRes(new scalarField(sf.size()));
	stabilise(tRes(), sf, s);
	return tRes;
}

tmp<scalarField> stabilise(const tmp<scalarField>& tsf, const scalar s)
{
	tmp<scalarField> tRes = reuseTmp<scalar, scalar>::New(tsf);
	stabilise(tRes(), tsf(), s);
	reuseTmp<scalar, scalar>::clear(tsf);
	return tRes;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

BINARY_TYPE_OPERATOR(scalar, scalar, scalar, +, add)
BINARY_TYPE_OPERATOR(scalar, scalar, scalar, -, subtract)

BINARY_OPERATOR(scalar, scalar, scalar, *, multiply)
BINARY_OPERATOR(scalar, scalar, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(scalar, scalar, scalar, /, divide)

BINARY_FUNCTION(scalar, scalar, scalar, pow)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, pow)

BINARY_FUNCTION(scalar, scalar, scalar, atan2)
BINARY_TYPE_FUNCTION(scalar, scalar, scalar, atan2)

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

UNARY_FUNCTION(scalar, scalar, pow3)
UNARY_FUNCTION(scalar, scalar, pow4)
UNARY_FUNCTION(scalar, scalar, pow5)
UNARY_FUNCTION(scalar, scalar, pow6)
UNARY_FUNCTION(scalar, scalar, sqrt)
UNARY_FUNCTION(scalar, scalar, cbrt)
UNARY_FUNCTION(scalar, scalar, sign)
UNARY_FUNCTION(scalar, scalar, pos)
UNARY_FUNCTION(scalar, scalar, neg)
UNARY_FUNCTION(scalar, scalar, exp)
UNARY_FUNCTION(scalar, scalar, log)
UNARY_FUNCTION(scalar, scalar, log10)
UNARY_FUNCTION(scalar, scalar, sin)
UNARY_FUNCTION(scalar, scalar, cos)
UNARY_FUNCTION(scalar, scalar, tan)
UNARY_FUNCTION(scalar, scalar, asin)
UNARY_FUNCTION(scalar, scalar, acos)
UNARY_FUNCTION(scalar, scalar, atan)
UNARY_FUNCTION(scalar, scalar, sinh)
UNARY_FUNCTION(scalar, scalar, cosh)
UNARY_FUNCTION(scalar, scalar, tanh)
UNARY_FUNCTION(scalar, scalar, asinh)
UNARY_FUNCTION(scalar, scalar, acosh)
UNARY_FUNCTION(scalar, scalar, atanh)
UNARY_FUNCTION(scalar, scalar, erf)
UNARY_FUNCTION(scalar, scalar, erfc)
UNARY_FUNCTION(scalar, scalar, lgamma)

#ifndef __STRICT_ANSI__

UNARY_FUNCTION(scalar, scalar, j0)
UNARY_FUNCTION(scalar, scalar, j1)
UNARY_FUNCTION(scalar, scalar, y0)
UNARY_FUNCTION(scalar, scalar, y1)


#define BesselFunc(func)                                                      \
void func(scalarField& res, const int n, const UList<scalar>& sf)             \
{                                                                             \
	TFOR_ALL_F_OP_FUNC_S_F(scalar, res, =, ::Foam::func, int, n, scalar, sf)  \
}                                                                             \
					                                                          \
tmp<scalarField> func(const int n, const UList<scalar>& sf)                   \
{                                                                             \
	tmp<scalarField> tRes(new scalarField(sf.size()));                        \
	func(tRes(), n, sf);                                                      \
	return tRes;                                                              \
}                                                                             \
					                                                          \
tmp<scalarField> func(const int n, const tmp<scalarField>& tsf)               \
{                                                                             \
	tmp<scalarField> tRes = reuseTmp<scalar, scalar>::New(tsf);               \
	func(tRes(), n, tsf());                                                   \
	reuseTmp<scalar, scalar>::clear(tsf);                                     \
	return tRes;                                                              \
}

BesselFunc(jn)
BesselFunc(yn)

#undef BesselFunc

#endif //__STRICT_ANSI__

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
