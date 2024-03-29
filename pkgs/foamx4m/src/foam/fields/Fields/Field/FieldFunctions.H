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

\*---------------------------------------------------------------------------*/

#include "scalarField.H"

#define TEMPLATE template<class Type>
#include "FieldFunctionsM.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Type>
void component
(
	Field<typename Field<Type>::cmptType>& res,
	const UList<Type>& f,
	const direction d
);


template<class Type>
void T(Field<Type>& res, const UList<Type>& f);


template<class Type, int r>
void pow
(
	Field<typename powProduct<Type, r>::type>& res,
	const UList<Type>& vf
);


template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
	const UList<Type>& f,
	typename powProduct<Type, r>::type
	  = pTraits<typename powProduct<Type, r>::type>::zero
);

template<class Type, int r>
tmp<Field<typename powProduct<Type, r>::type> >
pow
(
	const tmp<Field<Type> >& tf,
	typename powProduct<Type, r>::type
	  = pTraits<typename powProduct<Type, r>::type>::zero
);


template<class Type>
void sqr
(
	Field<typename outerProduct<Type, Type>::type>& res,
	const UList<Type>& vf
);

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const UList<Type>& f);

template<class Type>
tmp<Field<typename outerProduct<Type, Type>::type> >
sqr(const tmp<Field<Type> >& tf);


template<class Type>
void magSqr(Field<scalar>& res, const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > magSqr(const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > magSqr(const tmp<Field<Type> >& tf);


template<class Type>
void mag(Field<scalar>& res, const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > mag(const UList<Type>& f);

template<class Type>
tmp<Field<scalar> > mag(const tmp<Field<Type> >& tf);


template<class Type>
void cmptMax(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptMax(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> >
cmptMax(const tmp<Field<Type> >& tf);


template<class Type>
void cmptMin(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptMin(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> >
cmptMin(const tmp<Field<Type> >& tf);


template<class Type>
void cmptAv(Field<typename Field<Type>::cmptType>& res, const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const UList<Type>& f);

template<class Type>
tmp<Field<typename Field<Type>::cmptType> > cmptAv(const tmp<Field<Type> >& tf);


template<class Type>
void cmptMag(Field<Type>& res, const UList<Type>& f);

template<class Type>
tmp<Field<Type> > cmptMag(const UList<Type>& f);

template<class Type>
tmp<Field<Type> > cmptMag(const tmp<Field<Type> >& tf);


template<class Type>
void cmptSign(Field<Type>& res, const UList<Type>& f);

template<class Type>
tmp<Field<Type> > cmptSign(const UList<Type>& f);

template<class Type>
tmp<Field<Type> > cmptSign(const tmp<Field<Type> >& tf);

#define TMP_UNARY_FUNCTION(ReturnType, Func)                                  \
					                                                          \
template<class Type>                                                          \
ReturnType Func(const tmp<Field<Type> >& tf1);

template<class Type>
Type max(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, max)

template<class Type>
Type min(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, min)

template<class Type>
Type sum(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, sum)


template<class Type>
scalar sumProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
Type sumCmptProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
scalar sumSqr(const UList<Type>& f);

TMP_UNARY_FUNCTION(scalar, sumSqr)

template<class Type>
scalar sumMag(const UList<Type>& f);

TMP_UNARY_FUNCTION(scalar, sumMag)

template<class Type>
Type sumCmptMag(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, sumCmptMag)

template<class Type>
Type average(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, average)


#define G_UNARY_FUNCTION(ReturnType, gFunc, Func, rFunc)                      \
					                                                          \
template<class Type>                                                          \
ReturnType gFunc(const UList<Type>& f);                                       \
TMP_UNARY_FUNCTION(ReturnType, gFunc)

G_UNARY_FUNCTION(Type, gMax, max, max)
G_UNARY_FUNCTION(Type, gMin, min, min)
G_UNARY_FUNCTION(Type, gSum, sum, sum)
G_UNARY_FUNCTION(scalar, gSumSqr, sumSqr, sum)
G_UNARY_FUNCTION(scalar, gSumMag, sumMag, sum)
G_UNARY_FUNCTION(Type, gSumCmptMag, sumCmptMag, sum)

#undef G_UNARY_FUNCTION

template<class Type>
scalar gSumProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
Type gSumCmptProd(const UList<Type>& f1, const UList<Type>& f2);

template<class Type>
Type gAverage(const UList<Type>& f);

TMP_UNARY_FUNCTION(Type, gAverage)

#undef TMP_UNARY_FUNCTION


BINARY_FUNCTION(Type, Type, Type, max)
BINARY_FUNCTION(Type, Type, Type, min)
BINARY_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_FUNCTION(Type, Type, Type, cmptDivide)

BINARY_TYPE_FUNCTION(Type, Type, Type, max)
BINARY_TYPE_FUNCTION(Type, Type, Type, min)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptMultiply)
BINARY_TYPE_FUNCTION(Type, Type, Type, cmptDivide)


// * * * * * * * * * * * * * * * Global operators  * * * * * * * * * * * * * //

UNARY_OPERATOR(Type, Type, -, negate)

BINARY_OPERATOR(Type, Type, scalar, *, multiply)
BINARY_OPERATOR(Type, scalar, Type, *, multiply)
BINARY_OPERATOR(Type, Type, scalar, /, divide)

BINARY_TYPE_OPERATOR_SF(Type, scalar, Type, *, multiply)
BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, *, multiply)

BINARY_TYPE_OPERATOR_FS(Type, Type, scalar, /, divide)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define PRODUCT_OPERATOR(product, Op, OpFunc)                                 \
					                                                          \
template<class Type1, class Type2>                                            \
void OpFunc                                                                   \
(                                                                             \
	Field<typename product<Type1, Type2>::type>& res,                         \
	const UList<Type1>& f1,                                                   \
	const UList<Type2>& f2                                                    \
);                                                                            \
					                                                          \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const UList<Type1>& f1, const UList<Type2>& f2);                  \
					                                                          \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const UList<Type1>& f1, const tmp<Field<Type2> >& tf2);           \
					                                                          \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<Field<Type1> >& tf1, const UList<Type2>& f2);           \
					                                                          \
template<class Type1, class Type2>                                            \
tmp<Field<typename product<Type1, Type2>::type> >                             \
operator Op(const tmp<Field<Type1> >& tf1, const tmp<Field<Type2> >& tf2);    \
					                                                          \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
void OpFunc                                                                   \
(                                                                             \
	Field<typename product<Type, Form>::type>& res,                           \
	const UList<Type>& f1,                                                    \
	const VectorSpace<Form,Cmpt,nCmpt>& vs                                    \
);                                                                            \
					                                                          \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator Op(const UList<Type>& f1, const VectorSpace<Form,Cmpt,nCmpt>& vs);   \
					                                                          \
template<class Type, class Form, class Cmpt, int nCmpt>                       \
tmp<Field<typename product<Type, Form>::type> >                               \
operator Op(const tmp<Field<Type> >&tf1,const VectorSpace<Form,Cmpt,nCmpt>&vs);\
					                                                          \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
void OpFunc                                                                   \
(                                                                             \
	Field<typename product<Form, Type>::type>& res,                           \
	const VectorSpace<Form,Cmpt,nCmpt>& vs,                                   \
	const UList<Type>& f1                                                     \
);                                                                            \
					                                                          \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>& vs, const UList<Type>& f1);   \
					                                                          \
template<class Form, class Cmpt, int nCmpt, class Type>                       \
tmp<Field<typename product<Form, Type>::type> >                               \
operator Op(const VectorSpace<Form,Cmpt,nCmpt>&vs,const tmp<Field<Type> >&tf1);

PRODUCT_OPERATOR(typeOfSum, +, add)
PRODUCT_OPERATOR(typeOfSum, -, subtract)

PRODUCT_OPERATOR(outerProduct, *, outer)
PRODUCT_OPERATOR(crossProduct, ^, cross)
PRODUCT_OPERATOR(innerProduct, &, dot)
PRODUCT_OPERATOR(scalarProduct, &&, dotdot)

#undef PRODUCT_OPERATOR


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "undefFieldFunctionsM.H"

// ************************************************************************* //
