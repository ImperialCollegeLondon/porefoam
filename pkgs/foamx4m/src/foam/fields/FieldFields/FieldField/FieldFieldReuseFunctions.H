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

#ifndef FieldFieldReuseFunctions_H
#define FieldFieldReuseFunctions_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<template<class> class Field, class TypeR, class Type1>
class reuseTmpFieldField
{
public:

	static tmp<FieldField<Field, TypeR> > New
	(
		const tmp<FieldField<Field, Type1> >& tf1
	)
	{
		return tmp<FieldField<Field, TypeR> >
		(
			FieldField<Field, TypeR>::NewCalculatedType(tf1())
		);
	}

	static void clear(const tmp<FieldField<Field, Type1> >& tf1)
	{
		tf1.clear();
	}
};


template<template<class> class Field, class TypeR>
class reuseTmpFieldField<Field, TypeR, TypeR>
{
public:

	static tmp<FieldField<Field, TypeR> > New
	(
		const tmp<FieldField<Field, TypeR> >& tf1
	)
	{
		if (tf1.isTmp())
		{
			return tf1;
		}
		else
		{
			return tmp<FieldField<Field, TypeR> >
			(
				FieldField<Field, TypeR>::NewCalculatedType(tf1())
			);
		}
	}

	static void clear(const tmp<FieldField<Field, TypeR> >& tf1)
	{
		if (tf1.isTmp())
		{
			tf1.ptr();
		}
	}
};


template
<
	template<class> class Field,
	class TypeR,
	class Type1,
	class Type12,
	class Type2
>
class reuseTmpTmpFieldField
{
public:

	static tmp<FieldField<Field, TypeR> > New
	(
		const tmp<FieldField<Field, Type1> >& tf1,
		const tmp<FieldField<Field, Type2> >& tf2
	)
	{
		return tmp<FieldField<Field, TypeR> >
		(
			FieldField<Field, TypeR>::NewCalculatedType(tf1())
		);
	}

	static void clear
	(
		const tmp<FieldField<Field, Type1> >& tf1,
		const tmp<FieldField<Field, Type2> >& tf2
	)
	{
		tf1.clear();
		tf2.clear();
	}
};


template<template<class> class Field, class TypeR, class Type1, class Type12>
class reuseTmpTmpFieldField<Field, TypeR, Type1, Type12, TypeR>
{
public:

	static tmp<FieldField<Field, TypeR> > New
	(
		const tmp<FieldField<Field, Type1> >& tf1,
		const tmp<FieldField<Field, TypeR> >& tf2
	)
	{
		if (tf2.isTmp())
		{
			return tf2;
		}
		else
		{
			return tmp<FieldField<Field, TypeR> >
			(
				FieldField<Field, TypeR>::NewCalculatedType(tf1())
			);
		}
	}

	static void clear
	(
		const tmp<FieldField<Field, Type1> >& tf1,
		const tmp<FieldField<Field, TypeR> >& tf2
	)
	{
		tf1.clear();
		if (tf2.isTmp())
		{
			tf2.ptr();
		}
	}
};


template<template<class> class Field, class TypeR, class Type2>
class reuseTmpTmpFieldField<Field, TypeR, TypeR, TypeR, Type2>
{
public:

	static tmp<FieldField<Field, TypeR> > New
	(
		const tmp<FieldField<Field, TypeR> >& tf1,
		const tmp<FieldField<Field, Type2> >& tf2
	)
	{
		if (tf1.isTmp())
		{
			return tf1;
		}
		else
		{
			return tmp<FieldField<Field, TypeR> >
			(
				FieldField<Field, TypeR>::NewCalculatedType(tf1())
			);
		}
	}

	static void clear
	(
		const tmp<FieldField<Field, TypeR> >& tf1,
		const tmp<FieldField<Field, Type2> >& tf2
	)
	{
		if (tf1.isTmp())
		{
			tf1.ptr();
		}
		tf2.clear();
	}
};


template<template<class> class Field, class TypeR>
class reuseTmpTmpFieldField<Field, TypeR, TypeR, TypeR, TypeR>
{
public:

	static tmp<FieldField<Field, TypeR> > New
	(
		const tmp<FieldField<Field, TypeR> >& tf1,
		const tmp<FieldField<Field, TypeR> >& tf2
	)
	{
		if (tf1.isTmp())
		{
			return tf1;
		}
		else if (tf2.isTmp())
		{
			return tf2;
		}
		else
		{
			return tmp<FieldField<Field, TypeR> >
			(
				FieldField<Field, TypeR>::NewCalculatedType(tf1())
			);
		}
	}

	static void clear
	(
		const tmp<FieldField<Field, TypeR> >& tf1,
		const tmp<FieldField<Field, TypeR> >& tf2
	)
	{
		if (tf1.isTmp())
		{
			tf1.ptr();
			tf2.clear();
		}
		else if (tf2.isTmp())
		{
			tf1.clear();
			tf2.ptr();
		}
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
