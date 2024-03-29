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
	Foam::nullObject

Description
	Singleton null-object class and instance

SourceFiles
	nullObjectI.H
	nullObject.C

\*---------------------------------------------------------------------------*/

#ifndef nullObject_H
#define nullObject_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class NullObject
{
	//- Private constructor
	NullObject()
	{}

	//- Prevent copy-construction
	NullObject(const NullObject&);

	//- Prevent assignment
	NullObject& operator=(const NullObject&);

public:

	//- The unique null object
	static const NullObject nullObject;
};


//- Pointer to the unique nullObject
extern const NullObject* nullObjectPtr;


//- Return reference to the nullObject of type T
template <class T>
inline const T& NullObjectRef();

//- Return pointer to the nullObject of type T
template <class T>
inline const T* NullObjectPtr();


//- Return true if t is a reference to the nullObject of type T
template <class T>
inline bool isNull(const T& t);

//- Return true if t is not a reference to the nullObject of type T
template <class T>
inline bool notNull(const T& t);


//- Return true if t is a pointer to the nullObject of type T
template <class T>
inline bool isNull(const T* t);

//- Return true if t is not a pointer to the nullObject of type T
template <class T>
inline bool notNull(const T* t);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "nullObjectI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
