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
	Foam::autoPtr

Description
	An auto-pointer similar to the STL auto_ptr but with automatic casting
	to a reference to the type and with pointer allocation checking on access.

SourceFiles
	autoPtrI.H

\*---------------------------------------------------------------------------*/

#ifndef autoPtr_H
#define autoPtr_H

#include <cstddef>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class autoPtr
{
	// Public data

		//- Pointer to object
		mutable T* ptr_;


public:

	typedef T Type;


	// Constructors

		//- Store object pointer
		inline explicit autoPtr(T* = 0);

		//- Construct as copy by transfering pointer to this autoPtr and
		//  setting the arguments pointer to nullptr
		inline autoPtr(const autoPtr<T>&);

		//- Construct either by transferring pointer or cloning. Should
		//  only be called with type that supports cloning.
		inline autoPtr(const autoPtr<T>&, const bool reuse);

	// Destructor

		//- Delete object if pointer is not nullptr
		inline ~autoPtr();


	// Member Functions

		// Check

			//- Return true if the autoPtr is empty (ie, no pointer set).
			inline bool empty() const;

			//- Return true if the autoPtr valid (ie, the pointer is set).
			inline bool valid() const;


		// Edit

			//- Return object pointer for reuse
			inline T* ptr();

			//- Set pointer to that given.
			//  If object pointer already set issue a FatalError.
			inline void set(T*);

			//- If object pointer already set, delete object and set to
			//  given pointer
			inline void reset(T* = 0);

			//- If object pointer points to valid object:
			//  delete object and set pointer to nullptr
			inline void clear();


	// Member operators

		//- Return reference to the object data
		inline T& operator()();

		//- Return const reference to the object data
		inline const T& operator()() const;

		// inline T& operator*();
		// inline const T& operator*() const;

		inline operator const T&() const;

		//- Return object pointer
		inline T* operator->();

		//- Return const object pointer
		inline const T* operator->() const;

		//- Take over object pointer from parameter
		inline void operator=(const autoPtr<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "autoPtrI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
