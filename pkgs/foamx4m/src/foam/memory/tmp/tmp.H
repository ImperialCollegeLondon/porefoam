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
	Foam::tmp

Description
	A class for managing temporary objects

SourceFiles
	tmpI.H

\*---------------------------------------------------------------------------*/

#ifndef tmp_H
#define tmp_H

#include "refCount.H"
#include <cstddef>

#if defined(__GNUC__) && !defined(__INTEL_COMPILER) && !defined(darwin) && !defined(ARM_CLANG)
#	define ConstructFromTmp
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class tmp
{
	// Private data

		//- Flag for whether object is a temporary or a constant object
		bool isTmp_;

		//- Pointer to temporary object
		mutable T* ptr_;

		//- Const reference to constant object
		const T& ref_;


public:

	// Constructors

		//- Store object pointer
		inline explicit tmp(T* = 0);

		//- Store object const reference
		inline tmp(const T&);

		//- Construct copy and increment reference count
		inline tmp(const tmp<T>&);

		//- Construct copy transferring content of temporary if required
		inline tmp(const tmp<T>&, bool allowTransfer);


	//- Destructor, delete object when reference count == 0
	inline ~tmp();


	// Member Functions

		// Access

			//- Return true if this is really a temporary object
			inline bool isTmp() const;

			//- Return true if this temporary object empty,
			//  ie, a temporary without allocation
			inline bool empty() const;

			//- Is this temporary object valid,
			//  ie, it is a reference or a temporary that has been allocated
			inline bool valid() const;

		// Edit

			//- Return tmp pointer for reuse
			inline T* ptr() const;

			//- Return non-const reference. Fatal if the object is const.
			inline T& ref() const;

			//- If object pointer points to valid object:
			//  delete object and set pointer to nullptr
			inline void clear() const;


	// Member operators

		//- Dereference operator
		inline T& operator()();

		//- Const dereference operator
		inline const T& operator()() const;

		//- Const cast to the underlying type reference
		inline operator const T&() const;

		//- Return object pointer
		inline T* operator->();

		//- Return const object pointer
		inline const T* operator->() const;

		//- Assignment operator
		inline void operator=(const tmp<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "tmpI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
