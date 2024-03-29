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
	Foam::SubList

Description
	A List obtained as a section of another List.

	Since the SubList is itself unallocated, no storage is allocated or
	de-allocated during its use.  To achieve this behaviour, SubList is
	derived from UList rather than List.

SourceFiles
	SubListI.H

\*---------------------------------------------------------------------------*/

#ifndef SubList_H
#define SubList_H

#include "List.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class SubList
:
	public UList<T>
{

public:

	// Static Member Functions

		//- Return a null SubList
		//  Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
		inline static const SubList<T>& null();


	// Constructors

		//- Construct from UList and sub-list size
		inline SubList
		(
			const UList<T>& list,
			const label subSize
		);

		//- Construct from UList, sub-list size and start index
		inline SubList
		(
			const UList<T>& list,
			const label subSize,
			const label startIndex
		);


	// Member functions

		//- Reset the SubList
		//- Construct from UList and SubList size
		inline void reset
		(
			const UList<T>& list,
			const label subSize
		);

		//- Construct from UList, start and end indices
		inline void reset
		(
			const UList<T>& list,
			const label subSize,
			const label startIndex
		);


	// Member operators

		//- Allow cast to a const List<T>&
		inline operator const Foam::List<T>&() const;

		//- Assignment of all entries to the given value
		inline void operator=(const T&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "SubListI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
