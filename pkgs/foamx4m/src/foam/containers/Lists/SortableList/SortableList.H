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
	Foam::SortableList

Description
	A list that is sorted upon construction or when explicitly requested
	with the sort() method.

	Uses the Foam::stableSort() algorithm.

SourceFiles
	SortableList.C

\*---------------------------------------------------------------------------*/

#ifndef SortableList_H
#define SortableList_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class SortableList
:
	public List<T>
{
	// Private data

		//- Original indices
		labelList indices_;

		//- Resize, fill and sort the parameter according to the list values
		void sortIndices(List<label>&) const;

public:

	// Constructors

		//- Null constructor, sort later (eg, after assignment or transfer)
		SortableList();

		//- Construct from UList, sorting immediately.
		explicit SortableList(const UList<T>&);

		//- Construct from transferred List, sorting immediately.
		explicit SortableList(const Xfer<List<T> >&);

		//- Construct given size. Sort later on.
		//  The indices remain empty until the list is sorted
		explicit SortableList(const label size);

		//- Construct given size and initial value. Sort later on.
		//  The indices remain empty until the list is sorted
		SortableList(const label size, const T&);

		//- Construct as copy.
		SortableList(const SortableList<T>&);


	// Member Functions

		//- Return the list of sorted indices. Updated every sort.
		const labelList& indices() const
		{
			return indices_;
		}

		//- Return non-const access to the sorted indices. Updated every sort.
		labelList& indices()
		{
			return indices_;
		}

		//- Clear the list and the indices
		void clear();

		//- Clear the indices and return a reference to the underlying List
		List<T>& shrink();

		//- (stable) sort the list (if changed after construction time)
		//  also resizes the indices as required
		void sort();

		//- Reverse (stable) sort the list
		void reverseSort();

		//- Transfer contents to the Xfer container as a plain List
		inline Xfer<List<T> > xfer();


	// Member Operators

		//- Assignment of all entries to the given value
		inline void operator=(const T&);

		//- Assignment from UList operator. Takes linear time.
		inline void operator=(const UList<T>&);

		//- Assignment operator. Takes linear time.
		inline void operator=(const SortableList<T>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "SortableList.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
