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

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList()
:
	List<T>(SizeInc),
	capacity_(SizeInc)
{
	List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
	const label nElem
)
:
	List<T>(nElem),
	capacity_(nElem)
{
	// we could also enforce SizeInc granularity when (!SizeMult || !SizeDiv)
	List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
	const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
:
	List<T>(lst),
	capacity_(lst.size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
	const UList<T>& lst
)
:
	List<T>(lst),
	capacity_(lst.size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
	const UIndirectList<T>& lst
)
:
	List<T>(lst),
	capacity_(lst.size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::DynamicList
(
	const Xfer<List<T> >& lst
)
:
	List<T>(lst),
	capacity_(List<T>::size())
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::label Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::capacity()
const
{
	return capacity_;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setCapacity
(
	const label nElem
)
{
	label nextFree = List<T>::size();
	capacity_ = nElem;

	if (nextFree > capacity_)
	{
		// truncate addressed sizes too
		nextFree = capacity_;
	}
	// we could also enforce SizeInc granularity when (!SizeMult || !SizeDiv)

	List<T>::setSize(capacity_);
	List<T>::size(nextFree);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::reserve
(
	const label nElem
)
{
	// allocate more capacity?
	if (nElem > capacity_)
	{
// TODO: convince the compiler that division by zero does not occur
//        if (SizeInc && (!SizeMult || !SizeDiv))
//        {
//            // resize with SizeInc as the granularity
//            capacity_ = nElem;
//            unsigned pad = SizeInc - (capacity_ % SizeInc);
//            if (pad != SizeInc)
//            {
//                capacity_ += pad;
//            }
//        }
//        else
		{
			capacity_ = max
			(
				nElem,
				label(SizeInc + capacity_ * SizeMult / SizeDiv)
			);
		}

		// adjust allocated size, leave addressed size untouched
		label nextFree = List<T>::size();
		List<T>::setSize(capacity_);
		List<T>::size(nextFree);
	}
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
	const label nElem
)
{
	// allocate more capacity?
	if (nElem > capacity_)
	{
// TODO: convince the compiler that division by zero does not occur
//        if (SizeInc && (!SizeMult || !SizeDiv))
//        {
//            // resize with SizeInc as the granularity
//            capacity_ = nElem;
//            unsigned pad = SizeInc - (capacity_ % SizeInc);
//            if (pad != SizeInc)
//            {
//                capacity_ += pad;
//            }
//        }
//        else
		{
			capacity_ = max
			(
				nElem,
				label(SizeInc + capacity_ * SizeMult / SizeDiv)
			);
		}

		List<T>::setSize(capacity_);
	}

	// adjust addressed size
	List<T>::size(nElem);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::setSize
(
	const label nElem,
	const T& t
)
{
	label nextFree = List<T>::size();
	setSize(nElem);

	// set new elements to constant value
	while (nextFree < nElem)
	{
		this->operator[](nextFree++) = t;
	}
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::resize
(
	const label nElem
)
{
	this->setSize(nElem);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::resize
(
	const label nElem,
	const T& t
)
{
	this->setSize(nElem, t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clear()
{
	List<T>::size(0);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::clearStorage()
{
	List<T>::clear();
	capacity_ = 0;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>&
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::shrink()
{
	label nextFree = List<T>::size();
	if (capacity_ > nextFree)
	{
		// use the full list when resizing
		List<T>::size(capacity_);

		// the new size
		capacity_ = nextFree;
		List<T>::setSize(capacity_);
		List<T>::size(nextFree);
	}
	return *this;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer(List<T>& lst)
{
	capacity_ = lst.size();
	List<T>::transfer(lst);   // take over storage, clear addressing for lst.
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::transfer
(
	DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
	// take over storage as-is (without shrink), clear addressing for lst.
	capacity_ = lst.capacity_;
	lst.capacity_ = 0;

	List<T>::transfer(static_cast<List<T>&>(lst));
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline Foam::Xfer< Foam::List<T> >
Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::xfer()
{
	return xferMoveTo< List<T> >(*this);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
	const T& t
)
{
	label elemI = List<T>::size();
	setSize(elemI + 1);

	this->operator[](elemI) = t;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
	const UList<T>& lst
)
{
	if (this == &lst)
	{
		FatalErrorIn
		(
			"DynamicList<T, SizeInc, SizeMult, SizeDiv>::append"
			"(const UList<T>&)"
		)   << "attempted appending to self" << abort(FatalError);
	}

	label nextFree = List<T>::size();
	setSize(nextFree + lst.size());

	forAll(lst, elemI)
	{
		this->operator[](nextFree++) = lst[elemI];
	}
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::append
(
	const UIndirectList<T>& lst
)
{
	label nextFree = List<T>::size();
	setSize(nextFree + lst.size());

	forAll(lst, elemI)
	{
		this->operator[](nextFree++) = lst[elemI];
	}
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()
{
	label elemI = List<T>::size() - 1;

	if (elemI < 0)
	{
		FatalErrorIn
		(
			"Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::remove()"
		)   << "List is empty" << abort(FatalError);
	}

	const T& val = List<T>::operator[](elemI);

	List<T>::size(elemI);

	return val;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline T& Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator()
(
	const label elemI
)
{
	if (elemI >= List<T>::size())
	{
		setSize(elemI + 1);
	}

	return this->operator[](elemI);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
	const T& t
)
{
	UList<T>::operator=(t);
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
	const UList<T>& lst
)
{
	if (capacity_ >= lst.size())
	{
		// can copy w/o reallocating, match initial size to avoid reallocation
		List<T>::size(lst.size());
		List<T>::operator=(lst);
	}
	else
	{
		// make everything available for the copy operation
		List<T>::size(capacity_);

		List<T>::operator=(lst);
		capacity_ = List<T>::size();
	}
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
inline void Foam::DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator=
(
	const DynamicList<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
	if (this == &lst)
	{
		FatalErrorIn
		(
			"DynamicList<T, SizeInc, SizeMult, SizeDiv>::operator="
			"(const DynamicList<T, SizeInc, SizeMult, SizeDiv>&)"
		)   << "attempted assignment to self" << abort(FatalError);
	}

	if (capacity_ >= lst.size())
	{
		// can copy w/o reallocating, match initial size to avoid reallocation
		List<T>::size(lst.size());
		List<T>::operator=(lst);
	}
	else
	{
		// make everything available for the copy operation
		List<T>::size(capacity_);

		List<T>::operator=(lst);
		capacity_ = List<T>::size();
	}
}


// ************************************************************************* //
