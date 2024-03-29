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

#include "List.H"
#include "ListLoopM.H"

#include "FixedList.H"
#include "PtrList.H"
#include "SLList.H"
#include "IndirectList.H"
#include "UIndirectList.H"
#include "BiIndirectList.H"
#include "contiguous.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class Type>
const Foam::List<Type> Foam::List<Type>::zero;


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

// Construct with length specified
template<class T>
Foam::List<T>::List(const label s)
:
	UList<T>(nullptr, s)
{
	if (this->size_ < 0)
	{
		FatalErrorInFunction
			<< "bad size " << this->size_
			<< abort(FatalError);
	}

	alloc();
}


// Construct with length and single value specified
template<class T>
Foam::List<T>::List(const label s, const T& a)
:
	UList<T>(nullptr, s)
{
	if (this->size_ < 0)
	{
		FatalErrorInFunction
			<< "bad size " << this->size_
			<< abort(FatalError);
	}

	alloc();

	if (this->size_)
	{
		List_ACCESS(T, (*this), vp);
		List_FOR_ALL((*this), i)
			List_ELEM((*this), vp, i) = a;
		List_END_FOR_ALL
	}
}


// Construct as copy
template<class T>
Foam::List<T>::List(const List<T>& a)
:
	UList<T>(nullptr, a.size_)
{
	if (this->size_)
	{
		alloc();

		#ifdef USEMEMCPY
		if (contiguous<T>())
		{
			memcpy(this->v_, a.v_, this->byteSize());
		}
		else
		#endif
		{
			List_ACCESS(T, (*this), vp);
			List_CONST_ACCESS(T, a, ap);
			List_FOR_ALL((*this), i)
				List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
			List_END_FOR_ALL
		}
	}
}


// Construct by transferring the parameter contents
template<class T>
Foam::List<T>::List(const Xfer< List<T> >& lst)
{
	transfer(lst());
}


// Construct as copy or re-use as specified.
template<class T>
Foam::List<T>::List(List<T>& a, bool reuse)
:
	UList<T>(nullptr, a.size_)
{
	if (reuse)
	{
		this->v_ = a.v_;
		a.v_ = 0;
		a.size_ = 0;
	}
	else if (this->size_)
	{
		alloc();

		#ifdef USEMEMCPY
		if (contiguous<T>())
		{
			memcpy(this->v_, a.v_, this->byteSize());
		}
		else
		#endif
		{
			List_ACCESS(T, (*this), vp);
			List_CONST_ACCESS(T, a, ap);
			List_FOR_ALL((*this), i)
				List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
			List_END_FOR_ALL
		}
	}
}


// Construct as subset
template<class T>
Foam::List<T>::List(const UList<T>& a, const unallocLabelList& map)
:
	UList<T>(nullptr, map.size())
{
	if (this->size_)
	{
		// Note:cannot use List_ELEM since third argument has to be index.

		alloc();

		forAll(*this, i)
		{
			this->operator[](i) = a[map[i]];
		}
	}
}


// Construct given start and end iterators.
template<class T>
template<class InputIterator>
Foam::List<T>::List(InputIterator first, InputIterator last)
:
	List<T>(first, last, std::distance(first, last))
{}


// Construct as copy of FixedList<T, Size>
template<class T>
template<unsigned Size>
Foam::List<T>::List(const FixedList<T, Size>& lst)
:
	UList<T>(nullptr, Size)
{
	allocCopyList(lst);
}


// Construct as copy of PtrList<T>
template<class T>
Foam::List<T>::List(const PtrList<T>& lst)
:
	UList<T>(nullptr, lst.size())
{
	allocCopyList(lst);
}


// Construct as copy of SLList<T>
template<class T>
Foam::List<T>::List(const SLList<T>& lst)
:
	UList<T>(nullptr, lst.size())
{
	if (this->size_)
	{
		alloc();

		label i = 0;
		for
		(
			typename SLList<T>::const_iterator iter = lst.begin();
			iter != lst.end();
			++iter
		)
		{
			this->operator[](i++) = iter();
		}
	}
}


// Construct as copy of IndirectList<T>
template<class T>
Foam::List<T>::List(const IndirectList<T>& lst)
:
	UList<T>(nullptr, lst.size())
{
	allocCopyList(lst);
}


// Construct as copy of UIndirectList<T>
template<class T>
Foam::List<T>::List(const UIndirectList<T>& lst)
:
	UList<T>(nullptr, lst.size())
{
	allocCopyList(lst);
}


// Construct as copy of BiIndirectList<T>
template<class T>
Foam::List<T>::List(const BiIndirectList<T>& lst)
:
	UList<T>(nullptr, lst.size())
{
	allocCopyList(lst);
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

// Destroy list elements
template<class T>
Foam::List<T>::~List()
{
	if (this->v_)
	{
		delete[] this->v_;
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::List<T>::setSize(const label newSize)
{
	if (newSize < 0)
	{
		FatalErrorInFunction
			<< "bad set size " << newSize
			<< abort(FatalError);
	}

	if (newSize != this->size_)
	{
		if (newSize > 0)
		{
			T* nv = new T[label(newSize)];

			if (this->size_)
			{
				label i = min(this->size_, newSize);

#				ifdef USEMEMCPY
				if (contiguous<T>())
				{
					memcpy(nv, this->v_, i*sizeof(T));
				}
				else
#				endif
				{
					T* vv = &this->v_[i];
					T* av = &nv[i];
					while (i--) *--av = *--vv;
				}
			}

			clear();
			this->size_ = newSize;
			this->v_ = nv;
		}
		else
		{
			clear();
		}
	}
}


template<class T>
void Foam::List<T>::setSize(const label newSize, const T& a)
{
	label oldSize = this->size_;
	this->setSize(newSize);

	if (newSize > oldSize)
	{
		label i = newSize - oldSize;
		T* vv = &this->v_[newSize];
		while (i--) *--vv = a;
	}
}


// Transfer the contents of the argument List into this List
// and anull the argument list
template<class T>
void Foam::List<T>::transfer(List<T>& a)
{
	clear();
	this->size_ = a.size_;
	this->v_ = a.v_;

	a.size_ = 0;
	a.v_ = 0;
}


// Transfer the contents of the argument DynamicList into this List
// and anull the argument list
template<class T>
template<unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
void Foam::List<T>::transfer(DynamicList<T, SizeInc, SizeMult, SizeDiv>& a)
{
	// shrink the allocated space to the number of elements used
	a.shrink();
	transfer(static_cast<List<T>&>(a));
	a.clearStorage();
}


// Transfer the contents of the argument SortableList into this List
// and anull the argument list
template<class T>
void Foam::List<T>::transfer(SortableList<T>& a)
{
	// shrink away the sort indices
	a.shrink();
	transfer(static_cast<List<T>&>(a));
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// Assignment to UList operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const UList<T>& a)
{
	reAlloc(a.size_);

	if (this->size_)
	{
#		ifdef USEMEMCPY
		if (contiguous<T>())
		{
			memcpy(this->v_, a.v_, this->byteSize());
		}
		else
#		endif
		{
			List_ACCESS(T, (*this), vp);
			List_CONST_ACCESS(T, a, ap);
			List_FOR_ALL((*this), i)
				List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
			List_END_FOR_ALL
		}
	}
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const List<T>& a)
{
	if (this == &a)
	{
		FatalErrorInFunction
			<< "attempted assignment to self"
			<< abort(FatalError);
	}

	operator=(static_cast<const UList<T>&>(a));
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const SLList<T>& lst)
{
	reAlloc(lst.size());

	if (this->size_)
	{
		label i = 0;
		for
		(
			typename SLList<T>::const_iterator iter = lst.begin();
			iter != lst.end();
			++iter
		)
		{
			this->operator[](i++) = iter();
		}
	}
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const IndirectList<T>& lst)
{
	reAlloc(lst.size());
	copyList(lst);
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const UIndirectList<T>& lst)
{
	reAlloc(lst.size());
	copyList(lst);
}


// Assignment operator. Takes linear time.
template<class T>
void Foam::List<T>::operator=(const BiIndirectList<T>& lst)
{
	reAlloc(lst.size());
	copyList(lst);
}

// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "ListIO.C"

// ************************************************************************* //
