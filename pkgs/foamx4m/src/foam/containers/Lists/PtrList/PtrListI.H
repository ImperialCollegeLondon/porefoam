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

#include "error.H"

#include "autoPtr.H"
#include "tmp.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::PtrList<T>::size() const
{
	return ptrs_.size();
}


template<class T>
inline bool Foam::PtrList<T>::empty() const
{
	return ptrs_.empty();
}


template<class T>
inline void Foam::PtrList<T>::resize(const label newSize)
{
	this->setSize(newSize);
}


template<class T>
inline bool Foam::PtrList<T>::set(const label i) const
{
	return ptrs_[i] != nullptr;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set(const label i, T* ptr)
{
	autoPtr<T> old(ptrs_[i]);

	ptrs_[i] = ptr;

	return old;
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
	const label i,
	const autoPtr<T>& aptr
)
{
	return set(i, const_cast<autoPtr<T>&>(aptr).ptr());
}


template<class T>
inline Foam::autoPtr<T> Foam::PtrList<T>::set
(
	const label i,
	const tmp<T>& t
)
{
	return set(i, const_cast<tmp<T>&>(t).ptr());
}


template<class T>
inline Foam::Xfer<Foam::PtrList<T> > Foam::PtrList<T>::xfer()
{
	return xferMove(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
const T& Foam::PtrList<T>::operator[](const label i) const
{
	if (!ptrs_[i])
	{
		FatalErrorIn("PtrList::operator[] const")
			<< "hanging pointer, cannot dereference"
			<< abort(FatalError);
	}

	return *(ptrs_[i]);
}


template<class T>
T& Foam::PtrList<T>::operator[](const label i)
{
	if (!ptrs_[i])
	{
		FatalErrorIn("PtrList::operator[]")
			<< "hanging pointer, cannot dereference"
			<< abort(FatalError);
	}

	return *(ptrs_[i]);
}


template<class T>
const T* Foam::PtrList<T>::operator()(const label i) const
{
	return ptrs_[i];
}


// * * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * //

template<class T>
inline Foam::PtrList<T>::iterator::iterator(T** ptr)
:
	ptr_(ptr)
{}

template<class T>
inline bool Foam::PtrList<T>::iterator::operator==(const iterator& iter) const
{
	return ptr_ == iter.ptr_;
}

template<class T>
inline bool Foam::PtrList<T>::iterator::operator!=(const iterator& iter) const
{
	return ptr_ != iter.ptr_;
}

template<class T>
inline T& Foam::PtrList<T>::iterator::operator*()
{
	return **ptr_;
}

template<class T>
inline T& Foam::PtrList<T>::iterator::operator()()
{
	return operator*();
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::iterator::operator++()
{
	++ptr_;
	return *this;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::iterator::operator++(int)
{
	iterator tmp = *this;
	++ptr_;
	return tmp;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::iterator::operator--()
{
	--ptr_;
	return *this;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::iterator::operator--(int)
{
	iterator tmp = *this;
	--ptr_;
	return tmp;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::iterator::operator+=(label n)
{
	ptr_ += n;
	return *this;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::operator+(const typename PtrList<T>::iterator& iter, label n)
{
	typename PtrList<T>::iterator tmp = iter;
	return tmp += n;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::operator+(label n, const typename PtrList<T>::iterator& iter)
{
	typename PtrList<T>::iterator tmp = iter;
	return tmp += n;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::iterator::operator-=(label n)
{
	ptr_ -= n;
	return *this;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::operator-(const typename PtrList<T>::iterator& iter, label n)
{
	typename PtrList<T>::iterator tmp = iter;
	return tmp -= n;
}

template<class T>
inline Foam::label Foam::operator-
(
	const typename PtrList<T>::iterator& iter1,
	const typename PtrList<T>::iterator& iter2
)
{
	return (iter1.ptr_ - iter2.ptr_)/sizeof(T*);
}

template<class T>
inline T& Foam::PtrList<T>::iterator::operator[](label n)
{
	return *(*this + n);
}

template<class T>
inline bool Foam::PtrList<T>::iterator::operator<(const iterator& iter) const
{
	return ptr_ < iter.ptr_;
}

template<class T>
inline bool Foam::PtrList<T>::iterator::operator>(const iterator& iter) const
{
	return ptr_ > iter.ptr_;
}

template<class T>
inline bool Foam::PtrList<T>::iterator::operator<=(const iterator& iter) const
{
	return ptr_ <= iter.ptr_;
}

template<class T>
inline bool Foam::PtrList<T>::iterator::operator>=(const iterator& iter) const
{
	return ptr_ >= iter.ptr_;
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::begin()
{
	return ptrs_.begin();
}

template<class T>
inline typename Foam::PtrList<T>::iterator
Foam::PtrList<T>::end()
{
	return ptrs_.end();
}



// ************************************************************************* //
