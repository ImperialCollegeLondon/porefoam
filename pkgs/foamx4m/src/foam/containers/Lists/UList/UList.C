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

#include "UList.H"
#include "ListLoopM.H"
#include "contiguous.H"

#include <algorithm>

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::assign(const UList<T>& a)
{
	if (a.size_ != this->size_)
	{
		FatalErrorIn("UList<T>::assign(const UList<T>&)")
			<< "ULists have different sizes: "
			<< this->size_ << " " << a.size_
			<< abort(FatalError);
	}

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


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::operator=(const T& t)
{
	List_ACCESS(T, (*this), vp);
	List_FOR_ALL((*this), i)
		List_ELEM((*this), vp, i) = t;
	List_END_FOR_ALL
}


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T>
void Foam::UList<T>::swap(UList<T>& a)
{
	if (a.size_ != this->size_)
	{
		FatalErrorIn("UList<T>::swap(const UList<T>&)")
			<< "ULists have different sizes: "
			<< this->size_ << " " << a.size_
			<< abort(FatalError);
	}

	List_ACCESS(T, (*this), vp);
	List_ACCESS(T, a, ap);
	T tmp;
	List_FOR_ALL((*this), i)
		tmp = List_ELEM((*this), vp, i);
		List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
		List_ELEM(a, ap, i) = tmp;
	List_END_FOR_ALL
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
Foam::label Foam::UList<T>::byteSize() const
{
	if (!contiguous<T>())
	{
		FatalErrorIn("UList<T>::byteSize()")
			<< "Cannot return the binary size of a list of "
			   "non-primitive elements"
			<< abort(FatalError);
	}

	return this->size_*sizeof(T);
}


template<class T>
void Foam::sort(UList<T>& a)
{
	std::sort(a.begin(), a.end());
}


template<class T, class Cmp>
void Foam::sort(UList<T>& a, const Cmp& cmp)
{
	std::sort(a.begin(), a.end(), cmp);
}


template<class T>
void Foam::stableSort(UList<T>& a)
{
	std::stable_sort(a.begin(), a.end());
}


template<class T, class Cmp>
void Foam::stableSort(UList<T>& a, const Cmp& cmp)
{
	std::stable_sort(a.begin(), a.end(), cmp);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
bool Foam::UList<T>::operator==(const UList<T>& a) const
{
	if (this->size_ != a.size_)
	{
		return false;
	}

	bool equal = true;

	List_CONST_ACCESS(T, (*this), vp);
	List_CONST_ACCESS(T, (a), ap);

	List_FOR_ALL((*this), i)
		equal = equal && (List_ELEM((*this), vp, i) == List_ELEM((a), ap, i));
	List_END_FOR_ALL

	return equal;
}


template<class T>
bool Foam::UList<T>::operator!=(const UList<T>& a) const
{
	return !operator==(a);
}


template<class T>
bool Foam::UList<T>::operator<(const UList<T>& a) const
{
	for
	(
		const_iterator vi = begin(), ai = a.begin();
		vi < end() && ai < a.end();
		vi++, ai++
	)
	{
		if (*vi < *ai)
		{
			return true;
		}
		else if (*vi > *ai)
		{
			return false;
		}
	}

	if (this->size_ < a.size_)
	{
		return true;
	}
	else
	{
		return false;
	}
}


template<class T>
bool Foam::UList<T>::operator>(const UList<T>& a) const
{
	return a.operator<(*this);
}


template<class T>
bool Foam::UList<T>::operator<=(const UList<T>& a) const
{
	return !operator>(a);
}


template<class T>
bool Foam::UList<T>::operator>=(const UList<T>& a) const
{
	return !operator<(a);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "UListIO.C"

// ************************************************************************* //
