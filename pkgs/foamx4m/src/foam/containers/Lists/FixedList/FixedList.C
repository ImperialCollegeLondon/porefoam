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

#include "FixedList.H"
#include "ListLoopM.H"

// * * * * * * * * * * * * * * * Static Members  * * * * * * * * * * * * * * //

template<class T, unsigned Size>
const Foam::FixedList<T, Size> Foam::FixedList<T, Size>::zero;


// * * * * * * * * * * * * * * STL Member Functions  * * * * * * * * * * * * //

template<class T, unsigned Size>
void Foam::FixedList<T, Size>::swap(FixedList<T, Size>& a)
{
	List_ACCESS(T, (*this), vp);
	List_ACCESS(T, a, ap);
	T tmp;
	List_FOR_ALL((*this), i)
		tmp = List_ELEM((*this), vp, i);
		List_ELEM((*this), vp, i) = List_ELEM(a, ap, i);
		List_ELEM(a, ap, i) = tmp;
	List_END_FOR_ALL
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator==(const FixedList<T, Size>& a) const
{
	bool equal = true;

	List_CONST_ACCESS(T, (*this), vp);
	List_CONST_ACCESS(T, (a), ap);

	List_FOR_ALL((*this), i)
		equal = equal && (List_ELEM((*this), vp, i) == List_ELEM((a), ap, i));
	List_END_FOR_ALL

	return equal;
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator!=(const FixedList<T, Size>& a) const
{
	return !operator==(a);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator<(const FixedList<T, Size>& a) const
{
	for
	(
		const_iterator vi = cbegin(), ai = a.cbegin();
		vi < cend() && ai < a.cend();
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

	if (Size < a.Size)
	{
		return true;
	}
	else
	{
		return false;
	}
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator>(const FixedList<T, Size>& a) const
{
	return a.operator<(*this);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator<=(const FixedList<T, Size>& a) const
{
	return !operator>(a);
}


template<class T, unsigned Size>
bool Foam::FixedList<T, Size>::operator>=(const FixedList<T, Size>& a) const
{
	return !operator<(a);
}


// * * * * * * * * * * * * * * * *  IOStream operators * * * * * * * * * * * //

#include "FixedListIO.C"

// ************************************************************************* //
