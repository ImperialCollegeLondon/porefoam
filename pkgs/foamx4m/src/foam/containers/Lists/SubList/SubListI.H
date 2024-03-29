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

template<class T>
inline Foam::SubList<T>::SubList
(
	const UList<T>& list,
	const label subSize
)
:
	UList<T>(list.v_, subSize)
{
#	ifdef FULLDEBUG
	list.checkSize(subSize);
#	endif
}


template<class T>
inline Foam::SubList<T>::SubList
(
	const UList<T>& list,
	const label subSize,
	const label startIndex
)
:
	UList<T>(&(list.v_[startIndex]), subSize)
{
#	ifdef FULLDEBUG

	// Artificially allow the start of a zero-sized subList to be
	// one past the end of the original list.
	if (subSize)
	{
		list.checkStart(startIndex);
		list.checkSize(startIndex + subSize);
	}
	else
	{
		// Start index needs to fall between 0 and size.  One position
		// behind the last element is allowed
		list.checkSize(startIndex);
	}
#	endif
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline void Foam::SubList<T>::reset
(
	const UList<T>& list,
	const label subSize
)
{
	UList<T>::reset(list.v_, subSize);

#	ifdef FULLDEBUG
	list.checkSize(subSize);
#	endif
}


template<class T>
inline void Foam::SubList<T>::reset
(
	const UList<T>& list,
	const label subSize,
	const label startIndex
)
{
	UList<T>::reset(&(list.v_[startIndex]), subSize);

#	ifdef FULLDEBUG

	// Artificially allowing the start of a zero-sized subList to be
	// one past the end of the original list.
	if (subSize > 0)
	{
		list.checkStart(startIndex);
		list.checkSize(startIndex + subSize);
	}
	else
	{
		// Start index needs to fall between 0 and size.  One position
		// behind the last element is allowed
		list.checkSize(startIndex);
	}
#	endif
}


template<class T>
inline const Foam::SubList<T>& Foam::SubList<T>::null()
{
	return *reinterpret_cast< SubList<T>* >(0);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::SubList<T>::operator const Foam::List<T>&() const
{
	// Dangerous, but no null constructor.  Reconsider.  HJ, 10/Apr/2013
	return *reinterpret_cast< const List<T>* >(this);
}


template<class T>
inline void Foam::SubList<T>::operator=(const T& t)
{
	UList<T>::operator=(t);
}


// ************************************************************************* //
