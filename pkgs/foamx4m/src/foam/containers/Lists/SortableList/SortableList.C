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

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class T>
void Foam::SortableList<T>::sortIndices(List<label>& order) const
{
	// list lengths must be identical
	if (order.size() != this->size())
	{
		// avoid copying any elements, they are overwritten anyhow
		order.clear();
		order.setSize(this->size());
	}

	forAll(order, elemI)
	{
		order[elemI] = elemI;
	}

	Foam::stableSort(order, typename UList<T>::less(*this));
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class T>
Foam::SortableList<T>::SortableList()
{}


template<class T>
Foam::SortableList<T>::SortableList(const UList<T>& values)
:
	List<T>(values)
{
	sort();
}


template<class T>
Foam::SortableList<T>::SortableList(const Xfer<List<T> >& values)
:
	List<T>(values)
{
	sort();
}


template<class T>
Foam::SortableList<T>::SortableList(const label size)
:
	List<T>(size)
{}


template<class T>
Foam::SortableList<T>::SortableList(const label size, const T& val)
:
	List<T>(size, val)
{}


template<class T>
Foam::SortableList<T>::SortableList(const SortableList<T>& lst)
:
	List<T>(lst),
	indices_(lst.indices())
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


template<class T>
void Foam::SortableList<T>::clear()
{
	List<T>::clear();
	indices_.clear();
}


template<class T>
Foam::List<T>& Foam::SortableList<T>::shrink()
{
	indices_.clear();
	return static_cast<List<T>&>(*this);
}


template<class T>
void Foam::SortableList<T>::sort()
{
	sortIndices(indices_);

	List<T> lst(this->size());
	forAll(indices_, i)
	{
		lst[i] = this->operator[](indices_[i]);
	}

	List<T>::transfer(lst);
}


template<class T>
void Foam::SortableList<T>::reverseSort()
{
	sortIndices(indices_);

	List<T> lst(this->size());
	label endI = indices_.size();
	forAll(indices_, i)
	{
		lst[--endI] = this->operator[](indices_[i]);
	}

	List<T>::transfer(lst);
}


template<class T>
Foam::Xfer<Foam::List<T> > Foam::SortableList<T>::xfer()
{
	return xferMoveTo<List<T> >(*this);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline void Foam::SortableList<T>::operator=(const T& t)
{
	UList<T>::operator=(t);
}


template<class T>
inline void Foam::SortableList<T>::operator=(const UList<T>& rhs)
{
	List<T>::operator=(rhs);
	indices_.clear();
}


template<class T>
inline void Foam::SortableList<T>::operator=(const SortableList<T>& rhs)
{
	List<T>::operator=(rhs);
	indices_ = rhs.indices();
}


// ************************************************************************* //
