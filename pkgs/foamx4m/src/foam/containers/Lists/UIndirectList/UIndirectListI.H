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
inline Foam::UIndirectList<T>::UIndirectList
(
	const UList<T>& completeList,
	const UList<label>& addr
)
:
	completeList_(const_cast<UList<T>&>(completeList)),
	addressing_(addr)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
inline Foam::label Foam::UIndirectList<T>::size() const
{
	return addressing_.size();
}


template<class T>
inline bool Foam::UIndirectList<T>::empty() const
{
	return addressing_.empty();
}


template<class T>
inline const Foam::UList<T>& Foam::UIndirectList<T>::completeList() const
{
	return completeList_;
}


template<class T>
inline const Foam::List<Foam::label>& Foam::UIndirectList<T>::addressing() const
{
	return addressing_;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
inline Foam::List<T> Foam::UIndirectList<T>::operator()() const
{
	List<T> result(size());

	forAll(*this, i)
	{
		result[i] = operator[](i);
	}

	return result;
}


template<class T>
inline T& Foam::UIndirectList<T>::operator[](const label i)
{
	return completeList_[addressing_[i]];
}


template<class T>
inline const T& Foam::UIndirectList<T>::operator[](const label i) const
{
	return completeList_[addressing_[i]];
}


template<class T>
inline void Foam::UIndirectList<T>::operator=(const UList<T>& ae)
{
	if (addressing_.size() != ae.size())
	{
		FatalErrorIn("UIndirectList<T>::operator=(const UList<T>&)")
			<< "Addressing and list of addressed elements "
			   "have different sizes: "
			<< addressing_.size() << " " << ae.size()
			<< abort(FatalError);
	}

	forAll(addressing_, i)
	{
		completeList_[addressing_[i]] = ae[i];
	}
}


template<class T>
inline void Foam::UIndirectList<T>::operator=(const T& t)
{
	forAll(addressing_, i)
	{
		completeList_[addressing_[i]] = t;
	}
}


// ************************************************************************* //
