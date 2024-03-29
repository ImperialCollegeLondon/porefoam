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

#include "CompactIOList.H"
#include "labelList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::readFromStream()
{
	Istream& is = readStream(word::null);

	if (headerClassName() == IOList<T>::typeName)
	{
		is >> static_cast<List<T>&>(*this);
		close();
	}
	else if (headerClassName() == typeName)
	{
		is >> *this;
		close();
	}
	else
	{
		FatalIOErrorIn
		(
			"CompactIOList<T, BaseType>::readFromStream()",
			is
		)   << "unexpected class name " << headerClassName()
			<< " expected " << typeName << " or " << IOList<T>::typeName
			<< endl
			<< "    while reading object " << name()
			<< exit(FatalIOError);
	}
}


// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList(const IOobject& io)
:
	regIOobject(io)
{
	if
	(
		io.readOpt() == IOobject::MUST_READ
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readFromStream();
	}
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
	const IOobject& io,
	const label size
)
:
	regIOobject(io)
{
	if
	(
		io.readOpt() == IOobject::MUST_READ
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readFromStream();
	}
	else
	{
		List<T>::setSize(size);
	}
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
	const IOobject& io,
	const List<T>& list
)
:
	regIOobject(io)
{
	if
	(
		io.readOpt() == IOobject::MUST_READ
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readFromStream();
	}
	else
	{
		List<T>::operator=(list);
	}
}


template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::CompactIOList
(
	const IOobject& io,
	const Xfer<List<T> >& list
)
:
	regIOobject(io)
{
	List<T>::transfer(list());

	if
	(
		io.readOpt() == IOobject::MUST_READ
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readFromStream();
	}
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::CompactIOList<T, BaseType>::~CompactIOList()
{}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::writeObject
(
	IOstream::streamFormat fmt,
	IOstream::versionNumber ver,
	IOstream::compressionType cmp
) const
{
	if (fmt == IOstream::ASCII)
	{
		// Change type to be non-compact format type
		const word oldTypeName = typeName;

		const_cast<word&>(typeName) = IOList<T>::typeName;

		bool good = regIOobject::writeObject(fmt, ver, cmp);

		// Change type back
		const_cast<word&>(typeName) = oldTypeName;

		return good;
	}
	else
	{
		return regIOobject::writeObject(fmt, ver, cmp);
	}
}


template<class T, class BaseType>
bool Foam::CompactIOList<T, BaseType>::writeData(Ostream& os) const
{
	return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=
(
	const CompactIOList<T, BaseType>& rhs
)
{
	List<T>::operator=(rhs);
}


template<class T, class BaseType>
void Foam::CompactIOList<T, BaseType>::operator=(const List<T>& rhs)
{
	List<T>::operator=(rhs);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class T, class BaseType>
Foam::Istream& Foam::operator>>
(
	Foam::Istream& is,
	Foam::CompactIOList<T, BaseType>& L
)
{
	// Read compact
	const labelList start(is);
	const List<BaseType> elems(is);

	// Convert
	L.setSize(start.size()-1);

	forAll(L, i)
	{
		T& subList = L[i];

		label index = start[i];
		subList.setSize(start[i+1] - index);

		forAll(subList, j)
		{
			subList[j] = elems[index++];
		}
	}

	return is;
}


template<class T, class BaseType>
Foam::Ostream& Foam::operator<<
(
	Foam::Ostream& os,
	const Foam::CompactIOList<T, BaseType>& L
)
{
	// Keep ascii writing same.
	if (os.format() == IOstream::ASCII)
	{
		os << static_cast<const List<T>&>(L);
	}
	else
	{
		// Convert to compact format
		labelList start(L.size()+1);

		start[0] = 0;
		for (label i = 1; i < start.size(); i++)
		{
			start[i] = start[i-1]+L[i-1].size();
		}

		List<BaseType> elems(start[start.size()-1]);

		label elemI = 0;
		forAll(L, i)
		{
			const T& subList = L[i];

			forAll(subList, j)
			{
				elems[elemI++] = subList[j];
			}
		}
		os << start << elems;
	}

	return os;
}


// ************************************************************************* //
