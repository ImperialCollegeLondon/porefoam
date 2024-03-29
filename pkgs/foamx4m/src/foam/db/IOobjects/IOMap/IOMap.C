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

#include "IOMap.H"

// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io)
:
	regIOobject(io)
{
	if
	(
		(
			io.readOpt() == IOobject::MUST_READ
		 || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
		)
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readStream(typeName) >> *this;
		close();
	}
}

template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const label size)
:
	regIOobject(io)
{
	if
	(
		(
			io.readOpt() == IOobject::MUST_READ
		 || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
		)
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readStream(typeName) >> *this;
		close();
	}
	else
	{
		Map<T>::setSize(size);
	}
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const Map<T>& map)
:
	regIOobject(io)
{
	if
	(
		(
			io.readOpt() == IOobject::MUST_READ
		 || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
		)
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readStream(typeName) >> *this;
		close();
	}
	else
	{
		Map<T>::operator=(map);
	}
}


template<class T>
Foam::IOMap<T>::IOMap(const IOobject& io, const Xfer<Map<T> >& map)
:
	regIOobject(io)
{
	Map<T>::transfer(map());

	if
	(
		(
			io.readOpt() == IOobject::MUST_READ
		 || io.readOpt() == IOobject::MUST_READ_IF_MODIFIED
		)
	 || (io.readOpt() == IOobject::READ_IF_PRESENT && headerOk())
	 || (io.readOpt() == IOobject::READ_IF_PRESENT_IF_MODIFIED && headerOk())
	)
	{
		readStream(typeName) >> *this;
		close();
	}
}


// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class T>
Foam::IOMap<T>::~IOMap()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
bool Foam::IOMap<T>::writeData(Ostream& os) const
{
	return (os << *this).good();
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class T>
void Foam::IOMap<T>::operator=(const IOMap<T>& rhs)
{
	Map<T>::operator=(rhs);
}


template<class T>
void Foam::IOMap<T>::operator=(const Map<T>& rhs)
{
	Map<T>::operator=(rhs);
}


// ************************************************************************* //
