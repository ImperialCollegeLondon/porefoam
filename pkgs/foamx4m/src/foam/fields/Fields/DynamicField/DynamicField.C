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

#include "DynamicField.H"

// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::DynamicField(Istream& is)
:
	Field<T>(is),
	capacity_(Field<T>::size())
{}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::tmp<Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv> >
Foam::DynamicField<T, SizeInc, SizeMult, SizeDiv>::clone() const
{
	return tmp<DynamicField<T, SizeInc, SizeMult, SizeDiv> >
	(
		new DynamicField<T, SizeInc, SizeMult, SizeDiv>(*this)
	);
}


// * * * * * * * * * * * * * * * IOstream Operator * * * * * * * * * * * * * //

template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::Ostream& Foam::operator<<
(
	Ostream& os,
	const DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
	os << static_cast<const Field<T>&>(lst);
	return os;
}


template<class T, unsigned SizeInc, unsigned SizeMult, unsigned SizeDiv>
Foam::Istream& Foam::operator>>
(
	Istream& is,
	DynamicField<T, SizeInc, SizeMult, SizeDiv>& lst
)
{
	is >> static_cast<Field<T>&>(lst);
	lst.capacity_ = lst.Field<T>::size();

	return is;
}


// ************************************************************************* //
