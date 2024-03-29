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

#include "Constant.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
	const word& entryName,
	const dictionary& dict
)
:
	Function1<Type>(entryName),
	value_(pTraits<Type>::zero)
{
	Istream& is(dict.lookup(entryName));
	word entryType(is);
	is  >> value_;
}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant
(
	const word& entryName,
	Istream& is
)
:
	Function1<Type>(entryName),
	value_(pTraits<Type>(is))
{}


template<class Type>
Foam::Function1Types::Constant<Type>::Constant(const Constant<Type>& cnst)
:
	Function1<Type>(cnst),
	value_(cnst.value_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::Function1Types::Constant<Type>::~Constant()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
Type Foam::Function1Types::Constant<Type>::value(const scalar x) const
{
	return value_;
}


template<class Type>
Type Foam::Function1Types::Constant<Type>::integrate
(
	const scalar x1,
	const scalar x2
) const
{
	return (x2 - x1)*value_;
}


template<class Type>
void Foam::Function1Types::Constant<Type>::writeData(Ostream& os) const
{
	Function1<Type>::writeData(os);

	os  << token::SPACE << value_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
