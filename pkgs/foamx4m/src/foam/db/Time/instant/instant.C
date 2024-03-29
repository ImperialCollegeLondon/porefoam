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

#include "instantList.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::instant::typeName = "instant";

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::instant::instant()
{}

Foam::instant::instant(const scalar val, const word& tname)
:
	value_(val),
	name_(tname)
{}

Foam::instant::instant(const scalar val)
:
	value_(val),
	name_(Time::timeName(val))
{}

Foam::instant::instant(const word& tname)
:
	value_(atof(tname.c_str())),
	name_(tname)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::instant::equal(const scalar b) const
{
	return (value_ < b + SMALL  && value_ > b - SMALL);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

bool Foam::operator==(const instant& a, const instant& b)
{
	return a.equal(b.value_);
}


bool Foam::operator!=(const instant& a, const instant& b)
{
	return !operator==(a, b);
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, instant& I)
{
	is >> I.value_ >> I.name_;

	return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const instant& I)
{
   os << I.value_ << tab << I.name_;

   return os;
}


// ************************************************************************* //
