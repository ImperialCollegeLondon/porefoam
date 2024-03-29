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

#include "septernion.H"
#include "IOstreams.H"
#include "OStringStream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::septernion::typeName = "septernion";
const Foam::septernion Foam::septernion::zero(vector::zero, quaternion::zero);
const Foam::septernion Foam::septernion::I(vector::zero, quaternion::I);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::septernion::septernion(Istream& is)
{
	operator>>(is, *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const septernion& s)
{
	OStringStream buf;
	buf << '(' << s.t() << ',' << s.r() << ')';
	return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, septernion& s)
{
	// Read beginning of septernion
	is.readBegin("septernion");

	is  >> s.t() >> s.r();

	// Read end of septernion
	is.readEnd("septernion");

	// Check state of Istream
	is.check("operator>>(Istream&, septernion&)");

	return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const septernion& s)
{
	os  << token::BEGIN_LIST
		<< s.t() << token::SPACE << s.r()
		<< token::END_LIST;

	return os;
}


// ************************************************************************* //
