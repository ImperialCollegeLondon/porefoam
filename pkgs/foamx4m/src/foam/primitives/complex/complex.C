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

#include "complex.H"
#include "IOstreams.H"

#include <sstream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const char* const Foam::complex::typeName = "complex";
const Foam::complex Foam::complex::zero(0, 0);
const Foam::complex Foam::complex::one(1, 1);

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::complex::complex(Istream& is)
{
	is >> *this;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::word Foam::name(const complex& c)
{
	std::ostringstream buf;
	buf << '(' << c.Re() << ',' << c.Im() << ')';
	return buf.str();
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, complex& c)
{
	// Read beginning of complex
	is.readBegin("complex");

	is  >> c.re >> c.im;

	// Read end of complex
	is.readEnd("complex");

	// Check state of Istream
	is.check("operator>>(Istream&, complex&)");

	return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const complex& c)
{
	os  << token::BEGIN_LIST
		<< c.re << token::SPACE << c.im
		<< token::END_LIST;

	return os;
}


// ************************************************************************* //
