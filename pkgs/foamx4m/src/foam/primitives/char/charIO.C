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

Description
	Reads a char from an input stream, for a given version
	number and File format. If an ascii File is being read, then the line
	numbers are counted and an erroneous read is reported.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "char.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, char& c)
{
	is.read(c);
	is.check("Istream& operator>>(Istream&, char&)");
	return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const char c)
{
	os.write(c);
	os.check("Ostream& operator<<(Ostream&, const char)");
	return os;
}


char Foam::readChar(Istream& is)
{
   char c;
   is.read(c);
   return c;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const char* s)
{
	os.write(s);
	os.check("Ostream& operator<<(Ostream&, const char*)");
	return os;
}


// ************************************************************************* //
