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
	Reads an bool from an input stream, for a given version number and file
	format. If an ASCII file is being read, then the line numbers are counted
	and an erroneous read is reported.

\*---------------------------------------------------------------------------*/

#include "error.H"

#include "bool.H"
#include "Switch.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::Istream& Foam::operator>>(Istream& is, bool& b)
{
	if (is.good())
	{
		b = Switch(is);
	}

	return is;
}


Foam::Ostream& Foam::operator<<(Ostream& os, const bool b)
{
	// we could also write as text string without any difficulty
	// os << Switch::asText(b);
	os.write(label(b));
	os.check("Ostream& operator<<(Ostream&, const bool)");
	return os;
}


bool Foam::readBool(Istream& is)
{
	bool val;
	is >> val;

	return val;
}


// ************************************************************************* //
