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

\*---------------------------------------------------------------------------*/

#include "refineCell.H"
#include "error.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null
Foam::refineCell::refineCell()
:
	cellNo_(-1),
	direction_(vector::zero)
{}


// from components
Foam::refineCell::refineCell(const label cellI, const vector& direction)
:
	cellNo_(cellI),
	direction_(direction)
{
	scalar magDir = mag(direction_);

	if (magDir < SMALL)
	{
		FatalErrorIn("refineCell(const label, const vector&)")
			<< "(almost)zero vector as direction for cell " << cellNo_
			<< abort(FatalError);
	}
	else if (mag(magDir - 1) > SMALL)
	{
		// Normalize
		direction_ /= mag(direction_);
	}
}


// from Istream
Foam::refineCell::refineCell(Istream& is)
:
	cellNo_(readLabel(is)),
	direction_(is)
{
	scalar magDir = mag(direction_);

	if (magDir < SMALL)
	{
		FatalErrorIn("refineCell(Istream&)")
			<< "(almost)zero vector as direction for cell " << cellNo_
			<< abort(FatalError);
	}
	else if (mag(magDir - 1) > SMALL)
	{
		// Normalize
		direction_ /= mag(direction_);
	}
}


// * * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const refineCell& r)
{
	if (os.format() == IOstream::ASCII)
	{
		os << r.cellNo() << token::SPACE << r.direction();
	}
	else
	{
		os << r.cellNo() << r.direction();
	}
	return os;
}


// ************************************************************************* //

