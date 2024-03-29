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

#include "word.H"
#include "Ostream.H"
#include "token.H"
#include "keyType.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Decrement the indent level
void Foam::Ostream::decrIndent()
{
	if (indentLevel_ == 0)
	{
		cerr<< "Ostream::decrIndent() : attempt to decrement 0 indent level"
			<< std::endl;
	}
	else
	{
		indentLevel_--;
	}
}


// Write keyType
// write regular expression as quoted string
// write plain word as word (unquoted)
Foam::Ostream& Foam::Ostream::write(const keyType& kw)
{
	return writeQuoted(kw, kw.isPattern());
}


// Write the keyword followed by appropriate indentation
Foam::Ostream& Foam::Ostream::writeKeyword(const keyType& kw)
{
	indent();
	write(kw);

	label nSpaces = entryIndentation_ - label(kw.size());

	// pattern is surrounded by quotes
	if (kw.isPattern())
	{
		nSpaces -= 2;
	}

	// could also increment by indentSize_ ...
	if (nSpaces < 1)
	{
		nSpaces = 1;
	}

	while (nSpaces--)
	{
		write(char(token::SPACE));
	}

	return *this;
}


// ************************************************************************* //
