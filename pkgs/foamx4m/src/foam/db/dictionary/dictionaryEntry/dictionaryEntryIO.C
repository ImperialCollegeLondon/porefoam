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

#include "keyType.H"
#include "dictionaryEntry.H"
#include "IOstreams.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dictionaryEntry::dictionaryEntry
(
	const dictionary& parentDict,
	Istream& is
)
:
	entry(keyType(is)),
	dictionary(parentDict, is)
{
	is.fatalCheck
	(
		"dictionaryEntry::dictionaryEntry"
		"(const dictionary& parentDict, Istream&)"
	);
}


Foam::dictionaryEntry::dictionaryEntry
(
	const keyType& key,
	const dictionary& parentDict,
	Istream& is
)
:
	entry(key),
	dictionary(key, parentDict, is)
{
	is.fatalCheck
	(
		"dictionaryEntry::dictionaryEntry"
		"(const keyType&, const dictionary& parentDict, Istream&)"
	);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::dictionaryEntry::write(Ostream& os) const
{
	// write keyword with indent but without trailing spaces
	os.indent();
	os.write(keyword());
	dictionary::write(os);
}


// * * * * * * * * * * * * * * Ostream operator  * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const dictionaryEntry& de)
{
	de.write(os);
	return os;
}


#if defined (__GNUC__)
template<>
#endif
Foam::Ostream& Foam::operator<<
(
	Ostream& os,
	const InfoProxy<dictionaryEntry>& ip
)
{
	const dictionaryEntry& e = ip.t_;

	os  << "    dictionaryEntry '" << e.keyword() << "'" << endl;

	return os;
}


// ************************************************************************* //
