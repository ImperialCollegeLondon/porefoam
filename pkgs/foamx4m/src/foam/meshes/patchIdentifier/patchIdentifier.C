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

#include "patchIdentifier.H"
#include "dictionary.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchIdentifier::patchIdentifier
(
	const word& name,
	const label index,
	const word& physicalType,
	const wordList& inGroups
)
:
	name_(name),
	index_(index),
	physicalType_(physicalType),
	inGroups_(inGroups)
{}


Foam::patchIdentifier::patchIdentifier
(
	const word& name,
	const dictionary& dict,
	const label index
)
:
	name_(name),
	index_(index)
{
	dict.readIfPresent("physicalType", physicalType_);
	dict.readIfPresent("inGroups", inGroups_);
}


Foam::patchIdentifier::patchIdentifier
(
	const patchIdentifier& p,
	const label index
)
:
	name_(p.name_),
	index_(index),
	physicalType_(p.physicalType_),
	inGroups_(p.inGroups_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::patchIdentifier::~patchIdentifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::patchIdentifier::inGroup(const word& name) const
{
	return findIndex(inGroups_, name) != -1;
}


void Foam::patchIdentifier::write(Ostream& os) const
{
	if (physicalType_.size())
	{
		os.writeKeyword("physicalType") << physicalType_
			<< token::END_STATEMENT << nl;
	}
	if (inGroups_.size())
	{
		os.writeKeyword("inGroups") << inGroups_
			<< token::END_STATEMENT << nl;
	}
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const patchIdentifier& pi)
{
	pi.write(os);
	os.check("Ostream& operator<<(Ostream&, const patchIdentifier&)");
	return os;
}


// ************************************************************************* //
