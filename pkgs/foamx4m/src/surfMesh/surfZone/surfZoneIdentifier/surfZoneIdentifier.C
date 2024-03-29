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

#include "surfZoneIdentifier.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIdentifier::surfZoneIdentifier()
:
	name_(word::null),
	index_(0),
	geometricType_(word::null)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
	const word& name,
	const label index,
	const word& geometricType
)
:
	name_(name),
	index_(index),
	geometricType_(geometricType)
{}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
	const word& name,
	const dictionary& dict,
	const label index
)
:
	name_(name),
	index_(index)
{
	dict.readIfPresent("geometricType", geometricType_);
}


Foam::surfZoneIdentifier::surfZoneIdentifier
(
	const surfZoneIdentifier& p,
	const label index
)
:
	name_(p.name()),
	index_(index),
	geometricType_(p.geometricType())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfZoneIdentifier::~surfZoneIdentifier()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


void Foam::surfZoneIdentifier::write(Ostream& os) const
{
	if (geometricType_.size())
	{
		os.writeKeyword("geometricType") << geometricType_
			<< token::END_STATEMENT << nl;
	}
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

// bool Foam::surfZoneIdentifier::operator!=
// (
//	 const surfZoneIdentifier& p
// ) const
// {
//	 return !(*this == p);
// }
//
//
// bool Foam::surfZoneIdentifier::operator==
// (
//	 const surfZoneIdentifier& p
// ) const
// {
//	 return geometricType() == p.geometricType() && name() == p.name();
// }


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Foam::Istream& Foam::operator>>(Istream& is, surfZoneIdentifier& p)
// {
//	 is >> p.name_ >> p.geometricType_;
//
//	 return is;
// }


Foam::Ostream& Foam::operator<<(Ostream& os, const surfZoneIdentifier& p)
{
	p.write(os);
	os.check
	(
		"Ostream& operator<<(Ostream&, const surfZoneIdentifier&)"
	);
	return os;
}


// ************************************************************************* //
