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

#include "surfZoneIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfZoneIOList, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::surfZoneIOList::surfZoneIOList
(
	const IOobject& io
)
:
	surfZoneList(),
	regIOobject(io)
{
	Foam::string functionName =
		"surfZoneIOList::surfZoneIOList"
		"(const IOobject& io)";


	if (readOpt() == IOobject::MUST_READ)
	{
		surfZoneList& zones = *this;

		Istream& is = readStream(typeName);

		PtrList<entry> dictEntries(is);
		zones.setSize(dictEntries.size());

		label faceI = 0;
		forAll(zones, zoneI)
		{
			const dictionary& dict = dictEntries[zoneI].dict();

			label zoneSize = readLabel(dict.lookup("nFaces"));
			label startFaceI = readLabel(dict.lookup("startFace"));

			zones[zoneI] = surfZone
			(
				dictEntries[zoneI].keyword(),
				zoneSize,
				startFaceI,
				zoneI
			);

			word geoType;
			if (dict.readIfPresent("geometricType", geoType))
			{
				zones[zoneI].geometricType() = geoType;
			}

			if (startFaceI != faceI)
			{
				FatalErrorIn(functionName)
					<< "surfZones are not ordered. Start of zone " << zoneI
					<< " does not correspond to sum of preceding zones." << nl
					<< "while reading " << io.objectPath() << endl
					<< exit(FatalError);
			}

			faceI += zoneSize;
		}

		// Check state of IOstream
		is.check(functionName.c_str());

		close();
	}
}


Foam::surfZoneIOList::surfZoneIOList
(
	const IOobject& io,
	const surfZoneList& zones
)
:
	surfZoneList(zones),
	regIOobject(io)
{}


Foam::surfZoneIOList::surfZoneIOList
(
	const IOobject& io,
	const Xfer<surfZoneList>& zones
)
:
	surfZoneList(zones),
	regIOobject(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfZoneIOList::~surfZoneIOList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// writeData member function required by regIOobject
bool Foam::surfZoneIOList::writeData(Ostream& os) const
{
	os  << *this;
	return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

namespace Foam
{

	Ostream& operator<<(Ostream& os, const surfZoneIOList& L)
	{
		os  << L.size() << nl << token::BEGIN_LIST << incrIndent << nl;

		forAll(L, i)
		{
			L[i].writeDict(os);
		}

		os  << decrIndent << token::END_LIST;

		return os;
	}

}
// ************************************************************************* //
