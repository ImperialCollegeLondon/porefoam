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

#include "zoneToPoint.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(zoneToPoint, 0);

addToRunTimeSelectionTable(topoSetSource, zoneToPoint, word);

addToRunTimeSelectionTable(topoSetSource, zoneToPoint, istream);

}


Foam::topoSetSource::addToUsageTable Foam::zoneToPoint::usage_
(
	zoneToPoint::typeName,
	"\n    Usage: zoneToPoint zone\n\n"
	"    Select all points in the pointZone."
	" Note:accepts wildcards for zone.\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::zoneToPoint::combine(topoSet& set, const bool add) const
{
	bool hasMatched = false;

	forAll(mesh_.pointZones(), i)
	{
		const pointZone& zone = mesh_.pointZones()[i];

		if (zoneName_.match(zone.name()))
		{
			const labelList& pointLabels = mesh_.pointZones()[i];

			Info<< "    Found matching zone " << zone.name()
				<< " with " << pointLabels.size() << " points." << endl;

			hasMatched = true;

			forAll(pointLabels, i)
			{
				// Only do active points
				if (pointLabels[i] < mesh_.nPoints())
				{
					addOrDelete(set, pointLabels[i], add);
				}
			}
		}
	}

	if (!hasMatched)
	{
		WarningIn("zoneToPoint::combine(topoSet&, const bool)")
			<< "Cannot find any pointZone named " << zoneName_ << endl
			<< "Valid names are " << mesh_.pointZones().names() << endl;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::zoneToPoint::zoneToPoint
(
	const polyMesh& mesh,
	const word& zoneName
)
:
	topoSetSource(mesh),
	zoneName_(zoneName)
{}


// Construct from dictionary
Foam::zoneToPoint::zoneToPoint
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	zoneName_(dict.lookup("name"))
{}


// Construct from Istream
Foam::zoneToPoint::zoneToPoint
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	zoneName_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::zoneToPoint::~zoneToPoint()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zoneToPoint::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding all points of pointZone " << zoneName_ << " ..."
			<< endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing all points of pointZone " << zoneName_ << " ..."
			<< endl;

		combine(set, false);
	}
}


// ************************************************************************* //
