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

#include "boxToFace.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(boxToFace, 0);

addToRunTimeSelectionTable(topoSetSource, boxToFace, word);

addToRunTimeSelectionTable(topoSetSource, boxToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::boxToFace::usage_
(
	boxToFace::typeName,
	"\n    Usage: boxToFace ((minx miny minz) (maxx maxy maxz))\n\n"
	"    Select all face with faceCentre within bounding box\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::boxToFace::combine(topoSet& set, const bool add) const
{
	const pointField& ctrs = mesh_.faceCentres();

	forAll(ctrs, faceI)
	{
		if (bb_.contains(ctrs[faceI]))
		{
			addOrDelete(set, faceI, add);
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::boxToFace::boxToFace
(
	const polyMesh& mesh,
	const treeBoundBox& bb
)
:
	topoSetSource(mesh),
	bb_(bb)
{}


// Construct from dictionary
Foam::boxToFace::boxToFace
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	bb_(dict.lookup("box"))
{}


// Construct from Istream
Foam::boxToFace::boxToFace
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	bb_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::boxToFace::~boxToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::boxToFace::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding faces with centre within box " << bb_ << endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing faces with centre within box " << bb_ << endl;

		combine(set, false);
	}
}


// ************************************************************************* //
