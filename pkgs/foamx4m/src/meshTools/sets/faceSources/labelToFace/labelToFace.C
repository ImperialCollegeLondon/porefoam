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

#include "labelToFace.H"
#include "polyMesh.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(labelToFace, 0);

addToRunTimeSelectionTable(topoSetSource, labelToFace, word);

addToRunTimeSelectionTable(topoSetSource, labelToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::labelToFace::usage_
(
	labelToFace::typeName,
	"\n    Usage: labelToFace (i0 i1 .. in)\n\n"
	"    Select faces by label\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::labelToFace::combine(topoSet& set, const bool add) const
{
	forAll(labels_, labelI)
	{
		addOrDelete(set, labels_[labelI], add);
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::labelToFace::labelToFace
(
	const polyMesh& mesh,
	const labelList& labels
)
:
	topoSetSource(mesh),
	labels_(labels)
{}


// Construct from dictionary
Foam::labelToFace::labelToFace
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	labels_(dict.lookup("value"))
{}


// Construct from Istream
Foam::labelToFace::labelToFace
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	labels_(checkIs(is))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::labelToFace::~labelToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::labelToFace::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding faces mentioned in dictionary" << " ..." << endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing faces mentioned dictionary" << " ..." << endl;

		combine(set, false);
	}
}


// ************************************************************************* //
