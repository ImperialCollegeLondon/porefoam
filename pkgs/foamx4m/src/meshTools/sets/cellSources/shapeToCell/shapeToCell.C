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

#include "shapeToCell.H"
#include "polyMesh.H"
#include "mathematicalConstants.H"
#include "hexMatcher.H"
#include "cellFeatures.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(shapeToCell, 0);

addToRunTimeSelectionTable(topoSetSource, shapeToCell, word);

addToRunTimeSelectionTable(topoSetSource, shapeToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::shapeToCell::usage_
(
	shapeToCell::typeName,
	"\n    Usage: shapeToCell tet|pyr|prism|hex|tetWedge|wedge|splitHex\n\n"
	"    Select all cells of given cellShape.\n"
	"    (splitHex hardcoded with internal angle < 10 degrees)\n"
);


// Angle for polys to be considered splitHexes.
Foam::scalar Foam::shapeToCell::featureCos =
	Foam::cos(10.0 * mathematicalConstant::pi/180.0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::shapeToCell::combine(topoSet& set, const bool add) const
{
	if (type_ == "splitHex")
	{
		for (label cellI = 0; cellI < mesh_.nCells(); cellI++)
		{
			cellFeatures superCell(mesh_, featureCos, cellI);

			if (hexMatcher().isA(superCell.faces()))
			{
				addOrDelete(set, cellI, add);
			}
		}
	}
	else
	{
		const cellModel& wantedModel = *(cellModeller::lookup(type_));

		const cellShapeList& cellShapes = mesh_.cellShapes();

		forAll(cellShapes, cellI)
		{
			if (cellShapes[cellI].model() == wantedModel)
			{
				addOrDelete(set, cellI, add);
			}
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::shapeToCell::shapeToCell
(
	const polyMesh& mesh,
	const word& type
)
:
	topoSetSource(mesh),
	type_(type)
{
	if (!cellModeller::lookup(type_) && (type_ != "splitHex"))
	{
		FatalErrorIn
		(
			"shapeToCell::shapeToCell(const polyMesh&, const word&)"
		)   << "Illegal cell type " << type_ << exit(FatalError);
	}
}


// Construct from dictionary
Foam::shapeToCell::shapeToCell
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	type_(dict.lookup("type"))
{
	if (!cellModeller::lookup(type_) && (type_ != "splitHex"))
	{
		FatalErrorIn
		(
			"shapeToCell::shapeToCell(const polyMesh&, const dictionary&)"
		)   << "Illegal cell type " << type_ << exit(FatalError);
	}
}


// Construct from Istream
Foam::shapeToCell::shapeToCell
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	type_(checkIs(is))
{
	if (!cellModeller::lookup(type_) && (type_ != "splitHex"))
	{
		FatalErrorIn
		(
			"shapeToCell::shapeToCell(const polyMesh&, Istream&)"
		)   << "Illegal cell type " << type_ << exit(FatalError);
	}
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::shapeToCell::~shapeToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::shapeToCell::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding all cells of type " << type_ << " ..." << endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing all cells of type " << type_ << " ..." << endl;

		combine(set, false);
	}
}


// ************************************************************************* //
