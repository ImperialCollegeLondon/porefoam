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

#include "cylinderToCell.H"
#include "polyMesh.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(cylinderToCell, 0);
	addToRunTimeSelectionTable(topoSetSource, cylinderToCell, word);
	addToRunTimeSelectionTable(topoSetSource, cylinderToCell, istream);
}


Foam::topoSetSource::addToUsageTable Foam::cylinderToCell::usage_
(
	cylinderToCell::typeName,
	"\n    Usage: cylinderToCell (p1X p1Y p1Z) (p2X p2Y p2Z) radius\n\n"
	"    Select all cells with cell centre within bounding cylinder\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cylinderToCell::combine(topoSet& set, const bool add) const
{
	const vector axis = p2_ - p1_;
	const scalar rad2 = sqr(radius_);
	const scalar magAxis2 = magSqr(axis);

	const pointField& ctrs = mesh_.cellCentres();

	forAll(ctrs, cellI)
	{
		vector d = ctrs[cellI] - p1_;
		scalar magD = d & axis;

		if ((magD > 0) && (magD < magAxis2))
		{
			scalar d2 = (d & d) - sqr(magD)/magAxis2;
			if (d2 < rad2)
			{
				addOrDelete(set, cellI, add);
			}
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::cylinderToCell::cylinderToCell
(
	const polyMesh& mesh,
	const vector& p1,
	const vector& p2,
	const scalar radius
)
:
	topoSetSource(mesh),
	p1_(p1),
	p2_(p2),
	radius_(radius)
{}


Foam::cylinderToCell::cylinderToCell
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	p1_(dict.lookup("p1")),
	p2_(dict.lookup("p2")),
	radius_(readScalar(dict.lookup("radius")))
{}


// Construct from Istream
Foam::cylinderToCell::cylinderToCell
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	p1_(checkIs(is)),
	p2_(checkIs(is)),
	radius_(readScalar(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cylinderToCell::~cylinderToCell()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cylinderToCell::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding cells with centre within cylinder, with p1 = "
			<< p1_ << ", p2 = " << p2_ << " and radius = " << radius_ << endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing cells with centre within sphere, with p1 = "
			<< p1_ << ", p2 = " << p2_ << " and radius = " << radius_ << endl;

		combine(set, false);
	}
}


// ************************************************************************* //
