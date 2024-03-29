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

#include "cellToFace.H"
#include "polyMesh.H"
#include "cellSet.H"
#include "foamTime.H"
#include "syncTools.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(cellToFace, 0);

addToRunTimeSelectionTable(topoSetSource, cellToFace, word);

addToRunTimeSelectionTable(topoSetSource, cellToFace, istream);

}


Foam::topoSetSource::addToUsageTable Foam::cellToFace::usage_
(
	cellToFace::typeName,
	"\n    Usage: cellToFace <cellSet> all|both\n\n"
	"    Select -all : all faces of cells in the cellSet\n"
	"           -both: faces where both neighbours are in the cellSet\n\n"
);

template<>
const char* Foam::NamedEnum<Foam::cellToFace::cellAction, 2>::names[] =
{
	"all",
	"both"
};

const Foam::NamedEnum<Foam::cellToFace::cellAction, 2>
	Foam::cellToFace::cellActionNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::cellToFace::combine(topoSet& set, const bool add) const
{
	// Load the set
	if (!exists(mesh_.time().path()/topoSet::localPath(mesh_, setName_)))
	{
		SeriousError<< "Cannot load set "
			<< setName_ << endl;
	}

	cellSet loadedSet(mesh_, setName_);

	if (option_ == ALL)
	{
		// Add all faces from cell
		for
		(
			cellSet::const_iterator iter = loadedSet.begin();
			iter != loadedSet.end();
			++iter
		)
		{
			label cellI = iter.key();

			const labelList& cFaces = mesh_.cells()[cellI];

			forAll(cFaces, cFaceI)
			{
				addOrDelete(set, cFaces[cFaceI], add);
			}
		}
	}
	else if (option_ == BOTH)
	{
		// Add all faces whose both neighbours are in set.

		label nInt = mesh_.nInternalFaces();
		const labelList& own = mesh_.faceOwner();
		const labelList& nei = mesh_.faceNeighbour();
		const polyBoundaryMesh& patches = mesh_.boundaryMesh();


		// Check all internal faces
		for (label faceI = 0; faceI < nInt; faceI++)
		{
			if (loadedSet.found(own[faceI]) && loadedSet.found(nei[faceI]))
			{
				addOrDelete(set, faceI, add);
			}
		}


		// Get coupled cell status
		boolList neiInSet(mesh_.nFaces()-nInt, false);

		forAll(patches, patchI)
		{
			const polyPatch& pp = patches[patchI];

			if (pp.coupled())
			{
				label faceI = pp.start();
				forAll(pp, i)
				{
					neiInSet[faceI-nInt] = loadedSet.found(own[faceI]);
					faceI++;
				}
			}
		}
		syncTools::swapBoundaryFaceList(mesh_, neiInSet, false);


		// Check all boundary faces
		forAll(patches, patchI)
		{
			const polyPatch& pp = patches[patchI];

			if (pp.coupled())
			{
				label faceI = pp.start();
				forAll(pp, i)
				{
					if (loadedSet.found(own[faceI]) && neiInSet[faceI-nInt])
					{
					    addOrDelete(set, faceI, add);
					}
					faceI++;
				}
			}
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from componenta
Foam::cellToFace::cellToFace
(
	const polyMesh& mesh,
	const word& setName,
	const cellAction option
)
:
	topoSetSource(mesh),
	setName_(setName),
	option_(option)
{}


// Construct from dictionary
Foam::cellToFace::cellToFace
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	setName_(dict.lookup("set")),
	option_(cellActionNames_.read(dict.lookup("option")))
{}


// Construct from Istream
Foam::cellToFace::cellToFace
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	setName_(checkIs(is)),
	option_(cellActionNames_.read(checkIs(is)))
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::cellToFace::~cellToFace()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::cellToFace::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ((action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding faces according to cellSet " << setName_
			<< " ..." << endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing faces according to cellSet " << setName_
			<< " ..." << endl;

		combine(set, false);
	}
}


// ************************************************************************* //
