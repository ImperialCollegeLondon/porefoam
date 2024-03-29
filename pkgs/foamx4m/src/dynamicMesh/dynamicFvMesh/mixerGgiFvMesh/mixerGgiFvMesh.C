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

#include "mixerGgiFvMesh.H"
#include "foamTime.H"
#include "regionSplit.H"
#include "ggiPolyPatch.H"
#include "polyPatchID.H"
#include "addToRunTimeSelectionTable.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(mixerGgiFvMesh, 0);
	addToRunTimeSelectionTable(dynamicFvMesh, mixerGgiFvMesh, IOobject);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::mixerGgiFvMesh::addZonesAndModifiers()
{
	// Add zones and modifiers for motion action

	if (cellZones().size() > 0)
	{
		Info<< "void mixerGgiFvMesh::addZonesAndModifiers() : "
			<< "Zones and modifiers already present.  Skipping."
			<< endl;

		return;
	}

	Info<< "Time = " << time().timeName() << endl
		<< "Adding zones and modifiers to the mesh" << endl;

	// Add zones
	List<pointZone*> pz(0);
	List<faceZone*> fz(0);
	List<cellZone*> cz(1);

	// Copy the face zones associated with the GGI interfaces
	if (faceZones().size() > 0)
	{
		// Copy face zones
		Info << "Copying existing face zones" << endl;

		fz.setSize(faceZones().size());

		forAll (faceZones(), i)
		{
			fz[i] = faceZones()[i].clone(faceZones()).ptr();
		}
	}

	regionSplit rs(*this);

	// Get the region of the cell containing the origin.
	label originRegion = rs[findNearestCell(cs().origin())];

	labelList movingCells(nCells());
	label nMovingCells = 0;

	forAll(rs, cellI)
	{
		if (rs[cellI] == originRegion)
		{
			movingCells[nMovingCells] = cellI;
			nMovingCells++;
		}
	}

	movingCells.setSize(nMovingCells);
	Info << "Number of cells in the moving region: " << nMovingCells << endl;

	cz[0] = new cellZone
	(
		"movingCells",
		movingCells,
		0,
		cellZones()
	);

	Info << "Adding point, face and cell zones" << endl;
	removeZones();
	addZones(pz, fz, cz);

	// Write mesh
	syncUpdateMesh();
	write();
}


void Foam::mixerGgiFvMesh::calcMovingMasks() const
{
	if (debug)
	{
		Info<< "void mixerGgiFvMesh::calcMovingMasks() const : "
			<< "Calculating point and cell masks"
			<< endl;
	}

	if (movingPointsMaskPtr_)
	{
		FatalErrorIn("void mixerGgiFvMesh::calcMovingMasks() const")
			<< "point mask already calculated"
			<< abort(FatalError);
	}

	// Set the point mask
	movingPointsMaskPtr_ = new scalarField(allPoints().size(), 0);
	scalarField& movingPointsMask = *movingPointsMaskPtr_;

	const cellList& c = cells();
	const faceList& f = allFaces();

	label movingCellsID = cellZones().findZoneID("movingCells");

	if (movingCellsID < 0)
	{
		FatalErrorIn("void mixerGgiFvMesh::calcMovingMasks() const")
			<< "Cannot find moving cell zone ID"
			<< abort(FatalError);
	}

	const labelList& cellAddr = cellZones()[movingCellsID];

	forAll (cellAddr, cellI)
	{
		const cell& curCell = c[cellAddr[cellI]];

		forAll (curCell, faceI)
		{
			// Mark all the points as moving
			const face& curFace = f[curCell[faceI]];

			forAll (curFace, pointI)
			{
				movingPointsMask[curFace[pointI]] = 1;
			}
		}
	}

	// Grab the ggi patches on the moving side
	wordList movingPatches(dict_.subDict("slider").lookup("moving"));

	forAll (movingPatches, patchI)
	{
		const label movingSliderID =
			boundaryMesh().findPatchID(movingPatches[patchI]);

		if (movingSliderID < 0)
		{
			FatalErrorIn("void mixerGgiFvMesh::calcMovingMasks() const")
				<< "Moving slider named " << movingPatches[patchI]
				<< " not found.  Valid patch names: "
				<< boundaryMesh().names() << abort(FatalError);
		}

		const ggiPolyPatch& movingGgiPatch =
			refCast<const ggiPolyPatch>(boundaryMesh()[movingSliderID]);

		const labelList& movingSliderAddr = movingGgiPatch.zone();

		forAll (movingSliderAddr, faceI)
		{
			const face& curFace = f[movingSliderAddr[faceI]];

			forAll (curFace, pointI)
			{
				movingPointsMask[curFace[pointI]] = 1;
			}
		}
	}

	// Grab the ggi patches on the static side
	wordList staticPatches(dict_.subDict("slider").lookup("static"));

	forAll (staticPatches, patchI)
	{
		const label staticSliderID =
			boundaryMesh().findPatchID(staticPatches[patchI]);

		if (staticSliderID < 0)
		{
			FatalErrorIn("void mixerGgiFvMesh::calcMovingMasks() const")
				<< "Static slider named " << staticPatches[patchI]
				<< " not found.  Valid patch names: "
				<< boundaryMesh().names() << abort(FatalError);
		}

		const ggiPolyPatch& staticGgiPatch =
			refCast<const ggiPolyPatch>(boundaryMesh()[staticSliderID]);

		const labelList& staticSliderAddr = staticGgiPatch.zone();

		forAll (staticSliderAddr, faceI)
		{
			const face& curFace = f[staticSliderAddr[faceI]];

			forAll (curFace, pointI)
			{
				movingPointsMask[curFace[pointI]] = 0;
			}
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::mixerGgiFvMesh::mixerGgiFvMesh
(
	const IOobject& io
)
:
	dynamicFvMesh(io),
	dict_
	(
		IOdictionary
		(
			IOobject
			(
				"dynamicMeshDict",
				time().constant(),
				*this,
				IOobject::MUST_READ_IF_MODIFIED,
				IOobject::NO_WRITE
			)
		).subDict(typeName + "Coeffs")
	),
	cs_
	(
		"coordinateSystem",
		dict_.subDict("coordinateSystem")
	),
	rpm_(readScalar(dict_.lookup("rpm"))),
	movingPointsMaskPtr_(nullptr)
{
	// Make sure the coordinate system does not operate in degrees
	// Bug fix, HJ, 3/Oct/2011
	if (!cs_.inDegrees())
	{
		WarningIn("mixerGgiFvMesh::mixerGgiFvMesh(const IOobject& io)")
			<< "Mixer coordinate system is set to operate in radians.  "
			<< "Changing to rad for correct calculation of angular velocity."
			<< nl
			<< "To remove this message please add entry" << nl << nl
			<< "inDegrees true;" << nl << nl
			<< "to the specification of the coordinate system"
			<< endl;

		cs_.inDegrees() = true;
	}

	addZonesAndModifiers();

	Info<< "Mixer mesh:" << nl
		<< "    origin: " << cs().origin() << nl
		<< "    axis  : " << cs().axis() << nl
		<< "    rpm   : " << rpm_ << endl;
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::mixerGgiFvMesh::~mixerGgiFvMesh()
{
	deleteDemandDrivenData(movingPointsMaskPtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Return moving points mask.  Moving points marked with 1
const Foam::scalarField& Foam::mixerGgiFvMesh::movingPointsMask() const
{
	if (!movingPointsMaskPtr_)
	{
		calcMovingMasks();
	}

	return *movingPointsMaskPtr_;
}


bool Foam::mixerGgiFvMesh::update()
{
	// Rotational speed needs to be converted from rpm
	movePoints
	(
		cs_.globalPosition
		(
			cs_.localPosition(allPoints())
		  + vector(0, rpm_*360.0*time().deltaT().value()/60.0, 0)*
			movingPointsMask()
		)
	);

	// The mesh is not morphing, but flux re-calculation is required
	// HJ, 25/Jan/2016
	return true;
}


// ************************************************************************* //
