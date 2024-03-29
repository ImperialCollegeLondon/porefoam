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

Description
	Cell layer addition/removal mesh modifier

\*---------------------------------------------------------------------------*/

#include "layerAdditionRemoval.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "foamTime.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(layerAdditionRemoval, 0);
	addToRunTimeSelectionTable
	(
		polyMeshModifier,
		layerAdditionRemoval,
		dictionary
	);
}


const Foam::debug::tolerancesSwitch
Foam::layerAdditionRemoval::motionDelta_
(
	"layerAdditionRemoval::motionDelta",
	0.01
);

const Foam::debug::tolerancesSwitch
Foam::layerAdditionRemoval::addDelta_
(
	"layerAdditionRemoval::addDelta",
	0.3
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::layerAdditionRemoval::checkDefinition()
{
	if (!faceZoneID_.active())
	{
		FatalErrorIn
		(
			"void Foam::layerAdditionRemoval::checkDefinition()"
		)   << "Master face zone named " << faceZoneID_.name()
			<< " cannot be found."
			<< abort(FatalError);
	}

	const polyMesh& mesh = topoChanger().mesh();

	if (debug > 1)
	{
		fileName fvPath(mesh.time().path()/"VTK");
		mkDir(fvPath);

		Info<< "Writing VTK files master face zone"
			<< "Layer addition/removal " << name()
			<< ", face zone " << faceZoneID_.name()
			<< "into " << fvPath
			<< endl;

		primitiveFacePatch::writeVTK
		(
			fvPath/fileName(faceZoneID_.name() + "FaceZone"),
			mesh.faceZones()[faceZoneID_.index()]().localFaces(),
			mesh.faceZones()[faceZoneID_.index()]().localPoints()
		);

		primitiveFacePatch::writeVTKNormals
		(
			fvPath/fileName(faceZoneID_.name() + "FaceZoneNormals"),
			mesh.faceZones()[faceZoneID_.index()]().localFaces(),
			mesh.faceZones()[faceZoneID_.index()]().localPoints()
		);
	}

	if
	(
		minLayerThickness_ < VSMALL
	 || maxLayerThickness_ < minLayerThickness_
	)
	{
		FatalErrorIn
		(
			"void Foam::layerAdditionRemoval::checkDefinition()"
		)   << "Incorrect layer thickness definition."
			<< abort(FatalError);
	}

	// Check size of zones
	label globalZoneSize =
		returnReduce
		(
		   mesh.faceZones()[faceZoneID_.index()].size(),
			sumOp<label>()
		);

	if (globalZoneSize == 0)
	{
		FatalErrorIn
		(
			"void Foam::layerAdditionRemoval::checkDefinition()"
		)   << "Face extrusion zone contains no faces.  Please check your "
			<< "mesh definition."
			<< abort(FatalError);
	}

	if (debug)
	{
		Pout<< "Cell layer addition/removal object " << name() << " :" << nl
			<< "    faceZoneID: " << faceZoneID_ << endl;
	}
}

Foam::scalar Foam::layerAdditionRemoval::readOldThickness
(
	const dictionary& dict
)
{
	if (dict.found("oldLayerThickness"))
	{
		return readScalar(dict.lookup("oldLayerThickness"));
	}
	else
	{
		return -1.0;
	}
}


void Foam::layerAdditionRemoval::clearAddressing() const
{
	// Layer removal data
	deleteDemandDrivenData(pointsPairingPtr_);
	deleteDemandDrivenData(facesPairingPtr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::layerAdditionRemoval::layerAdditionRemoval
(
	const word& name,
	const label index,
	const polyTopoChanger& mme,
	const word& zoneName,
	const scalar minThickness,
	const scalar maxThickness,
	const label cellZone
)
:
	polyMeshModifier(name, index, mme, true),
	faceZoneID_(zoneName, mme.mesh().faceZones()),
	minLayerThickness_(minThickness),
	maxLayerThickness_(maxThickness),
	oldLayerThickness_(-1.0),
	pointsPairingPtr_(nullptr),
	facesPairingPtr_(nullptr),
	triggerRemoval_(-1),
	triggerAddition_(-1),
	cellZone_(cellZone)
{
	checkDefinition();
}


// Construct from dictionary
Foam::layerAdditionRemoval::layerAdditionRemoval
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyTopoChanger& mme
)
:
	polyMeshModifier(name, index, mme, Switch(dict.lookup("active"))),
	faceZoneID_(dict.lookup("faceZoneName"), mme.mesh().faceZones()),
	minLayerThickness_(readScalar(dict.lookup("minLayerThickness"))),
	maxLayerThickness_(readScalar(dict.lookup("maxLayerThickness"))),
	oldLayerThickness_(readOldThickness(dict)),
	pointsPairingPtr_(nullptr),
	facesPairingPtr_(nullptr),
	triggerRemoval_(-1),
	triggerAddition_(-1),
	cellZone_(-1)
{
	checkDefinition();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::layerAdditionRemoval::~layerAdditionRemoval()
{
	clearAddressing();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void  Foam::layerAdditionRemoval::setRemoval()
{
	triggerRemoval_ = topoChanger().morphIndex();
}


void  Foam::layerAdditionRemoval::setAddition()
{
	triggerAddition_ = topoChanger().morphIndex();
}


bool Foam::layerAdditionRemoval::changeTopology() const
{
	// Protect from multiple calculation in the same time-step
	if (triggerRemoval_ > -1 || triggerAddition_ > -1)
	{
		return true;
	}

	// Go through all the cells in the master layer and calculate
	// approximate layer thickness as the ratio of the cell volume and
	// face area in the face zone.
	// Layer addition:
	//	 When the max thickness exceeds the threshold, trigger refinement.
	// Layer removal:
	//	 When the min thickness falls below the threshold, trigger removal.

	const faceZone& fz = topoChanger().mesh().faceZones()[faceZoneID_.index()];
	const labelList& mc = fz.masterCells();

	if (debug)
	{
		// Check master cell addressing
		if (min(mc) < 0)
		{
			const polyMesh& mesh = topoChanger().mesh();

			fileName fvPath(mesh.time().path()/"VTK");
			mkDir(fvPath);

			Info<< "Writing VTK files master face zone"
				<< "Layer addition/removal " << name()
				<< ", face zone " << faceZoneID_.name()
				<< "into " << fvPath
				<< endl;

			primitiveFacePatch::writeVTK
			(
				fvPath/fileName(faceZoneID_.name() + "FaceZone"),
				mesh.faceZones()[faceZoneID_.index()]().localFaces(),
				mesh.faceZones()[faceZoneID_.index()]().localPoints()
			);

			primitiveFacePatch::writeVTKNormals
			(
				fvPath/fileName(faceZoneID_.name() + "FaceZoneNormals"),
				mesh.faceZones()[faceZoneID_.index()]().localFaces(),
				mesh.faceZones()[faceZoneID_.index()]().localPoints()
			);

			FatalErrorIn("bool layerAdditionRemoval::changeTopology() const")
				<< "Error in master cell addressing for face zone "
					<< faceZoneID_.name()
					<< abort(FatalError);
		}
	}

	const scalarField& V = topoChanger().mesh().cellVolumes();
	const vectorField& S = topoChanger().mesh().faceAreas();

	if (min(V) < -VSMALL)
	{
		FatalErrorIn("bool layerAdditionRemoval::changeTopology() const")
			<< "negative cell volume. Error in mesh motion before "
			<< "topological change."
			<< abort(FatalError);
	}

	scalar minDelta = GREAT;
	scalar maxDelta = 0;
	scalar avgDelta = 0;
	scalar nAvg = 0;

	if (fz.size())
	{
		forAll (fz, faceI)
		{
			scalar curDelta = V[mc[faceI]]/mag(S[fz[faceI]]);
			avgDelta += curDelta;
			minDelta = min(minDelta, curDelta);
			maxDelta = max(maxDelta, curDelta);
		}

		nAvg += fz.size();
	}

	// If the patch is empty on a processor in a parallel simulation,
	// original values will be preserved.  HJ, 7/Mar/2011

	reduce(minDelta, minOp<scalar>());
	reduce(maxDelta, maxOp<scalar>());
	reduce(avgDelta, sumOp<scalar>());
	reduce(nAvg, sumOp<scalar>());

	avgDelta /= nAvg;

	if (debug)
	{
		Info<< "bool layerAdditionRemoval::changeTopology() const "
			<< " for object " << name() << " : " << nl
			<< "Layer thickness: min: " << minDelta
			<< " max: " << maxDelta << " avg: " << avgDelta
			<< " old thickness: " << oldLayerThickness_ << nl
			<< "Removal threshold: " << minLayerThickness_
			<< " addition threshold: " << maxLayerThickness_ << endl;
	}

	bool topologicalChange = false;

	// If the thickness is decreasing and crosses the min thickness,
	// trigger removal
	if (oldLayerThickness_ < 0)
	{
		if (debug)
		{
			Info<< "First step. No addition/removal" << endl;
		}

		// No topological changes allowed before first mesh motion
		// HJ, 8/Oct/2002
		oldLayerThickness_ = avgDelta;

		topologicalChange = false;
	}
	// New criterion to avoid round-off triggering layer addition/removal
	// HJ, 30/Mar/2018
	else if
	(
		(oldLayerThickness_ - avgDelta) > motionDelta_()*minLayerThickness_
	)
	{
		// Layers moving towards removal
		if (minDelta < minLayerThickness_)
		{
			// Check layer pairing
			if (setLayerPairing())
			{
				// A mesh layer detected.  Check that collapse is valid
				if (validCollapse())
				{
					// At this point, info about moving the old mesh
					// in a way to collapse the cells in the removed
					// layer is available.  Not sure what to do with
					// it.  HJ, 3/Nov/2003

					if (debug)
					{
						Info<< "bool layerAdditionRemoval::changeTopology() "
							<< " const for object " << name() << " : "
							<< "Triggering layer removal" << endl;
					}

					triggerRemoval_ = topoChanger().morphIndex();

					// Old thickness looses meaning.
					// Set it up to indicate layer removal
					oldLayerThickness_ = GREAT;

					topologicalChange = true;
				}
				else
				{
					// No removal, clear addressing
					clearAddressing();
				}
			}
		}
		else
		{
			oldLayerThickness_ = avgDelta;
		}
	}
	// New criterion to avoid round-off triggering layer addition/removal
	// HJ, 30/Mar/2018
	else if
	(
		(avgDelta - oldLayerThickness_) > motionDelta_()*minLayerThickness_
	)
	{
		// Layers moving towards addition
		if (maxDelta > maxLayerThickness_)
		{
			if (debug)
			{
				Info<< "bool layerAdditionRemoval::changeTopology() const "
					<< " for object " << name() << " : "
					<< "Triggering layer addition" << endl;
			}

			triggerAddition_ = topoChanger().morphIndex();

			// Old thickness looses meaning.
			// Set it up to indicate layer removal
			oldLayerThickness_ = 0;

			topologicalChange = true;
		}
		else
		{
			oldLayerThickness_ = avgDelta;
		}
	}
	// else the motion change is smaller than the tolerance and the layer
	// interface is practically static.  HJ, 30/Mar/2018

	return topologicalChange;
}


void Foam::layerAdditionRemoval::setRefinement(polyTopoChange& ref) const
{
	// Insert the layer addition/removal instructions
	// into the topological change

	if (triggerRemoval_ == topoChanger().morphIndex())
	{
		removeCellLayer(ref);

		// Clear addressing.  This also resets the addition/removal data
		if (debug)
		{
			Info<< "layerAdditionRemoval::setRefinement(polyTopoChange& ref) "
				<< " for object " << name() << " : "
				<< "Clearing addressing after layer removal. " << endl;
		}

		triggerRemoval_ = -1;
		clearAddressing();
	}

	if (triggerAddition_ == topoChanger().morphIndex())
	{
		addCellLayer(ref);

		// Clear addressing.  This also resets the addition/removal data
		if (debug)
		{
			Info<< "layerAdditionRemoval::setRefinement(polyTopoChange& ref) "
				<< " for object " << name() << " : "
				<< "Clearing addressing after layer addition. " << endl;
		}

		triggerAddition_ = -1;
		clearAddressing();
	}
}


void Foam::layerAdditionRemoval::updateMesh(const mapPolyMesh&)
{
	if (debug)
	{
		Info<< "layerAdditionRemoval::updateMesh(const mapPolyMesh&) "
			<< " for object " << name() << " : "
			<< "Clearing addressing on external request. ";

		if (pointsPairingPtr_ || facesPairingPtr_)
		{
			Info << "Pointers set." << endl;
		}
		else
		{
			Info << "Pointers not set." << endl;
		}
	}

	// Mesh has changed topologically.  Update local topological data
	faceZoneID_.update(topoChanger().mesh().faceZones());

	clearAddressing();
}


void Foam::layerAdditionRemoval::setMinLayerThickness(const scalar t) const
{
	if
	(
		t < VSMALL
	 || maxLayerThickness_ < t
	)
	{
		FatalErrorIn
		(
			"void layerAdditionRemoval::setMinLayerThickness("
			"const scalar t) const"
		)   << "Incorrect layer thickness definition."
			<< abort(FatalError);
	}

	minLayerThickness_ = t;
}


void Foam::layerAdditionRemoval::setMaxLayerThickness(const scalar t) const
{
	if (t < minLayerThickness_)
	{
		FatalErrorIn
		(
			"void layerAdditionRemoval::setMaxLayerThickness("
			"const scalar t) const"
		)   << "Incorrect layer thickness definition."
			<< abort(FatalError);
	}

	maxLayerThickness_ = t;
}


void Foam::layerAdditionRemoval::write(Ostream& os) const
{
	os  << nl << type() << nl
		<< name()<< nl
		<< faceZoneID_ << nl
		<< minLayerThickness_ << nl
		<< oldLayerThickness_ << nl
		<< maxLayerThickness_ << endl;
}


void Foam::layerAdditionRemoval::writeDict(Ostream& os) const
{
	os  << nl << name() << nl << token::BEGIN_BLOCK << nl
		<< "    type " << type()
		<< token::END_STATEMENT << nl
		<< "    faceZoneName " << faceZoneID_.name()
		<< token::END_STATEMENT << nl
		<< "    minLayerThickness " << minLayerThickness_
		<< token::END_STATEMENT << nl
		<< "    maxLayerThickness " << maxLayerThickness_
		<< token::END_STATEMENT << nl
		<< "    oldLayerThickness " << oldLayerThickness_
		<< token::END_STATEMENT << nl
		<< "    active " << active()
		<< token::END_STATEMENT << nl
		<< token::END_BLOCK << endl;
}


// ************************************************************************* //
