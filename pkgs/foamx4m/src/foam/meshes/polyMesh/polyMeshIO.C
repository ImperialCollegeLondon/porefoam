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

#include "polyMesh.H"
#include "foamTime.H"
#include "IOstream.H"
#include "cellIOList.H"
#include "meshObjectBase.H"
#include "mapPolyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::polyMesh::polyMesh
(
	const IOobject& io,
	Istream& is,
	const bool syncPar
)
:
	objectRegistry(io),
	primitiveMesh(),
	allPoints_
	(
		IOobject
		(
			"points",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		pointField(is)
	),
	// To be re-sliced later.  HJ, 19/Oct/2008
	points_(allPoints_, allPoints_.size()),
	allFaces_
	(
		IOobject
		(
			"faces",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		faceList(is)
	),
	// To be re-sliced later.  HJ, 19/Oct/2008
	faces_(allFaces_, allFaces_.size()),
	owner_
	(
		IOobject
		(
			"owner",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		labelList(is)
	),
	neighbour_
	(
		IOobject
		(
			"neighbour",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		labelList(is)
	),
	syncPar_(syncPar),
	clearedPrimitives_(false),
	boundary_
	(
		IOobject
		(
			"boundary",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		*this,
		is
	),
	bounds_(allPoints_, syncPar),
	geometricD_(Vector<label>::zero),
	solutionD_(Vector<label>::zero),
	comm_(Pstream::worldComm),
	pointZones_
	(
		IOobject
		(
			"pointZones",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		*this,
		is
	),
	faceZones_
	(
		IOobject
		(
			"faceZones",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		*this,
		is
	),
	cellZones_
	(
		IOobject
		(
			"cellZones",
			instance(),
			meshSubDir,
			*this,
			IOobject::NO_READ,
			IOobject::NO_WRITE
		),
		*this,
		is
	),
	globalMeshDataPtr_(nullptr),
	moving_(false),
	changing_(false),
	curMotionTimeIndex_(time().timeIndex()),
	oldAllPointsPtr_(nullptr),
	oldPointsPtr_(nullptr)
{
	// Check if the faces and cells are valid
	forAll (allFaces_, faceI)
	{
		const face& curFace = allFaces_[faceI];

		if (min(curFace) < 0 || max(curFace) > allPoints_.size())
		{
			FatalErrorIn
			(
				"polyMesh::polyMesh\n"
				"(\n"
				"    const IOobject& io,\n"
				"    const pointField& points,\n"
				"    const faceList& faces,\n"
				"    const cellList& cells\n"
				")\n"
			)   << "Face " << faceI << "contains vertex labels out of range: "
				<< curFace << " Max point index = " << allPoints_.size()
				<< abort(FatalError);
		}
	}

	// Set the primitive mesh
	initMesh();
}



// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::polyMesh::setInstance(const fileName& inst)
{
	if (debug)
	{
		Info<< "void polyMesh::setInstance(const fileName& inst) : "
			<< "Resetting file instance to " << inst << endl;
	}

	allPoints_.instance() = inst;
	allFaces_.instance() = inst;
	owner_.instance() = inst;
	neighbour_.instance() = inst;
	boundary_.instance() = inst;
	pointZones_.instance() = inst;
	faceZones_.instance() = inst;
	cellZones_.instance() = inst;

	setMotionWriteOpt(IOobject::AUTO_WRITE);
	setTopoWriteOpt(IOobject::AUTO_WRITE);
}


void Foam::polyMesh::setMotionWriteOpt(IOobject::writeOption wOpt)
{
	if (debug)
	{
		Info<< "void polyMesh::setMotionWriteOpt(IOobject::writeOption) "
			<< "Setting motion writeOpt to " << wOpt << endl;
	}

	allPoints_.writeOpt() = wOpt;
}


void Foam::polyMesh::setTopoWriteOpt(IOobject::writeOption wOpt)
{
	if (debug)
	{
		Info<< "void polyMesh::setTopoWriteOpt(IOobject::writeOption) "
			<< "Setting topo writeOpt to " << wOpt << endl;
	}

	allFaces_.writeOpt() = wOpt;
	owner_.writeOpt() = wOpt;
	neighbour_.writeOpt() = wOpt;
	boundary_.writeOpt() = wOpt;

	pointZones_.writeOpt() = wOpt;
	faceZones_.writeOpt() = wOpt;
	cellZones_.writeOpt() = wOpt;

}


Foam::polyMesh::readUpdateState Foam::polyMesh::readUpdate()
{
	if (debug)
	{
		Info<< "polyMesh::readUpdateState polyMesh::readUpdate() : "
			<< "Updating mesh based on saved data." << endl;
	}

	// Find the point and cell instance
	fileName pointsInst(time().findInstance(meshDir(), "points"));
	fileName facesInst(time().findInstance(meshDir(), "faces"));

	if (debug)
	{
		Info<< "Faces instance: old = " << facesInstance()
			<< " new = " << facesInst << nl
			<< "Points instance: old = " << pointsInstance()
			<< " new = " << pointsInst << endl;
	}

	if (facesInst != facesInstance())
	{
		// Topological change
		if (debug)
		{
			Info << "Topological change" << endl;
		}

		clearOut();

		// Set instance to new instance. Note that points instance can differ
		// from from faces instance.
		setInstance(facesInst);
		allPoints_.instance() = pointsInst;

		allPoints_ = pointIOField
		(
			IOobject
			(
				"points",
				pointsInst,
				meshSubDir,
				*this,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		allFaces_ = faceIOList
		(
			IOobject
			(
				"faces",
				facesInst,
				meshSubDir,
				*this,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		owner_ = labelIOList
		(
			IOobject
			(
				"owner",
				facesInst,
				meshSubDir,
				*this,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE,
				false
			)
		);

		neighbour_ = labelIOList
		(
			IOobject
			(
				"neighbour",
				facesInst,
				meshSubDir,
				*this,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE,
				false
			)
		);

		// Reset the boundary patches
		polyBoundaryMesh newBoundary
		(
			IOobject
			(
				"boundary",
				facesInst,
				meshSubDir,
				*this,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			),
			*this
		);

		// Check that patch types and names are unchanged
		bool boundaryChanged = false;

		if (newBoundary.size() != boundary_.size())
		{
			boundaryChanged = true;
		}
		else
		{
			wordList newTypes = newBoundary.types();
			wordList newNames = newBoundary.names();

			wordList oldTypes = boundary_.types();
			wordList oldNames = boundary_.names();

			forAll (oldTypes, patchI)
			{
				if
				(
					oldTypes[patchI] != newTypes[patchI]
				 || oldNames[patchI] != newNames[patchI]
				)
				{
					boundaryChanged = true;
					break;
				}
			}
		}

		if (boundaryChanged)
		{
			WarningIn("polyMesh::readUpdateState polyMesh::readUpdate()")
				<< "Number of patches has changed.  This may have "
				<< "unexpected consequences.  Proceed with care." << endl;

			boundary_.clear();
			boundary_.setSize(newBoundary.size());

			forAll (newBoundary, patchI)
			{
				boundary_.set(patchI, newBoundary[patchI].clone(boundary_));
			}
		}
		else
		{
			forAll (boundary_, patchI)
			{
				boundary_[patchI] = polyPatch
				(
					newBoundary[patchI].name(),
					newBoundary[patchI].size(),
					newBoundary[patchI].start(),
					patchI,
					boundary_
				);
			}
		}


		// Boundary is set so can use initMesh now (uses boundary_ to
		// determine internal and active faces)

		if (exists(owner_.objectPath()))
		{
			initMesh();
		}
		else
		{
			cellIOList cells
			(
				IOobject
				(
					"cells",
					facesInst,
					meshSubDir,
					*this,
					IOobject::MUST_READ,
					IOobject::NO_WRITE,
					false
				)
			);

			// Recalculate the owner/neighbour addressing and reset the
			// primitiveMesh
			initMesh(cells);
		}

		// Zones
		pointZoneMesh newPointZones
		(
			IOobject
			(
				"pointZones",
				facesInst,
				meshSubDir,
				*this,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE,
				false
			),
			*this
		);

		label oldSize = pointZones_.size();

		if (newPointZones.size() <= pointZones_.size())
		{
			pointZones_.setSize(newPointZones.size());
		}

		// Reset existing ones
		forAll (pointZones_, czI)
		{
			pointZones_[czI] = newPointZones[czI];
		}

		// Extend with extra ones
		pointZones_.setSize(newPointZones.size());

		for (label czI = oldSize; czI < newPointZones.size(); czI++)
		{
			pointZones_.set(czI, newPointZones[czI].clone(pointZones_));
		}

		pointZones_.setSize(newPointZones.size());
		forAll (pointZones_, pzI)
		{
			pointZones_[pzI] = newPointZones[pzI];
		}


		faceZoneMesh newFaceZones
		(
			IOobject
			(
				"faceZones",
				facesInst,
				meshSubDir,
				*this,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE,
				false
			),
			*this
		);

		oldSize = faceZones_.size();

		if (newFaceZones.size() <= faceZones_.size())
		{
			faceZones_.setSize(newFaceZones.size());
		}

		// Reset existing ones
		forAll (faceZones_, fzI)
		{
			faceZones_[fzI].resetAddressing
			(
				newFaceZones[fzI],
				newFaceZones[fzI].flipMap()
			);
		}

		// Extend with extra ones
		faceZones_.setSize(newFaceZones.size());

		for (label fzI = oldSize; fzI < newFaceZones.size(); fzI++)
		{
			faceZones_.set(fzI, newFaceZones[fzI].clone(faceZones_));
		}


		cellZoneMesh newCellZones
		(
			IOobject
			(
				"cellZones",
				facesInst,
				meshSubDir,
				*this,
				IOobject::READ_IF_PRESENT,
				IOobject::NO_WRITE,
				false
			),
			*this
		);

		oldSize = cellZones_.size();

		if (newCellZones.size() <= cellZones_.size())
		{
			cellZones_.setSize(newCellZones.size());
		}

		// Reset existing ones
		forAll (cellZones_, czI)
		{
			cellZones_[czI] = newCellZones[czI];
		}

		// Extend with extra ones
		cellZones_.setSize(newCellZones.size());

		for (label czI = oldSize; czI < newCellZones.size(); czI++)
		{
			cellZones_.set(czI, newCellZones[czI].clone(cellZones_));
		}

		// Even if number of patches stayed same still recalculate boundary
		// data.

		// Calculate topology for the patches (processor-processor comms etc.)
		boundary_.updateMesh();

		// Calculate the geometry for the patches (transformation tensors etc.)
		boundary_.calcGeometry();

		// Derived info
		bounds_ = boundBox(allPoints_);
		geometricD_ = Vector<label>::zero;
		solutionD_ = Vector<label>::zero;

		// Instantiate a dummy mapPolyMesh
		autoPtr<mapPolyMesh> mapPtr(new mapPolyMesh(*this));

		// Execute dummy  topo change on all mesh objects
		meshObjectBase::allUpdateTopology(*this, mapPtr());

		if (boundaryChanged)
		{
			return polyMesh::TOPO_PATCH_CHANGE;
		}
		else
		{
			return polyMesh::TOPO_CHANGE;
		}
	}
	else if (pointsInst != pointsInstance())
	{
		// Points moved
		if (debug)
		{
			Info << "Point motion" << endl;
		}

		clearGeom();

		allPoints_.instance() = pointsInst;

		allPoints_ = pointIOField
		(
			IOobject
			(
				"points",
				pointsInst,
				meshSubDir,
				*this,
				IOobject::MUST_READ,
				IOobject::NO_WRITE,
				false
			)
		);

		// Reset points, mesh is not moved
		points_ = pointField::subField(allPoints_, nPoints());

		// Derived info
		bounds_ = boundBox(allPoints_);

		// Rotation can cause direction vector to change
		geometricD_ = Vector<label>::zero;
		solutionD_ = Vector<label>::zero;

		// Move points in all mesh objects
		meshObjectBase::allMovePoints<polyMesh>(*this);

		return polyMesh::POINTS_MOVED;
	}
	else
	{
		if (debug)
		{
			Info << "No change" << endl;
		}

		return polyMesh::UNCHANGED;
	}
}


// * * * * * * * * * * * * IOstream operators  * * * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const polyMesh& m)
{
	os  << m.allPoints() << nl
		<< m.allFaces() << nl
		<< m.faceOwner() << nl
		<< m.faceNeighbour() << nl
		<< m.boundaryMesh() << nl
		<< m.pointZones() << nl
		<< m.faceZones() << nl
		<< m.cellZones() << nl
		<< endl;

	return os;
}


// ************************************************************************* //
