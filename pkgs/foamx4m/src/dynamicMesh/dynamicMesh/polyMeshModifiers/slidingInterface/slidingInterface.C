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


Class
	slidingInterface

Description
	Sliding interface

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "slidingInterface.H"
#include "polyTopoChanger.H"
#include "polyMesh.H"
#include "primitiveMesh.H"
#include "polyTopoChange.H"
#include "mapPolyMesh.H"
#include "addToRunTimeSelectionTable.H"
#include "triPointRef.H"
#include "plane.H"

// Index of debug signs:
// p - adjusting a projection point
// * - adjusting edge intersection

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(slidingInterface, 0);
	addToRunTimeSelectionTable
	(
		polyMeshModifier,
		slidingInterface,
		dictionary
	);
}


template<>
const char* Foam::NamedEnum<Foam::slidingInterface::typeOfMatch, 2>::names[] =
{
	"integral",
	"partial"
};


const Foam::NamedEnum<Foam::slidingInterface::typeOfMatch, 2>
Foam::slidingInterface::typeOfMatchNames_;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::slidingInterface::checkDefinition() const
{
	const polyMesh& mesh = topoChanger().mesh();

	if
	(
		!masterFaceZoneID_.active()
	 || !slaveFaceZoneID_.active()
	 || !cutPointZoneID_.active()
	 || !cutFaceZoneID_.active()
	 || !masterPatchID_.active()
	 || !slavePatchID_.active()
	)
	{
		FatalErrorIn
		(
			"void slidingInterface::checkDefinition()"
		)   << "Sliding interface object " << name() << " :" << nl
			<< "    master face zone: " << masterFaceZoneID_.index() << nl
			<< "    slave face zone: " << slaveFaceZoneID_.index() << nl
			<< "Not all zones and patches needed in the definition "
			<< "have been found.  Please check your mesh definition." << nl
			<< "Error code: "
			<< masterFaceZoneID_.active() << slaveFaceZoneID_.active()
			<< cutPointZoneID_.active() << cutFaceZoneID_.active()
			<< masterPatchID_.active() << slavePatchID_.active()
			<< abort(FatalError);
	}

	// Check the sizes and set up state
	if
	(
		mesh.faceZones()[masterFaceZoneID_.index()].size() == 0
	 || mesh.faceZones()[slaveFaceZoneID_.index()].size() == 0
	)
	{
		FatalErrorIn("void slidingInterface::checkDefinition()")
			<< "Sliding interface object " << name() << " :" << nl
			<< "    master face zone: " << masterFaceZoneID_.index() << nl
			<< "    slave face zone: " << slaveFaceZoneID_.index() << nl
			<< "Master or slave face zone contain no faces "
			<< "Please check your mesh definition."
			<< abort(FatalError);
	}

	if (debug)
	{
		if (!attached_)
		{
			// Check for points shared between master and slave.  If a point
			// is shared, the projection will be illegal
			const primitiveFacePatch& masterPatch =
				mesh.faceZones()[masterFaceZoneID_.index()]();

			const primitiveFacePatch& slavePatch =
				mesh.faceZones()[slaveFaceZoneID_.index()]();

			const labelList& smp = slavePatch.meshPoints();
			const pointField& slp = slavePatch.localPoints();

			label nSharedPoints = 0;

			forAll (smp, i)
			{
				if (masterPatch.whichPoint(smp[i]) != -1)
				{
					Warning<< "Shared point between master and slave: "
						<< smp[i] << " " << slp[i] << endl;

					nSharedPoints++;
				}
			}

			if (nSharedPoints > 0)
			{
				FatalErrorIn("void slidingInterface::checkDefinition()")
					<< "Sliding interface object " << name() << " :" << nl
					<< "    master face zone: "
					<< masterFaceZoneID_.index() << nl
					<< "    slave face zone: "
					<< slaveFaceZoneID_.index() << nl
					<< "Master and slave face zone share " << nSharedPoints
					<< " point.  This is not allowed." << nl
					<< "Please check mesh for topological errors."
					<< abort(FatalError);
			}
		}

		Pout<< "Sliding interface object " << name() << " :" << nl
			<< "    master face zone: " << masterFaceZoneID_.index() << nl
			<< "    slave face zone: " << slaveFaceZoneID_.index() << endl;
	}
}


void Foam::slidingInterface::clearOut() const
{
	clearPointProjection();
	clearAttachedAddressing();
	clearAddressing();
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


// Construct from components
Foam::slidingInterface::slidingInterface
(
	const word& name,
	const label index,
	const polyTopoChanger& mme,
	const word& masterFaceZoneName,
	const word& slaveFaceZoneName,
	const word& cutPointZoneName,
	const word& cutFaceZoneName,
	const word& masterPatchName,
	const word& slavePatchName,
	const typeOfMatch tom,
	const bool coupleDecouple,
	const intersection::algorithm algo
)
:
	polyMeshModifier(name, index, mme, true),
	masterFaceZoneID_
	(
		masterFaceZoneName,
		mme.mesh().faceZones()
	),
	slaveFaceZoneID_
	(
		slaveFaceZoneName,
		mme.mesh().faceZones()
	),
	cutPointZoneID_
	(
		cutPointZoneName,
		mme.mesh().pointZones()
	),
	cutFaceZoneID_
	(
		cutFaceZoneName,
		mme.mesh().faceZones()
	),
	masterPatchID_
	(
		masterPatchName,
		mme.mesh().boundaryMesh()
	),
	slavePatchID_
	(
		slavePatchName,
		mme.mesh().boundaryMesh()
	),
	matchType_(tom),
	coupleDecouple_(coupleDecouple),
	attached_(false),
	projectionAlgo_(algo),
	trigger_(false),
	cutFaceMasterPtr_(nullptr),
	cutFaceSlavePtr_(nullptr),
	masterFaceCellsPtr_(nullptr),
	slaveFaceCellsPtr_(nullptr),
	masterStickOutFacesPtr_(nullptr),
	slaveStickOutFacesPtr_(nullptr),
	retiredPointMapPtr_(nullptr),
	cutPointEdgePairMapPtr_(nullptr),
	slavePointPointHitsPtr_(nullptr),
	slavePointEdgeHitsPtr_(nullptr),
	slavePointFaceHitsPtr_(nullptr),
	masterPointEdgeHitsPtr_(nullptr),
	projectedSlavePointsPtr_(nullptr)
{
	checkDefinition();

	if (attached_)
	{
		FatalErrorIn
		(
			"Foam::slidingInterface::slidingInterface\n"
			"(\n"
			"    const word& name,\n"
			"    const label index,\n"
			"    const polyTopoChanger& mme,\n"
			"    const word& masterFaceZoneName,\n"
			"    const word& slaveFaceZoneName,\n"
			"    const word& cutFaceZoneName,\n"
			"    const word& cutPointZoneName,\n"
			"    const word& masterPatchName,\n"
			"    const word& slavePatchName,\n"
			"    const typeOfMatch tom,\n"
			"    const bool coupleDecouple\n"
			")"
		)   << "Creation of a sliding interface from components "
			<< "in attached state not supported."
			<< abort(FatalError);
	}
	else
	{
		calcAttachedAddressing();
	}
}


// Construct from components
Foam::slidingInterface::slidingInterface
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyTopoChanger& mme
)
:
	polyMeshModifier(name, index, mme, Switch(dict.lookup("active"))),
	masterFaceZoneID_
	(
		dict.lookup("masterFaceZoneName"),
		mme.mesh().faceZones()
	),
	slaveFaceZoneID_
	(
		dict.lookup("slaveFaceZoneName"),
		mme.mesh().faceZones()
	),
	cutPointZoneID_
	(
		dict.lookup("cutPointZoneName"),
		mme.mesh().pointZones()
	),
	cutFaceZoneID_
	(
		dict.lookup("cutFaceZoneName"),
		mme.mesh().faceZones()
	),
	masterPatchID_
	(
		dict.lookup("masterPatchName"),
		mme.mesh().boundaryMesh()
	),
	slavePatchID_
	(
		dict.lookup("slavePatchName"),
		mme.mesh().boundaryMesh()
	),
	matchType_(typeOfMatchNames_.read((dict.lookup("typeOfMatch")))),
	coupleDecouple_(dict.lookup("coupleDecouple")),
	attached_(dict.lookup("attached")),
	projectionAlgo_
	(
		intersection::algorithmNames_.read(dict.lookup("projection"))
	),
	trigger_(false),
	cutFaceMasterPtr_(nullptr),
	cutFaceSlavePtr_(nullptr),
	masterFaceCellsPtr_(nullptr),
	slaveFaceCellsPtr_(nullptr),
	masterStickOutFacesPtr_(nullptr),
	slaveStickOutFacesPtr_(nullptr),
	retiredPointMapPtr_(nullptr),
	cutPointEdgePairMapPtr_(nullptr),
	slavePointPointHitsPtr_(nullptr),
	slavePointEdgeHitsPtr_(nullptr),
	slavePointFaceHitsPtr_(nullptr),
	masterPointEdgeHitsPtr_(nullptr),
	projectedSlavePointsPtr_(nullptr)
{
	checkDefinition();

	// If the interface is attached, the master and slave face zone addressing
	// needs to be read in; otherwise it will be created
	if (attached_)
	{
		if (debug)
		{
			Pout<< "slidingInterface::slidingInterface(...) "
				<< " for object " << name << " : "
				<< "Interface attached.  Reading master and slave face zones "
				<< "and retired point lookup." << endl;
		}

		// The face zone addressing is written out in the definition dictionary
		masterFaceCellsPtr_ = new labelList(dict.lookup("masterFaceCells"));
		slaveFaceCellsPtr_ = new labelList(dict.lookup("slaveFaceCells"));

		masterStickOutFacesPtr_ =
			new labelList(dict.lookup("masterStickOutFaces"));
		slaveStickOutFacesPtr_ =
			new labelList(dict.lookup("slaveStickOutFaces"));

		retiredPointMapPtr_ = new Map<label>(dict.lookup("retiredPointMap"));
		cutPointEdgePairMapPtr_ =
			new Map<Pair<edge> >(dict.lookup("cutPointEdgePairMap"));
	}
	else
	{
		calcAttachedAddressing();
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::slidingInterface::~slidingInterface()
{
	clearOut();
}


void Foam::slidingInterface::clearAddressing() const
{
	deleteDemandDrivenData(cutFaceMasterPtr_);
	deleteDemandDrivenData(cutFaceSlavePtr_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faceZoneID& Foam::slidingInterface::masterFaceZoneID() const
{
	return masterFaceZoneID_;
}


const Foam::faceZoneID& Foam::slidingInterface::slaveFaceZoneID() const
{
	return slaveFaceZoneID_;
}


const Foam::pointZoneID& Foam::slidingInterface::cutPointZoneID() const
{
	return cutPointZoneID_;
}


const Foam::faceZoneID& Foam::slidingInterface::cutFaceZoneID() const
{
	return cutFaceZoneID_;
}


bool Foam::slidingInterface::changeTopology() const
{
	if (coupleDecouple_)
	{
		// Always changes.  If not attached, project points
		if (debug)
		{
			Pout<< "bool slidingInterface::changeTopology() const "
				<< "for object " << name() << " : "
				<< "Couple-decouple mode." << endl;
		}

		if (!attached_)
		{
			projectPoints();
		}
		else
		{
		}

		return true;
	}

	if
	(
		attached_
	 && !topoChanger().mesh().changing()
	)
	{
		// If the mesh is not moving or morphing and the interface is
		// already attached, the topology will not change
		return false;
	}
	else
	{
		// Check if the motion changes point projection
		return projectPoints();
	}
}


void Foam::slidingInterface::setRefinement(polyTopoChange& ref) const
{
	if (coupleDecouple_)
	{
		if (attached_)
		{
			// Attached, coupling
			decoupleInterface(ref);
		}
		else
		{
			// Detached, coupling
			coupleInterface(ref);
		}

		return;
	}

	if (trigger_)
	{
		if (attached_)
		{
			// Clear old coupling data
			clearCouple(ref);
		}

		coupleInterface(ref);

		trigger_ = false;
	}
}


void Foam::slidingInterface::modifyMotionPoints(pointField& motionPoints) const
{
	if (debug)
	{
		Pout<< "void slidingInterface::modifyMotionPoints("
			<< "pointField& motionPoints) const for object " << name() << " : "
			<< "Adjusting motion points." << endl;
	}

	const polyMesh& mesh = topoChanger().mesh();

	// Get point from the cut zone
	const labelList& cutPoints = mesh.pointZones()[cutPointZoneID_.index()];

	if (cutPoints.size() > 0 && !projectedSlavePointsPtr_)
	{
		return;
	}
	else
	{
		const pointField& projectedSlavePoints = *projectedSlavePointsPtr_;

		const Map<label>& rpm = retiredPointMap();

		const Map<Pair<edge> >& cpepm = cutPointEdgePairMap();

		const Map<label>& slaveZonePointMap =
			mesh.faceZones()[slaveFaceZoneID_.index()]().meshPointMap();

		const primitiveFacePatch& masterPatch =
			mesh.faceZones()[masterFaceZoneID_.index()]();
		const edgeList& masterEdges = masterPatch.edges();
		const pointField& masterLocalPoints = masterPatch.localPoints();

		const primitiveFacePatch& slavePatch =
			mesh.faceZones()[slaveFaceZoneID_.index()]();
		const edgeList& slaveEdges = slavePatch.edges();
		const pointField& slaveLocalPoints = slavePatch.localPoints();
		const vectorField& slavePointNormals = slavePatch.pointNormals();

		forAll (cutPoints, pointI)
		{
			// Try to find the cut point in retired points
			Map<label>::const_iterator rpmIter = rpm.find(cutPoints[pointI]);

			if (rpmIter != rpm.end())
			{
				if (debug)
				{
					Pout << "p";
				}

				// Cut point is a retired point
				motionPoints[cutPoints[pointI]] =
					projectedSlavePoints[slaveZonePointMap.find(rpmIter())()];
			}
			else
			{
				// A cut point is not a projected slave point.  Therefore, it
				// must be an edge-to-edge intersection.  HJ, 24/Jul/2003

				Map<Pair<edge> >::const_iterator cpepmIter =
					cpepm.find(cutPoints[pointI]);

				if (cpepmIter != cpepm.end())
				{
//					 Pout << "Need to re-create hit for point " << cutPoints[pointI] << " lookup: " << cpepmIter() << endl;

					// Note.
					// The edge cutting code is repeated in
					// slidingInterface::coupleInterface.  This is done for
					// efficiency reasons and avoids multiple creation of
					// cutting planes.  Please update both simultaneously.
					// HJ, 28/Jul/2003
					const edge& globalMasterEdge = cpepmIter().first();

					const label curMasterEdgeIndex =
						masterPatch.whichEdge
						(
							edge
							(
								masterPatch.whichPoint
								(
									globalMasterEdge.start()
								),
								masterPatch.whichPoint
								(
									globalMasterEdge.end()
								)
							)
						);

					const edge& cme = masterEdges[curMasterEdgeIndex];
//					 Pout << "curMasterEdgeIndex: " << curMasterEdgeIndex << " cme: " << cme << endl;
					const edge& globalSlaveEdge = cpepmIter().second();

					const label curSlaveEdgeIndex =
						slavePatch.whichEdge
						(
							edge
							(
								slavePatch.whichPoint
								(
									globalSlaveEdge.start()
								),
								slavePatch.whichPoint
								(
									globalSlaveEdge.end()
								)
							)
						);

					const edge& curSlaveEdge = slaveEdges[curSlaveEdgeIndex];
//					 Pout << "curSlaveEdgeIndex: " << curSlaveEdgeIndex << " curSlaveEdge: " << curSlaveEdge << endl;
					const point& a = projectedSlavePoints[curSlaveEdge.start()];
					const point& b = projectedSlavePoints[curSlaveEdge.end()];

					point c =
						0.5*
						(
							slaveLocalPoints[curSlaveEdge.start()]
						  + slavePointNormals[curSlaveEdge.start()]
						  + slaveLocalPoints[curSlaveEdge.end()]
						  + slavePointNormals[curSlaveEdge.end()]
						);

					// Create the plane
					plane cutPlane(a, b, c);

					linePointRef curSlaveLine =
						curSlaveEdge.line(slaveLocalPoints);
					const scalar curSlaveLineMag = curSlaveLine.mag();

					scalar cutOnMaster =
						cutPlane.lineIntersect
						(
							cme.line(masterLocalPoints)
						);

					if
					(
						cutOnMaster > edgeEndCutoffTol_()
					 && cutOnMaster < 1.0 - edgeEndCutoffTol_()
					)
					{
						// Master is cut, check the slave
						point masterCutPoint =
							masterLocalPoints[cme.start()]
						  + cutOnMaster*cme.vec(masterLocalPoints);

						pointHit slaveCut =
							curSlaveLine.nearestDist(masterCutPoint);

						if (slaveCut.hit())
						{
							// Strict checking of slave cut to avoid capturing
							// end points.  HJ, 15/Oct/2004
							scalar cutOnSlave =
								(
									(
										slaveCut.hitPoint()
									  - curSlaveLine.start()
									) & curSlaveLine.vec()
								)/sqr(curSlaveLineMag);

							// Calculate merge tolerance from the
							// target edge length
							scalar mergeTol =
								edgeCoPlanarTol_()*mag(b - a);

							if
							(
								cutOnSlave > edgeEndCutoffTol_()
							 && cutOnSlave < 1.0 - edgeEndCutoffTol_()
							 && slaveCut.distance() < mergeTol
							)
							{
								// Cut both master and slave.
								motionPoints[cutPoints[pointI]] =
									masterCutPoint;
							}
						}
						else
						{
							Pout<< "Missed slave edge!!!  This is an error.  "
								<< "Master edge: "
								<< cme.line(masterLocalPoints)
								<< " slave edge: " << curSlaveLine
								<< " point: " << masterCutPoint
								<< " weight: " <<
								(
									(
										slaveCut.missPoint()
									  - curSlaveLine.start()
									) & curSlaveLine.vec()
								)/sqr(curSlaveLineMag)
								<< endl;
						}
					}
					else
					{
						Pout<< "Missed master edge!!!  This is an error"
							<< endl;
					}
				}
				else
				{
					FatalErrorIn
					(
						"void slidingInterface::modifyMotionPoints"
						"(pointField&) const"
					)   << "Cut point " << cutPoints[pointI]
						<< " not recognised as either the projected "
						<< "or as intersection point.  Error in point "
						<< "projection or data mapping"
						<< abort(FatalError);
				}
			}
		}
		if (debug)
		{
			Pout << endl;
		}
	}
}


void Foam::slidingInterface::updateMesh(const mapPolyMesh& m)
{
	if (debug)
	{
		Pout<< "void slidingInterface::updateMesh(const mapPolyMesh& m)"
			<< " const for object " << name() << " : "
			<< "Updating topology." << endl;
	}

	// Mesh has changed topologically.  Update local topological data
	const polyMesh& mesh = topoChanger().mesh();

	masterFaceZoneID_.update(mesh.faceZones());
	slaveFaceZoneID_.update(mesh.faceZones());
	cutPointZoneID_.update(mesh.pointZones());
	cutFaceZoneID_.update(mesh.faceZones());

	masterPatchID_.update(mesh.boundaryMesh());
	slavePatchID_.update(mesh.boundaryMesh());

	if (!attached())
	{
		calcAttachedAddressing();
	}
	else
	{
		renumberAttachedAddressing(m);
	}
}


const Foam::pointField& Foam::slidingInterface::pointProjection() const
{
	if (!projectedSlavePointsPtr_)
	{
		projectPoints();
	}

	return *projectedSlavePointsPtr_;
}


void Foam::slidingInterface::write(Ostream& os) const
{
	os  << nl << type() << nl
		<< name()<< nl
		<< masterFaceZoneID_.name() << nl
		<< slaveFaceZoneID_.name() << nl
		<< cutPointZoneID_.name() << nl
		<< cutFaceZoneID_.name() << nl
		<< masterPatchID_.name() << nl
		<< slavePatchID_.name() << nl
		<< typeOfMatchNames_[matchType_] << nl
		<< coupleDecouple_ << nl
		<< attached_ << endl;
}


void Foam::slidingInterface::writeDict(Ostream& os) const
{
	os  << nl << name() << nl << token::BEGIN_BLOCK << nl
		<< "    type " << type() << token::END_STATEMENT << nl
		<< "    masterFaceZoneName " << masterFaceZoneID_.name()
		<< token::END_STATEMENT << nl
		<< "    slaveFaceZoneName " << slaveFaceZoneID_.name()
		<< token::END_STATEMENT << nl
		<< "    cutPointZoneName " << cutPointZoneID_.name()
		<< token::END_STATEMENT << nl
		<< "    cutFaceZoneName " << cutFaceZoneID_.name()
		<< token::END_STATEMENT << nl
		<< "    masterPatchName " << masterPatchID_.name()
		<< token::END_STATEMENT << nl
		<< "    slavePatchName " << slavePatchID_.name()
		<< token::END_STATEMENT << nl
		<< "    typeOfMatch " << typeOfMatchNames_[matchType_]
		<< token::END_STATEMENT << nl
		<< "    coupleDecouple " << coupleDecouple_
		<< token::END_STATEMENT << nl
		<< "    projection " << intersection::algorithmNames_[projectionAlgo_]
		<< token::END_STATEMENT << nl
		<< "    attached " << attached_
		<< token::END_STATEMENT << nl
		<< "    active " << active()
		<< token::END_STATEMENT << nl;

	if (attached_)
	{
		masterFaceCellsPtr_->writeEntry("masterFaceCells", os);
		slaveFaceCellsPtr_->writeEntry("slaveFaceCells", os);
		masterStickOutFacesPtr_->writeEntry("masterStickOutFaces", os);
		slaveStickOutFacesPtr_->writeEntry("slaveStickOutFaces", os);

		 os << "    retiredPointMap " << retiredPointMap()
			<< token::END_STATEMENT << nl
			<< "    cutPointEdgePairMap " << cutPointEdgePairMap()
			<< token::END_STATEMENT << nl;
	}

	os  << token::END_BLOCK << endl;
}


// ************************************************************************* //
