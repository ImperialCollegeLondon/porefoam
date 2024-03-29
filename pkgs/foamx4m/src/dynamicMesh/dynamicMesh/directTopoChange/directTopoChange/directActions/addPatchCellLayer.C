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

#include "addPatchCellLayer.H"
#include "polyMesh.H"
#include "directTopoChange.H"
#include "meshTools.H"
#include "mapPolyMesh.H"
#include "syncTools.H"
#include "polyAddPoint.H"
#include "polyAddFace.H"
#include "polyModifyFace.H"
#include "polyAddCell.H"
#include "wallPoint.H"
#include "globalIndex.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::addPatchCellLayer, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Calculate global faces per pp edge.
Foam::labelListList Foam::addPatchCellLayer::calcGlobalEdgeFaces
(
	const polyMesh& mesh,
	const globalIndex& globalFaces,
	const indirectPrimitivePatch& pp,
	const labelList& meshEdges
)
{
	//// Determine coupled edges just so we don't have to have storage
	//// for all non-coupled edges.
	//
	//PackedList<1> isCoupledEdge(mesh.nEdges(), 0);
	//
	//const polyBoundaryMesh& patches = mesh.boundaryMesh();
	//
	//forAll(patches, patchI)
	//{
	//    const polyPatch& pp = patches[patchI];
	//
	//    if (pp.coupled())
	//    {
	//        const labelList& meshEdges = pp.meshEdges();
	//
	//        forAll(meshEdges, i)
	//        {
	//            isCoupledEdge.set(meshEdges[i], 1);
	//        }
	//    }
	//}

	// From mesh edge to global face labels. Sized only for pp edges.
	labelListList globalEdgeFaces(mesh.nEdges());

	const labelListList& edgeFaces = pp.edgeFaces();

	forAll(edgeFaces, edgeI)
	{
		label meshEdgeI = meshEdges[edgeI];

		//if (isCoupledEdge.get(meshEdgeI) == 1)
		{
			const labelList& eFaces = edgeFaces[edgeI];

			// Store face and processor as unique tag.
			labelList& globalEFaces = globalEdgeFaces[meshEdgeI];
			globalEFaces.setSize(eFaces.size());
			forAll(eFaces, i)
			{
				globalEFaces[i] =
					globalFaces.toGlobal(pp.addressing()[eFaces[i]]);
			}
		}
	}

	// Synchronise across coupled edges.
	syncTools::syncEdgeList
	(
		mesh,
		globalEdgeFaces,
		uniqueEqOp(),
		labelList(),    // null value
		false           // no separation
	);

	// Extract pp part
	return IndirectList<labelList>(globalEdgeFaces, meshEdges)();
}


Foam::label Foam::addPatchCellLayer::nbrFace
(
	const labelListList& edgeFaces,
	const label edgeI,
	const label faceI
)
{
	const labelList& eFaces = edgeFaces[edgeI];

	if (eFaces.size() == 2)
	{
		return (eFaces[0] != faceI ? eFaces[0] : eFaces[1]);
	}
	else
	{
		return -1;
	}
}


void Foam::addPatchCellLayer::addVertex
(
	const label pointI,
	face& f,
	label& fp
)
{
	if (fp == 0)
	{
		f[fp++] = pointI;
	}
	else
	{
		if (f[fp-1] != pointI && f[0] != pointI)
		{
			f[fp++] = pointI;
		}
	}

	// Check for duplicates.
	if (debug)
	{
		label n = 0;
		for (label i = 0; i < fp; i++)
		{
			if (f[i] == pointI)
			{
				n++;

				if (n == 2)
				{
					f.setSize(fp);
					FatalErrorIn
					(
					    "addPatchCellLayer::addVertex(const label, face&"
					    ", label&)"
					)   << "Point " << pointI << " already present in face "
					    << f << abort(FatalError);
				}
			}
		}
	}
}


// Is edge to the same neighbour? (and needs extrusion and has not been
// dealt with already)
bool Foam::addPatchCellLayer::sameEdgeNeighbour
(
	const indirectPrimitivePatch& pp,
	const labelListList& globalEdgeFaces,
	const boolList& doneEdge,
	const label thisGlobalFaceI,
	const label nbrGlobalFaceI,
	const label edgeI
) const
{
	const edge& e = pp.edges()[edgeI];

	return
		!doneEdge[edgeI]                            // not yet handled
	 && (
			addedPoints_[e[0]].size() != 0          // is extruded
		 || addedPoints_[e[1]].size() != 0
		)
	 && (
			nbrFace(globalEdgeFaces, edgeI, thisGlobalFaceI)
		 == nbrGlobalFaceI  // is to same neighbour
		);
}


// Collect consecutive string of edges that connects the same two
// (possibly coupled) faces. Returns -1 if no unvisited edge can be found.
// Otherwise returns start and end index in face.
Foam::labelPair Foam::addPatchCellLayer::getEdgeString
(
	const indirectPrimitivePatch& pp,
	const labelListList& globalEdgeFaces,
	const boolList& doneEdge,
	const label patchFaceI,
	const label globalFaceI
) const
{
	const labelList& fEdges = pp.faceEdges()[patchFaceI];

	label startFp = -1;
	label endFp = -1;

	// Get edge that hasn't been done yet but needs extrusion
	forAll(fEdges, fp)
	{
		label edgeI = fEdges[fp];
		const edge& e = pp.edges()[edgeI];

		if
		(
			!doneEdge[edgeI]
		 && (
				addedPoints_[e[0]].size() != 0
			 || addedPoints_[e[1]].size() != 0
			)
		)
		{
			startFp = fp;
			break;
		}
	}

	if (startFp != -1)
	{
		// We found an edge that needs extruding but hasn't been done yet.
		// Now find the face on the other side
		label nbrGlobalFaceI = nbrFace
		(
			globalEdgeFaces,
			fEdges[startFp],
			globalFaceI
		);

		if (nbrGlobalFaceI == -1)
		{
			// Proper boundary edge. Only extrude single edge.
			endFp = startFp;
		}
		else
		{
			// Search back for edge
			// - which hasn't been handled yet
			// - with same neighbour
			// - that needs extrusion
			while(true)
			{
				label prevFp = fEdges.rcIndex(startFp);

				if
				(
					!sameEdgeNeighbour
					(
					    pp,
					    globalEdgeFaces,
					    doneEdge,
					    globalFaceI,
					    nbrGlobalFaceI,
					    fEdges[prevFp]
					)
				)
				{
					break;
				}
				startFp = prevFp;
			}

			// Search forward for end of string
			endFp = startFp;
			while(true)
			{
				label nextFp = fEdges.fcIndex(endFp);

				if
				(
					!sameEdgeNeighbour
					(
					    pp,
					    globalEdgeFaces,
					    doneEdge,
					    globalFaceI,
					    nbrGlobalFaceI,
					    fEdges[nextFp]
					)
				)
				{
					break;
				}
				endFp = nextFp;
			}
		}
	}

	return labelPair(startFp, endFp);
}


// Adds a side face i.e. extrudes a patch edge.
Foam::label Foam::addPatchCellLayer::addSideFace
(
	const indirectPrimitivePatch& pp,
	const labelList& patchID,           // prestored patch per pp face
	const labelListList& addedCells,    // per pp face the new extruded cell
	const face& newFace,
	const label ownFaceI,               // pp face that provides owner
	const label nbrFaceI,
	const label patchEdgeI,             // edge to add to
	const label meshEdgeI,              // corresponding mesh edge
	const label layerI,                 // layer
	const label numEdgeFaces,           // number of layers for edge
	directTopoChange& meshMod
) const
{
	// Edge to 'inflate' from
	label inflateEdgeI = -1;

	// Mesh faces using edge
	const labelList& meshFaces = mesh_.edgeFaces()[meshEdgeI];

	forAll(meshFaces, i)
	{
		if (mesh_.isInternalFace(meshFaces[i]))
		{
			// meshEdge uses internal faces so ok to inflate from it
			inflateEdgeI = meshEdgeI;
			break;
		}
	}


	// Get my mesh face and its zone.
	label meshFaceI = pp.addressing()[ownFaceI];
	label zoneI = mesh_.faceZones().whichZone(meshFaceI);

	label addedFaceI = -1;

	// Is patch edge external edge of indirectPrimitivePatch?
	if (nbrFaceI == -1)
	{
		// External edge so external face. Patch id is obtained from
		// any other patch connected to edge.

		const polyBoundaryMesh& patches = mesh_.boundaryMesh();

		// Loop over all faces connected to edge to inflate and
		// see if any boundary face (but not meshFaceI)
		label otherPatchID = patchID[ownFaceI];

		forAll(meshFaces, k)
		{
			label faceI = meshFaces[k];

			if
			(
				faceI != meshFaceI
			 && !mesh_.isInternalFace(faceI)
			)
			{
				otherPatchID = patches.whichPatch(faceI);
				break;
			}
		}

		// Determine if different number of layer on owner and neighbour side
		// (relevant only for coupled faces). See section for internal edge
		// below.

		label layerOwn;

		if (addedCells[ownFaceI].size() < numEdgeFaces)
		{
			label offset = numEdgeFaces - addedCells[ownFaceI].size();
			if (layerI <= offset)
			{
				layerOwn = 0;
			}
			else
			{
				layerOwn = layerI - offset;
			}
		}
		else
		{
			layerOwn = layerI;
		}


		//Pout<< "Added boundary face:" << newFace
		//    << " own:" << addedCells[ownFaceI][layerOwn]
		//    << " patch:" << otherPatchID
		//    << endl;

		addedFaceI = meshMod.setAction
		(
			polyAddFace
			(
				newFace,                    // face
				addedCells[ownFaceI][layerOwn],   // owner
				-1,                         // neighbour
				-1,                         // master point
				inflateEdgeI,               // master edge
				-1,                         // master face
				false,                      // flux flip
				otherPatchID,               // patch for face
				zoneI,                      // zone for face
				false                       // face zone flip
			)
		);
	}
	else
	{
		// When adding side faces we need to modify neighbour and owners
		// in region where layer mesh is stopped. Determine which side
		// has max number of faces and make sure layers match closest to
		// original pp if there are different number of layers.

		label layerNbr;
		label layerOwn;

		if (addedCells[ownFaceI].size() > addedCells[nbrFaceI].size())
		{
			label offset =
				addedCells[ownFaceI].size() - addedCells[nbrFaceI].size();

			layerOwn = layerI;

			if (layerI <= offset)
			{
				layerNbr = 0;
			}
			else
			{
				layerNbr = layerI - offset;
			}
		}
		else if (addedCells[nbrFaceI].size() > addedCells[ownFaceI].size())
		{
			label offset =
				addedCells[nbrFaceI].size() - addedCells[ownFaceI].size();

			layerNbr = layerI;

			if (layerI <= offset)
			{
				layerOwn = 0;
			}
			else
			{
				layerOwn = layerI - offset;
			}
		}
		else
		{
			// Same number of layers on both sides.
			layerNbr = layerI;
			layerOwn = layerI;
		}

		addedFaceI = meshMod.setAction
		(
			polyAddFace
			(
				newFace,                    // face
				addedCells[ownFaceI][layerOwn],   // owner
				addedCells[nbrFaceI][layerNbr],   // neighbour
				-1,                         // master point
				inflateEdgeI,               // master edge
				-1,                         // master face
				false,                      // flux flip
				-1,                         // patch for face
				zoneI,                      // zone for face
				false                       // face zone flip
			)
		);

	   //Pout<< "Added internal face:" << newFace
		//    << " own:" << addedCells[ownFaceI][layerOwn]
		//    << " nei:" << addedCells[nbrFaceI][layerNbr]
		//    << endl;
	}

	return addedFaceI;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh
Foam::addPatchCellLayer::addPatchCellLayer(const polyMesh& mesh)
:
	mesh_(mesh),
	addedPoints_(0),
	layerFaces_(0)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::labelListList Foam::addPatchCellLayer::addedCells
(
	const polyMesh& mesh,
	const labelListList& layerFaces
)
{
	labelListList layerCells(layerFaces.size());

	forAll(layerFaces, patchFaceI)
	{
		const labelList& faceLabels = layerFaces[patchFaceI];

		if (faceLabels.size() > 0)
		{
			labelList& added = layerCells[patchFaceI];
			added.setSize(faceLabels.size()-1);

			for (label i = 0; i < faceLabels.size()-1; i++)
			{
				added[i] = mesh.faceNeighbour()[faceLabels[i]];
			}
		}
	}
	return layerCells;
}


Foam::labelListList Foam::addPatchCellLayer::addedCells() const
{
	return addedCells(mesh_, layerFaces_);
}


void Foam::addPatchCellLayer::setRefinement
(
	const scalarField& expansionRatio,
	const indirectPrimitivePatch& pp,
	const labelList& nFaceLayers,
	const labelList& nPointLayers,
	const vectorField& firstLayerDisp,
	directTopoChange& meshMod
)
{
	if (debug)
	{
		Pout<< "addPatchCellLayer::setRefinement : Adding up to "
			<< max(nPointLayers)
			<< " layers of cells to indirectPrimitivePatch with "
			<< pp.nPoints() << " points" << endl;
	}

	if
	(
		pp.nPoints() != firstLayerDisp.size()
	 || pp.nPoints() != nPointLayers.size()
	 || pp.size() != nFaceLayers.size()
	)
	{
		FatalErrorIn
		(
			"addPatchCellLayer::setRefinement"
			"(const scalar, const indirectPrimitivePatch&"
			", const labelList&, const vectorField&, directTopoChange&)"
		)   << "Size of new points is not same as number of points used by"
			<< " the face subset" << endl
			<< "  patch.nPoints:" << pp.nPoints()
			<< "  displacement:" << firstLayerDisp.size()
			<< "  nPointLayers:" << nPointLayers.size() << nl
			<< " patch.nFaces:" << pp.size()
			<< "  nFaceLayers:" << nFaceLayers.size()
			<< abort(FatalError);
	}

	forAll(nPointLayers, i)
	{
		if (nPointLayers[i] < 0)
		{
			FatalErrorIn
			(
				"addPatchCellLayer::setRefinement"
				"(const scalar, const indirectPrimitivePatch&"
				", const labelList&, const vectorField&, directTopoChange&)"
			)   << "Illegal number of layers " << nPointLayers[i]
				<< " at patch point " << i << abort(FatalError);
		}
	}
	forAll(nFaceLayers, i)
	{
		if (nFaceLayers[i] < 0)
		{
			FatalErrorIn
			(
				"addPatchCellLayer::setRefinement"
				"(const scalar, const indirectPrimitivePatch&"
				", const labelList&, const vectorField&, directTopoChange&)"
			)   << "Illegal number of layers " << nFaceLayers[i]
				<< " at patch face " << i << abort(FatalError);
		}
	}

	const labelList& meshPoints = pp.meshPoints();

	// Precalculate mesh edges for pp.edges.
	labelList meshEdges(calcMeshEdges(mesh_, pp));

	if (debug)
	{
		// Check synchronisation
		// ~~~~~~~~~~~~~~~~~~~~~

		{
			labelList n(mesh_.nPoints(), 0);
			IndirectList<label>(n, meshPoints) = nPointLayers;
			syncTools::syncPointList
			(
				mesh_,
				n,
				maxEqOp<label>(),
				label(0),
				false
			);

			// Non-synced
			forAll(meshPoints, i)
			{
				label meshPointI = meshPoints[i];

				if (n[meshPointI] != nPointLayers[i])
				{
					FatalErrorIn
					(
					    "addPatchCellLayer::setRefinement"
					    "(const scalar, const indirectPrimitivePatch&"
					    ", const labelList&, const vectorField&"
					    ", directTopoChange&)"
					)   << "At mesh point:" << meshPointI
					    << " coordinate:" << mesh_.points()[meshPointI]
					    << " specified nLayers:" << nPointLayers[i] << endl
					    << "On coupled point a different nLayers:"
					    << n[meshPointI] << " was specified."
					    << abort(FatalError);
				}
			}


			// Check that nPointLayers equals the max layers of connected faces
			// (or 0). Anything else makes no sense.
			labelList nFromFace(mesh_.nPoints(), 0);
			forAll(nFaceLayers, i)
			{
				const face& f = pp[i];

				forAll(f, fp)
				{
					label pointI = f[fp];

					nFromFace[pointI] = max(nFromFace[pointI], nFaceLayers[i]);
				}
			}
			syncTools::syncPointList
			(
				mesh_,
				nFromFace,
				maxEqOp<label>(),
				label(0),
				false
			);

			forAll(nPointLayers, i)
			{
				label meshPointI = meshPoints[i];

				if
				(
					nPointLayers[i] > 0
				 && nPointLayers[i] != nFromFace[meshPointI]
				)
				{
					FatalErrorIn
					(
					    "addPatchCellLayer::setRefinement"
					    "(const scalar, const indirectPrimitivePatch&"
					    ", const labelList&, const vectorField&"
					    ", directTopoChange&)"
					)   << "At mesh point:" << meshPointI
					    << " coordinate:" << mesh_.points()[meshPointI]
					    << " specified nLayers:" << nPointLayers[i] << endl
					    << "but the max nLayers of surrounding faces is:"
					    << nFromFace[meshPointI]
					    << abort(FatalError);
				}
			}
		}

		{
			pointField d(mesh_.nPoints(), wallPoint::greatPoint);
			IndirectList<point>(d, meshPoints) = firstLayerDisp;
			syncTools::syncPointList
			(
				mesh_,
				d,
				minEqOp<vector>(),
				wallPoint::greatPoint,
				false
			);

			forAll(meshPoints, i)
			{
				label meshPointI = meshPoints[i];

				if (mag(d[meshPointI] - firstLayerDisp[i]) > SMALL)
				{
					FatalErrorIn
					(
					    "addPatchCellLayer::setRefinement"
					    "(const scalar, const indirectPrimitivePatch&"
					    ", const labelList&, const vectorField&"
					    ", directTopoChange&)"
					)   << "At mesh point:" << meshPointI
					    << " coordinate:" << mesh_.points()[meshPointI]
					    << " specified displacement:" << firstLayerDisp[i]
					    << endl
					    << "On coupled point a different displacement:"
					    << d[meshPointI] << " was specified."
					    << abort(FatalError);
				}
			}
		}

		// Check that edges of pp (so ones that become boundary faces)
		// connect to only one boundary face. Guarantees uniqueness of
		// patch that they go into so if this is a coupled patch both
		// sides decide the same.
		// ~~~~~~~~~~~~~~~~~~~~~~

		for (label edgeI = pp.nInternalEdges(); edgeI < pp.nEdges(); edgeI++)
		{
			const edge& e = pp.edges()[edgeI];

			if (nPointLayers[e[0]] > 0 || nPointLayers[e[1]] > 0)
			{
				// Edge is to become a face

				const labelList& eFaces = pp.edgeFaces()[edgeI];

				// First check: pp should be single connected.
				if (eFaces.size() != 1)
				{
					FatalErrorIn
					(
					    "addPatchCellLayer::setRefinement"
					    "(const scalar, const indirectPrimitivePatch&"
					    ", const labelList&, const vectorField&"
					    ", directTopoChange&)"
					)   << "boundary-edge-to-be-extruded:"
					    << pp.points()[meshPoints[e[0]]]
					    << pp.points()[meshPoints[e[0]]]
					    << " has more than two faces using it:" << eFaces
					    << abort(FatalError);
				}

				label myFaceI = pp.addressing()[eFaces[0]];

				label meshEdgeI = meshEdges[edgeI];

				// Mesh faces using edge
				const labelList& meshFaces = mesh_.edgeFaces()[meshEdgeI];

				// Check that there is only one patchface using edge.
				const polyBoundaryMesh& patches = mesh_.boundaryMesh();

				label bFaceI = -1;

				forAll(meshFaces, i)
				{
					label faceI = meshFaces[i];

					if (faceI != myFaceI)
					{
					    if (!mesh_.isInternalFace(faceI))
					    {
					        if (bFaceI == -1)
					        {
					            bFaceI = faceI;
					        }
					        else
					        {
					            FatalErrorIn
					            (
					                "addPatchCellLayer::setRefinement"
					                "(const scalar"
					                ", const indirectPrimitivePatch&"
					                ", const labelList&, const vectorField&"
					                ", directTopoChange&)"
					            )   << "boundary-edge-to-be-extruded:"
					                << pp.points()[meshPoints[e[0]]]
					                << pp.points()[meshPoints[e[0]]]
					                << " has more than two boundary faces"
					                << " using it:"
					                << bFaceI << " fc:"
					                << mesh_.faceCentres()[bFaceI]
					                << " patch:" << patches.whichPatch(bFaceI)
					                << " and " << faceI << " fc:"
					                << mesh_.faceCentres()[faceI]
					                << " patch:" << patches.whichPatch(faceI)
					                << abort(FatalError);
					        }
					    }
					}
				}
			}
		}
	}


	// From master point (in patch point label) to added points (in mesh point
	// label)
	addedPoints_.setSize(pp.nPoints());

	// Mark points that do not get extruded by setting size of addedPoints_ to 0
	label nTruncated = 0;

	forAll(nPointLayers, patchPointI)
	{
		if (nPointLayers[patchPointI] > 0)
		{
			addedPoints_[patchPointI].setSize(nPointLayers[patchPointI]);
		}
		else
		{
			nTruncated++;
		}
	}

	if (debug)
	{
		Pout<< "Not adding points at " << nTruncated << " out of "
			<< pp.nPoints() << " points" << endl;
	}


	//
	// Create new points
	//

	forAll(firstLayerDisp, patchPointI)
	{
		if (addedPoints_[patchPointI].size() > 0)
		{
			label meshPointI = meshPoints[patchPointI];

			label zoneI = mesh_.pointZones().whichZone(meshPointI);

			point pt = mesh_.points()[meshPointI];

			vector disp = firstLayerDisp[patchPointI];

			for (label i = 0; i < addedPoints_[patchPointI].size(); i++)
			{
				pt += disp;

				label addedVertI = meshMod.setAction
				(
					polyAddPoint
					(
					    pt,         // point
					    meshPointI, // master point
					    zoneI,      // zone for point
					    true        // supports a cell
					)
				);

				addedPoints_[patchPointI][i] = addedVertI;

				disp *= expansionRatio[patchPointI];
			}
		}
	}


	//
	// Add cells to all boundaryFaces
	//

	labelListList addedCells(pp.size());

	forAll(pp, patchFaceI)
	{
		if (nFaceLayers[patchFaceI] > 0)
		{
			addedCells[patchFaceI].setSize(nFaceLayers[patchFaceI]);

			label meshFaceI = pp.addressing()[patchFaceI];

			label ownZoneI = mesh_.cellZones().whichZone
			(
				mesh_.faceOwner()[meshFaceI]
			);

			for (label i = 0; i < nFaceLayers[patchFaceI]; i++)
			{
				// Note: add from cell (owner of patch face) or from face?
				// for now add from cell so we can map easily.
				addedCells[patchFaceI][i] = meshMod.setAction
				(
					polyAddCell
					(
					    -1,             // master point
					    -1,             // master edge
					    -1,             // master face
					    mesh_.faceOwner()[meshFaceI],   // master cell id
					    ownZoneI        // zone for cell
					)
				);

				//Pout<< "For patchFace:" << patchFaceI
				//    << " meshFace:" << pp.addressing()[patchFaceI]
				//    << " layer:" << i << " added cell:"
				//    << addedCells[patchFaceI][i]
				//    << endl;
			}
		}
	}


	const polyBoundaryMesh& patches = mesh_.boundaryMesh();

	// Precalculated patchID for each patch face
	labelList patchID(pp.size());

	forAll(pp, patchFaceI)
	{
		label meshFaceI = pp.addressing()[patchFaceI];

		patchID[patchFaceI] = patches.whichPatch(meshFaceI);
	}



	// Create faces on top of the original patch faces.
	// These faces are created from original patch faces outwards so in order
	// of increasing cell number. So orientation should be same as original
	// patch face for them to have owner<neighbour.

	layerFaces_.setSize(pp.size());

	forAll(pp.localFaces(), patchFaceI)
	{
		label meshFaceI = pp.addressing()[patchFaceI];

		if (addedCells[patchFaceI].size() > 0)
		{
			layerFaces_[patchFaceI].setSize(addedCells[patchFaceI].size() + 1);
			layerFaces_[patchFaceI][0] = meshFaceI;

			label zoneI = mesh_.faceZones().whichZone(meshFaceI);

			// Get duplicated vertices on the patch face.
			const face& f = pp.localFaces()[patchFaceI];

			face newFace(f.size());

			for (label i = 0; i < addedCells[patchFaceI].size(); i++)
			{
				forAll(f, fp)
				{
					if (addedPoints_[f[fp]].size() == 0)
					{
					    // Keep original point
					    newFace[fp] = meshPoints[f[fp]];
					}
					else
					{
					    // Get new outside point
					    label offset =
					        addedPoints_[f[fp]].size()
					      - addedCells[patchFaceI].size();
					    newFace[fp] = addedPoints_[f[fp]][i+offset];
					}
				}


				// Get new neighbour
				label nei;
				label patchI;

				if (i == addedCells[patchFaceI].size()-1)
				{
					// Top layer so is patch face.
					nei = -1;
					patchI = patchID[patchFaceI];
				}
				else
				{
					// Internal face between layer i and i+1
					nei = addedCells[patchFaceI][i+1];
					patchI = -1;
				}


				layerFaces_[patchFaceI][i+1] = meshMod.setAction
				(
					polyAddFace
					(
					    newFace,                    // face
					    addedCells[patchFaceI][i],  // owner
					    nei,                        // neighbour
					    -1,                         // master point
					    -1,                         // master edge
					    meshFaceI,                  // master face for addition
					    false,                      // flux flip
					    patchI,                     // patch for face
					    zoneI,                      // zone for face
					    false                       // face zone flip
					)
				);
				//Pout<< "Added inbetween face " << newFace
				//    << " own:" << addedCells[patchFaceI][i]
				//    << " nei:" << nei
				//    << " patch:" << patchI
				//    << endl;
			}
		}
	}

	//
	// Modify old patch faces to be on the inside
	//
	forAll(pp, patchFaceI)
	{
		if (addedCells[patchFaceI].size() > 0)
		{
			label meshFaceI = pp.addressing()[patchFaceI];

			label zoneI = mesh_.faceZones().whichZone(meshFaceI);

			meshMod.setAction
			(
				polyModifyFace
				(
					pp[patchFaceI],                 // modified face
					meshFaceI,                      // label of face
					mesh_.faceOwner()[meshFaceI],   // owner
					addedCells[patchFaceI][0],      // neighbour
					false,                          // face flip
					-1,                             // patch for face
					false,                          // remove from zone
					zoneI,                          // zone for face
					false                           // face flip in zone
				)
			);
			//Pout<< "Modified old patch face " << meshFaceI
			//    << " own:" << mesh_.faceOwner()[meshFaceI]
			//    << " nei:" << addedCells[patchFaceI][0]
			//    << endl;
		}
	}


	//
	// Create 'side' faces, one per edge that is being extended.
	//

	const labelListList& faceEdges = pp.faceEdges();
	const faceList& localFaces = pp.localFaces();
	const edgeList& edges = pp.edges();

	// Get number of layers per edge. This is 0 if edge is not extruded;
	// max of connected faces otherwise.
	labelList edgeLayers(pp.nEdges());

	{
		// Use list over mesh.nEdges() since syncTools does not yet support
		// partial list synchronisation.
		labelList meshEdgeLayers(mesh_.nEdges(), -1);

		forAll(meshEdges, edgeI)
		{
			const edge& e = edges[edgeI];

			label meshEdgeI = meshEdges[edgeI];

			if ((nPointLayers[e[0]] == 0) && (nPointLayers[e[1]] == 0))
			{
				meshEdgeLayers[meshEdgeI] = 0;
			}
			else
			{
				const labelList& eFaces = pp.edgeFaces()[edgeI];

				forAll(eFaces, i)
				{
					meshEdgeLayers[meshEdgeI] = max
					(
					    nFaceLayers[eFaces[i]],
					    meshEdgeLayers[meshEdgeI]
					);
				}
			}
		}

		syncTools::syncEdgeList
		(
			mesh_,
			meshEdgeLayers,
			maxEqOp<label>(),
			label(0),           // initial value
			false               // no separation
		);

		forAll(meshEdges, edgeI)
		{
			edgeLayers[edgeI] = meshEdgeLayers[meshEdges[edgeI]];
		}
	}


	// Global indices engine
	const globalIndex globalFaces(mesh_.nFaces());

	// Get for all pp edgeFaces a unique faceID
	labelListList globalEdgeFaces
	(
		 calcGlobalEdgeFaces
		 (
			mesh_,
			globalFaces,
			pp,
			meshEdges
		)
	);


	// Mark off which edges have been extruded
	boolList doneEdge(pp.nEdges(), false);


	// Create faces. Per face walk connected edges and find string of edges
	// between the same two faces and extrude string into a single face.
	forAll(pp, patchFaceI)
	{
		const labelList& fEdges = faceEdges[patchFaceI];

		forAll(fEdges, fp)
		{
			// Get string of edges that needs to be extruded as a single face.
			// Returned as indices in fEdges.
			labelPair indexPair
			(
				getEdgeString
				(
					pp,
					globalEdgeFaces,
					doneEdge,
					patchFaceI,
					globalFaces.toGlobal(pp.addressing()[patchFaceI])
				)
			);

			//Pout<< "Found unextruded edges in edges:" << fEdges
			//    << " start:" << indexPair[0]
			//    << " end:" << indexPair[1]
			//    << endl;

			const label startFp = indexPair[0];
			const label endFp = indexPair[1];

			if (startFp != -1)
			{
				// Extrude edges from indexPair[0] up to indexPair[1]
				// (note indexPair = indices of edges. There is one more vertex
				//  than edges)
				const face& f = localFaces[patchFaceI];

				labelList stringedVerts;
				if (endFp >= startFp)
				{
					stringedVerts.setSize(endFp-startFp+2);
				}
				else
				{
					stringedVerts.setSize(endFp+f.size()-startFp+2);
				}

				label fp = startFp;

				for (label i = 0; i < stringedVerts.size()-1; i++)
				{
					stringedVerts[i] = f[fp];
					doneEdge[fEdges[fp]] = true;
					fp = f.fcIndex(fp);
				}
				stringedVerts[stringedVerts.size()-1] = f[fp];


				// Now stringedVerts contains the vertices in order of face f.
				// This is consistent with the order if f becomes the owner cell
				// and nbrFaceI the neighbour cell. Note that the cells get
				// added in order of pp so we can just use face ordering and
				// because we loop in incrementing order as well we will
				// always have nbrFaceI > patchFaceI.

				label startEdgeI = fEdges[startFp];

				label meshEdgeI = meshEdges[startEdgeI];

				label numEdgeSideFaces = edgeLayers[startEdgeI];

				for (label i = 0; i < numEdgeSideFaces; i++)
				{
					label vEnd = stringedVerts[stringedVerts.size()-1];
					label vStart = stringedVerts[0];

					// calculate number of points making up a face
					label newFp = 2*stringedVerts.size();

					if (i == 0)
					{
					    // layer 0 gets all the truncation of neighbouring
					    // faces with more layers.
					    if (addedPoints_[vEnd].size() != 0)
					    {
					        newFp +=
					            addedPoints_[vEnd].size() - numEdgeSideFaces;
					    }
					    if (addedPoints_[vStart].size() != 0)
					    {
					        newFp +=
					            addedPoints_[vStart].size()  - numEdgeSideFaces;
					    }
					}

					face newFace(newFp);

					newFp = 0;

					// For layer 0 get pp points, for all other layers get
					// points of layer-1.
					if (i == 0)
					{
					    forAll(stringedVerts, stringedI)
					    {
					        label v = stringedVerts[stringedI];
					        addVertex(meshPoints[v], newFace, newFp);
					    }
					}
					else
					{
					    forAll(stringedVerts, stringedI)
					    {
					        label v = stringedVerts[stringedI];
					        if (addedPoints_[v].size() > 0)
					        {
					            label offset =
					                addedPoints_[v].size() - numEdgeSideFaces;
					            addVertex
					            (
					                addedPoints_[v][i+offset-1],
					                newFace,
					                newFp
					            );
					        }
					        else
					        {
					            addVertex(meshPoints[v], newFace, newFp);
					        }
					    }
					}

					// add points between stringed vertices (end)
					if (numEdgeSideFaces < addedPoints_[vEnd].size())
					{
					    if (i == 0 && addedPoints_[vEnd].size() != 0)
					    {
					        label offset =
					            addedPoints_[vEnd].size() - numEdgeSideFaces;
					        for (label ioff = 0; ioff < offset; ioff++)
					        {
					            addVertex
					            (
					                addedPoints_[vEnd][ioff],
					                newFace,
					                newFp
					            );
					        }
					    }
					}

					forAllReverse(stringedVerts, stringedI)
					{
					    label v = stringedVerts[stringedI];
					    if (addedPoints_[v].size() > 0)
					    {
					        label offset =
					            addedPoints_[v].size() - numEdgeSideFaces;
					        addVertex
					        (
					            addedPoints_[v][i+offset],
					            newFace,
					            newFp
					        );
					    }
					    else
					    {
					        addVertex(meshPoints[v], newFace, newFp);
					    }
					}


					// add points between stringed vertices (start)
					if (numEdgeSideFaces < addedPoints_[vStart].size())
					{
					    if (i == 0 && addedPoints_[vStart].size() != 0)
					    {
					        label offset =
					            addedPoints_[vStart].size() - numEdgeSideFaces;
					        for (label ioff = offset; ioff > 0; ioff--)
					        {
					            addVertex
					            (
					                addedPoints_[vStart][ioff-1],
					                newFace,
					                newFp
					            );
					        }
					    }
					}

					if (newFp >= 3)
					{
					    // Add face inbetween faces patchFaceI and nbrFaceI
					    // (possibly -1 for external edges)

					    newFace.setSize(newFp);

					    label nbrFaceI = nbrFace
					    (
					        pp.edgeFaces(),
					        startEdgeI,
					        patchFaceI
					    );

					    addSideFace
					    (
					        pp,
					        patchID,
					        addedCells,
					        newFace,
					        patchFaceI,
					        nbrFaceI,
					        startEdgeI,     // edge to inflate from
					        meshEdgeI,      // corresponding mesh edge
					        i,
					        numEdgeSideFaces,
					        meshMod
					    );
					}
				}
			}
		}
	}
}


void Foam::addPatchCellLayer::updateMesh
(
	const mapPolyMesh& morphMap,
	const labelList& faceMap,   // new to old patch faces
	const labelList& pointMap   // new to old patch points
)
{
	{
		labelListList newAddedPoints(pointMap.size());

		forAll(newAddedPoints, newPointI)
		{
			label oldPointI = pointMap[newPointI];

			const labelList& added = addedPoints_[oldPointI];

			labelList& newAdded = newAddedPoints[newPointI];
			newAdded.setSize(added.size());
			label newI = 0;

			forAll(added, i)
			{
				label newPointI = morphMap.reversePointMap()[added[i]];

				if (newPointI >= 0)
				{
					newAdded[newI++] = newPointI;
				}
			}
			newAdded.setSize(newI);
		}
		addedPoints_.transfer(newAddedPoints);
	}

	{
		labelListList newLayerFaces(faceMap.size());

		forAll(newLayerFaces, newFaceI)
		{
			label oldFaceI = faceMap[newFaceI];

			const labelList& added = layerFaces_[oldFaceI];

			labelList& newAdded = newLayerFaces[newFaceI];
			newAdded.setSize(added.size());
			label newI = 0;

			forAll(added, i)
			{
				label newFaceI = morphMap.reverseFaceMap()[added[i]];

				if (newFaceI >= 0)
				{
					newAdded[newI++] = newFaceI;
				}
			}
			newAdded.setSize(newI);
		}
		layerFaces_.transfer(newLayerFaces);
	}
}


Foam::labelList Foam::addPatchCellLayer::calcMeshEdges
(
	const primitiveMesh& mesh,
	const indirectPrimitivePatch& pp
)
{
	labelList meshEdges(pp.nEdges());

	forAll(meshEdges, patchEdgeI)
	{
		const edge& e = pp.edges()[patchEdgeI];

		label v0 = pp.meshPoints()[e[0]];
		label v1 = pp.meshPoints()[e[1]];
		meshEdges[patchEdgeI] = meshTools::findEdge
		(
			mesh.edges(),
			mesh.pointEdges()[v0],
			v0,
			v1
		);
	}
	return meshEdges;
}


// ************************************************************************* //
