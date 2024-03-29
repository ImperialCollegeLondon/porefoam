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

#include "intersectedSurface.H"
#include "surfaceIntersection.H"
#include "faceList.H"
#include "faceTriangulation.H"
#include "treeBoundBox.H"
#include "OFstream.H"
#include "error.H"
#include "meshTools.H"
#include "edgeSurface.H"
#include "DynamicList.H"
#include "transform.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::intersectedSurface, 0);

const Foam::label Foam::intersectedSurface::UNVISITED = 0;
const Foam::label Foam::intersectedSurface::STARTTOEND = 1;
const Foam::label Foam::intersectedSurface::ENDTOSTART = 2;
const Foam::label Foam::intersectedSurface::BOTH = STARTTOEND | ENDTOSTART;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Write whole pointField and edges to stream
void Foam::intersectedSurface::writeOBJ
(
	const pointField& points,
	const edgeList& edges,
	Ostream& os
)
{
	forAll(points, pointI)
	{
		const point& pt = points[pointI];

		os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
	}
	forAll(edges, edgeI)
	{
		const edge& e = edges[edgeI];

		os << "l " << e.start()+1 << ' ' << e.end()+1 << nl;
	}
}


// Write whole pointField and selected edges to stream
void Foam::intersectedSurface::writeOBJ
(
	const pointField& points,
	const edgeList& edges,
	const labelList& faceEdges,
	Ostream& os
)
{
	forAll(points, pointI)
	{
		const point& pt = points[pointI];

		os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
	}
	forAll(faceEdges, i)
	{
		const edge& e = edges[faceEdges[i]];

		os << "l " << e.start()+1 << ' ' << e.end()+1 << nl;
	}
}


// write local points and edges to stream
void Foam::intersectedSurface::writeLocalOBJ
(
	const pointField& points,
	const edgeList& edges,
	const labelList& faceEdges,
	const fileName& fName
)
{
	OFstream os(fName);

	labelList pointMap(points.size(), -1);

	label maxVertI = 0;

	forAll(faceEdges, i)
	{
		const edge& e = edges[faceEdges[i]];

		forAll(e, i)
		{
			label pointI = e[i];

			if (pointMap[pointI] == -1)
			{
				const point& pt = points[pointI];

				os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;

				pointMap[pointI] = maxVertI++;
			}
		}
	}

	forAll(faceEdges, i)
	{
		const edge& e = edges[faceEdges[i]];

		os << "l " << pointMap[e.start()]+1 << ' ' << pointMap[e.end()]+1
			<< nl;
	}
}


// Write whole pointField and face to stream
void Foam::intersectedSurface::writeOBJ
(
	const pointField& points,
	const face& f,
	Ostream& os
)
{
	forAll(points, pointI)
	{
		const point& pt = points[pointI];

		os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
	}

	os << 'f';

	forAll(f, fp)
	{
		os << ' ' << f[fp]+1;
	}
	os << nl;
}


// Print current visited state.
void Foam::intersectedSurface::printVisit
(
	const edgeList& edges,
	const labelList& edgeLabels,
	const Map<label>& visited
)
{
	Pout<< "Visited:" << nl;
	forAll(edgeLabels, i)
	{
		label edgeI = edgeLabels[i];

		const edge& e = edges[edgeI];

		label stat = visited[edgeI];

		if (stat == UNVISITED)
		{
			Pout<< "    edge:" << edgeI << "  verts:" << e
				<< "  unvisited" << nl;
		}
		else if (stat == STARTTOEND)
		{
			Pout<< "    edge:" << edgeI << "  from " << e[0]
				<< " to " << e[1] << nl;
		}
		else if (stat == ENDTOSTART)
		{
			Pout<< "    edge:" << edgeI << "  from " << e[1]
				<< " to " << e[0] << nl;
		}
		else
		{
			Pout<< "    edge:" << edgeI << "  both " << e
				<< nl;
		}
	}
	Pout<< endl;
}



// Check if the two vertices that f0 and f1 share are in the same order on
// both faces.
bool Foam::intersectedSurface::sameEdgeOrder
(
	const labelledTri& fA,
	const labelledTri& fB
)
{
	forAll(fA, fpA)
	{
		label fpB = findIndex(fB, fA[fpA]);

		if (fpB != -1)
		{
			// Get prev/next vertex on fA
			label vA1 = fA[(fpA + 1) % 3];
			label vAMin1 = fA[fpA ? fpA-1 : 2];

			// Get prev/next vertex on fB
			label vB1 = fB[(fpB + 1) % 3];
			label vBMin1 = fB[fpB ? fpB-1 : 2];

			if (vA1 == vB1 || vAMin1 == vBMin1)
			{
				return true;
			}
			else if (vA1 == vBMin1 || vAMin1 == vB1)
			{
				// shared vertices in opposite order.
				return false;
			}
			else
			{
				FatalErrorIn("intersectedSurface::sameEdgeOrder")
					<< "Triangle:" << fA << " and triangle:" << fB
					<< " share a point but not an edge"
					<< abort(FatalError);
			}
		}
	}

	FatalErrorIn("intersectedSurface::sameEdgeOrder")
		<< "Triangle:" << fA << " and triangle:" << fB
		<< " do not share an edge"
		<< abort(FatalError);

	return false;
}


void Foam::intersectedSurface::incCount
(
	Map<label>& visited,
	const label key,
	const label offset
)
{
	Map<label>::iterator iter = visited.find(key);

	if (iter == visited.end())
	{
		visited.insert(key, offset);
	}
	else
	{
		iter() += offset;
	}
}


// Calculate point to edge addressing for the face given by the edge
// subset faceEdges. Constructs facePointEdges which for every point
// gives a list of edge labels connected to it.
Foam::Map<Foam::DynamicList<Foam::label> >
Foam::intersectedSurface::calcPointEdgeAddressing
(
	const edgeSurface& eSurf,
	const label faceI
)
{
	const pointField& points = eSurf.points();
	const edgeList& edges = eSurf.edges();

	const labelList& fEdges = eSurf.faceEdges()[faceI];

	Map<dynamicLabelList > facePointEdges(4*fEdges.size());

	forAll(fEdges, i)
	{
		label edgeI = fEdges[i];

		const edge& e = edges[edgeI];

		// Add e.start to point-edges
		Map<dynamicLabelList >::iterator iter =
			facePointEdges.find(e.start());

		if (iter == facePointEdges.end())
		{
			dynamicLabelList oneEdge;
			oneEdge.append(edgeI);
			facePointEdges.insert(e.start(), oneEdge);
		}
		else
		{
			iter().append(edgeI);
		}

		// Add e.end to point-edges
		Map<dynamicLabelList >::iterator iter2 =
			facePointEdges.find(e.end());

		if (iter2 == facePointEdges.end())
		{
			dynamicLabelList oneEdge;
			oneEdge.append(edgeI);
			facePointEdges.insert(e.end(), oneEdge);
		}
		else
		{
			iter2().append(edgeI);
		}
	}

	// Shrink it
	for
	(
		Map<dynamicLabelList >::iterator iter = facePointEdges.begin();
		iter != facePointEdges.end();
		++iter
	)
	{
		iter().shrink();

		// Check on dangling points.
		if (iter().empty())
		{
			FatalErrorIn
			(
				"intersectedSurface::calcPointEdgeAddressing"
				"(const edgeSurface&, const label)"
			)   << "Point:" << iter.key() << " used by too few edges:"
				<< iter() << abort(FatalError);
		}
	}

	if (debug & 2)
	{
		// Print facePointEdges
		Pout<< "calcPointEdgeAddressing: face consisting of edges:" << endl;
		forAll(fEdges, i)
		{
			label edgeI = fEdges[i];
			const edge& e = edges[edgeI];
			Pout<< "    " << edgeI << ' ' << e << points[e.start()]
				<< points[e.end()] << endl;
		}

		Pout<< "    Constructed point-edge adressing:" << endl;
		for
		(
			Map<dynamicLabelList >::iterator iter = facePointEdges.begin();
			iter != facePointEdges.end();
			++iter
		)
		{
			Pout<< "    vertex " << iter.key() << " is connected to edges "
				<< iter() << endl;
		}
		Pout<<endl;
	}

	return facePointEdges;
}


// Find next (triangle or cut) edge label coming from point prevVertI on
// prevEdgeI doing a right handed walk (i.e. following right wall).
// Note: normal is provided externally. Could be deducted from angle between
// two triangle edges but these could be in line.
Foam::label Foam::intersectedSurface::nextEdge
(
	const edgeSurface& eSurf,
	const Map<label>& visited,
	const label faceI,
	const vector& n,
	const Map<dynamicLabelList >& facePointEdges,
	const label prevEdgeI,
	const label prevVertI
)
{
	const pointField& points = eSurf.points();
	const edgeList& edges = eSurf.edges();
	const labelList& fEdges = eSurf.faceEdges()[faceI];


	// Edges connected to prevVertI
	const dynamicLabelList& connectedEdges = facePointEdges[prevVertI];

	if (connectedEdges.size() <= 1)
	{
		// Problem. Point not connected.
		{
			Pout<< "Writing face:" << faceI << " to face.obj" << endl;
			OFstream str("face.obj");
			writeOBJ(points, edges, fEdges, str);

			Pout<< "Writing connectedEdges edges to faceEdges.obj" << endl;
			writeLocalOBJ(points, edges, connectedEdges, "faceEdges.obj");
		}

		FatalErrorIn
		(
			"intersectedSurface::nextEdge(const pointField&, const edgeList&"
			", const vector&, Map<dynamicLabelList >, const label"
			", const label)"
		)   << "Problem: prevVertI:" << prevVertI << " on edge " << prevEdgeI
			<< " has less than 2 connected edges."
			<< " connectedEdges:" << connectedEdges << abort(FatalError);

		return -1;
	}

	if (connectedEdges.size() == 2)
	{
		// Simple case. Take other edge
		if (connectedEdges[0] == prevEdgeI)
		{
			if (debug & 2)
			{
				Pout<< "Stepped from edge:" << edges[prevEdgeI]
					<< " to edge:" << edges[connectedEdges[1]]
					<< " chosen from candidates:" << connectedEdges << endl;
			}
			return connectedEdges[1];
		}
		else
		{
			if (debug & 2)
			{
			   Pout<< "Stepped from edge:" << edges[prevEdgeI]
				   << " to edge:" << edges[connectedEdges[0]]
				   << " chosen from candidates:" << connectedEdges << endl;
			}
			return connectedEdges[0];
		}
	}


	// Multiple choices. Look at angle between edges.

	const edge& prevE = edges[prevEdgeI];

	// x-axis of coordinate system
	vector e0 = n ^ (points[prevE.otherVertex(prevVertI)] - points[prevVertI]);
	e0 /= mag(e0) + VSMALL;

	// Get y-axis of coordinate system
	vector e1 = n ^ e0;

	if (mag(mag(e1) - 1) > SMALL)
	{
		{
			Pout<< "Writing face:" << faceI << " to face.obj" << endl;
			OFstream str("face.obj");
			writeOBJ(points, edges, fEdges, str);

			Pout<< "Writing connectedEdges edges to faceEdges.obj" << endl;
			writeLocalOBJ(points, edges, connectedEdges, "faceEdges.obj");
		}

		FatalErrorIn("intersectedSurface::nextEdge")
			<< "Unnormalized normal e1:" << e1
			<< " formed from cross product of e0:" << e0 << " n:" << n
			<< abort(FatalError);
	}


	//
	// Determine maximum angle over all connected (unvisited) edges.
	//

	scalar maxAngle = -GREAT;
	label maxEdgeI = -1;

	forAll(connectedEdges, connI)
	{
		label edgeI = connectedEdges[connI];

		if (edgeI != prevEdgeI)
		{
			label stat = visited[edgeI];

			const edge& e = edges[edgeI];

			// Find out whether walk of edge from prevVert would be acceptible.
			if
			(
				stat == UNVISITED
			 || (
					stat == STARTTOEND
				 && prevVertI == e[1]
				)
			 || (
					stat == ENDTOSTART
				 && prevVertI == e[0]
				)
			)
			{
				// Calculate angle of edge with respect to base e0, e1
				vector vec =
					n ^ (points[e.otherVertex(prevVertI)] - points[prevVertI]);

				vec /= mag(vec) + VSMALL;

				scalar angle = pseudoAngle(e0, e1, vec);

				if (angle > maxAngle)
				{
					maxAngle = angle;
					maxEdgeI = edgeI;
				}
			}
		}
	}


	if (maxEdgeI == -1)
	{
		// No unvisited edge found
		{
			Pout<< "Writing face:" << faceI << " to face.obj" << endl;
			OFstream str("face.obj");
			writeOBJ(points, edges, fEdges, str);

			Pout<< "Writing connectedEdges edges to faceEdges.obj" << endl;
			writeLocalOBJ(points, edges, connectedEdges, "faceEdges.obj");
		}

		FatalErrorIn
		(
			"intersectedSurface::nextEdge(const pointField&, const edgeList&"
			", const Map<label>&, const vector&"
			", const Map<dynamicLabelList >&"
			", const label, const label"
		)   << "Trying to step from edge " << edges[prevEdgeI]
			<< ", vertex " << prevVertI
			<< " but cannot find 'unvisited' edges among candidates:"
			<< connectedEdges
			<< abort(FatalError);
	}

	if (debug & 2)
	{
		Pout<< "Stepped from edge:" << edges[prevEdgeI]
			<< " to edge:" << maxEdgeI << " angle:" << edges[maxEdgeI]
			<< " with angle:" << maxAngle
			<< endl;
	}

	return maxEdgeI;
}


// Create (polygonal) face by walking full circle starting from startVertI on
// startEdgeI.
// Uses nextEdge(..) to do the walking. Returns face. Updates visited with
// the labels of visited edges.
Foam::face Foam::intersectedSurface::walkFace
(
	const edgeSurface& eSurf,
	const label faceI,
	const vector& n,
	const Map<dynamicLabelList >& facePointEdges,

	const label startEdgeI,
	const label startVertI,

	Map<label>& visited
)
{
	const pointField& points = eSurf.points();
	const edgeList& edges = eSurf.edges();

	// Overestimate size of face
	face f(eSurf.faceEdges()[faceI].size());

	label fp = 0;

	label vertI = startVertI;
	label edgeI = startEdgeI;

	while(true)
	{
		const edge& e = edges[edgeI];

		if (debug & 2)
		{
			Pout<< "Now at:" << endl
				<< "    edge:" << edgeI << " vertices:" << e
				<< " positions:" << points[e.start()] << ' ' << points[e.end()]
				<< "    vertex:" << vertI << endl;
		}

		// Mark edge as visited
		if (e[0] == vertI)
		{
			visited[edgeI] |= STARTTOEND;
		}
		else
		{
			visited[edgeI] |= ENDTOSTART;
		}


		// Store face vertex
		f[fp++] = vertI;

		// step to other vertex
		vertI = e.otherVertex(vertI);

		if (vertI == startVertI)
		{
			break;
		}

		// step from vertex to next edge
		edgeI = nextEdge
		(
			eSurf,
			visited,
			faceI,
			n,
			facePointEdges,
			edgeI,
			vertI
		);
	}

	f.setSize(fp);

	return f;
}


void Foam::intersectedSurface::findNearestVisited
(
	const edgeSurface& eSurf,
	const label faceI,
	const Map<dynamicLabelList >& facePointEdges,
	const Map<label>& pointVisited,
	const point& pt,
	const label excludePointI,

	label& minVertI,
	scalar& minDist
)
{
	minVertI = -1;
	minDist = GREAT;

	forAllConstIter(Map<label>, pointVisited, iter)
	{
		label pointI = iter.key();

		if (pointI != excludePointI)
		{
			label nVisits = iter();

			if (nVisits == 2*facePointEdges[pointI].size())
			{
				// Fully visited (i.e. both sides of all edges)
				scalar dist = mag(eSurf.points()[pointI] - pt);

				if (dist < minDist)
				{
					minDist = dist;
					minVertI = pointI;
				}
			}
		}
	}

	if (minVertI == -1)
	{
		const labelList& fEdges = eSurf.faceEdges()[faceI];

		SeriousErrorIn("intersectedSurface::findNearestVisited")
			<< "Dumping face edges to faceEdges.obj" << endl;

		writeLocalOBJ(eSurf.points(), eSurf.edges(), fEdges, "faceEdges.obj");

		FatalErrorIn("intersectedSurface::findNearestVisited")
			<< "No fully visited edge found for pt " << pt
			<< abort(FatalError);
	}
}


// Sometimes there are bunches of edges that are not connected to the
// other edges in the face. This routine tries to connect the loose
// edges up to the 'proper' edges by adding two extra edges between a
// properly connected edge and an unconnected one. Since at this level the
// adressing is purely in form of points and a cloud of edges this can
// be easily done.
// (edges are to existing points. Don't want to introduce new vertices here
// since then also neighbouring face would have to be split)
Foam::faceList Foam::intersectedSurface::resplitFace
(
	const triSurface& surf,
	const label faceI,
	const Map<dynamicLabelList >& facePointEdges,
	const Map<label>& visited,
	edgeSurface& eSurf
)
{
	{
		// Dump face for debugging.
		Pout<< "Writing face:" << faceI << " to face.obj" << endl;
		OFstream str("face.obj");
		writeOBJ(eSurf.points(), eSurf.edges(), eSurf.faceEdges()[faceI], str);
	}


	// Count the number of times point has been visited so we
	// can compare number to facePointEdges.
	Map<label> pointVisited(2*facePointEdges.size());


	{
		OFstream str0("visitedNone.obj");
		OFstream str1("visitedOnce.obj");
		OFstream str2("visitedTwice.obj");
		forAll(eSurf.points(), pointI)
		{
			const point& pt = eSurf.points()[pointI];

			str0 << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
			str1 << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
			str2 << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
		}


	forAllConstIter(Map<label>, visited, iter)
	{
		label edgeI = iter.key();

		const edge& e = eSurf.edges()[edgeI];

		label stat = iter();

		if (stat == STARTTOEND || stat == ENDTOSTART)
		{
			incCount(pointVisited, e[0], 1);
			incCount(pointVisited, e[1], 1);

			str1 << "l " << e[0]+1 << ' ' << e[1]+1 << nl;
		}
		else if (stat == BOTH)
		{
			incCount(pointVisited, e[0], 2);
			incCount(pointVisited, e[1], 2);

			str2 << "l " << e[0]+1 << ' ' << e[1]+1 << nl;
		}
		else if (stat == UNVISITED)
		{
			incCount(pointVisited, e[0], 0);
			incCount(pointVisited, e[1], 0);

			str0 << "l " << e[0]+1 << ' ' << e[1]+1 << nl;
		}
	}
	}


	{
		forAllConstIter(Map<label>, pointVisited, iter)
		{
			label pointI = iter.key();

			label nVisits = iter();

			Pout<< "point:" << pointI << "  nVisited:" << nVisits
				<< "  pointEdges:" << facePointEdges[pointI].size() << endl;
		}
	}


	// Find nearest point pair where one is not fully visited and
	// the other is.

	label visitedVert0 = -1;
	label unvisitedVert0 = -1;

	{
		scalar minDist = GREAT;

		forAllConstIter(Map<dynamicLabelList >, facePointEdges, iter)
		{
			label pointI = iter.key();

			label nVisits = pointVisited[pointI];

			const dynamicLabelList& pEdges = iter();

			if (nVisits < 2*pEdges.size())
			{
				// Not fully visited. Find nearest fully visited.

				scalar nearDist;
				label nearVertI;

				findNearestVisited
				(
					eSurf,
					faceI,
					facePointEdges,
					pointVisited,
					eSurf.points()[pointI],
					-1,                         // Do not exclude vertex
					nearVertI,
					nearDist
				);


				if (nearDist < minDist)
				{
					minDist = nearDist;
					visitedVert0 = nearVertI;
					unvisitedVert0 = pointI;
				}
			}
		}
	}


	// Find second intersection.
	label visitedVert1 = -1;
	label unvisitedVert1 = -1;
	{
		scalar minDist = GREAT;

		forAllConstIter(Map<dynamicLabelList >, facePointEdges, iter)
		{
			label pointI = iter.key();

			if (pointI != unvisitedVert0)
			{
				label nVisits = pointVisited[pointI];

				const dynamicLabelList& pEdges = iter();

				if (nVisits < 2*pEdges.size())
				{
					// Not fully visited. Find nearest fully visited.

					scalar nearDist;
					label nearVertI;

					findNearestVisited
					(
					    eSurf,
					    faceI,
					    facePointEdges,
					    pointVisited,
					    eSurf.points()[pointI],
					    visitedVert0,           // vertex to exclude
					    nearVertI,
					    nearDist
					);


					if (nearDist < minDist)
					{
					    minDist = nearDist;
					    visitedVert1 = nearVertI;
					    unvisitedVert1 = pointI;
					}
				}
			}
		}
	}


	Pout<< "resplitFace : adding intersection from " << visitedVert0
		<< " to " << unvisitedVert0 << endl
		<< " and from " << visitedVert1
		<< " to " << unvisitedVert1 << endl;


//    // Add the new intersection edges to the edgeSurface
//    edgeList additionalEdges(2);
//    additionalEdges[0] = edge(visitedVert0, unvisitedVert0);
//    additionalEdges[1] = edge(visitedVert1, unvisitedVert1);

	// Add the new intersection edges to the edgeSurface
	edgeList additionalEdges(1);
	additionalEdges[0] = edge(visitedVert0, unvisitedVert0);

	eSurf.addIntersectionEdges(faceI, additionalEdges);

	fileName newFName("face_" + Foam::name(faceI) + "_newEdges.obj");
	Pout<< "Dumping face:" << faceI << " to " << newFName << endl;
	writeLocalOBJ
	(
		eSurf.points(),
		eSurf.edges(),
		eSurf.faceEdges()[faceI],
		newFName
	);

	// Retry splitFace. Use recursion since is rare situation.
	return splitFace(surf, faceI, eSurf);
}


Foam::faceList Foam::intersectedSurface::splitFace
(
	const triSurface& surf,
	const label faceI,
	edgeSurface& eSurf
)
{
	// Alias
	const pointField& points = eSurf.points();
	const edgeList& edges = eSurf.edges();
	const labelList& fEdges = eSurf.faceEdges()[faceI];

	// Create local (for the face only) point-edge connectivity.
	Map<dynamicLabelList > facePointEdges
	(
		calcPointEdgeAddressing
		(
			eSurf,
			faceI
		)
	);

	// Order in which edges have been walked. Initialize outside edges.
	Map<label> visited(fEdges.size()*2);

	forAll(fEdges, i)
	{
		label edgeI = fEdges[i];

		if (eSurf.isSurfaceEdge(edgeI))
		{
			// Edge is edge from original surface so an outside edge for
			// the current face.
			label surfEdgeI = eSurf.parentEdge(edgeI);

			label owner = surf.edgeOwner()[surfEdgeI];

			if
			(
				owner == faceI
			 || sameEdgeOrder
				(
					surf.localFaces()[owner],
					surf.localFaces()[faceI]
				)
			)
			{
				// Edge is in same order as current face.
				// Mark off the opposite order.
				visited.insert(edgeI, ENDTOSTART);
			}
			else
			{
				// Edge is in same order as current face.
				// Mark off the opposite order.
				visited.insert(edgeI, STARTTOEND);
			}
		}
		else
		{
			visited.insert(edgeI, UNVISITED);
		}
	}



	// Storage for faces
	DynamicList<face> faces(fEdges.size());

	while (true)
	{
		// Find starting edge:
		// - unvisted triangle edge
		// - once visited intersection edge
		// Give priority to triangle edges.
		label startEdgeI = -1;
		label startVertI = -1;

		forAll(fEdges, i)
		{
			label edgeI = fEdges[i];

			const edge& e = edges[edgeI];

			label stat = visited[edgeI];

			if (stat == STARTTOEND)
			{
				startEdgeI = edgeI;
				startVertI = e[1];

				if (eSurf.isSurfaceEdge(edgeI))
				{
					break;
				}
			}
			else if (stat == ENDTOSTART)
			{
				startEdgeI = edgeI;
				startVertI = e[0];

				if (eSurf.isSurfaceEdge(edgeI))
				{
					break;
				}
			}
		}

		if (startEdgeI == -1)
		{
			break;
		}

		//Pout<< "splitFace: starting walk from edge:" << startEdgeI
		//    << ' ' << edges[startEdgeI] << " vertex:" << startVertI << endl;

		//// Print current visited state.
		//printVisit(eSurf.edges(), fEdges, visited);

		//{
		//    Pout<< "Writing face:" << faceI << " to face.obj" << endl;
		//    OFstream str("face.obj");
		//    writeOBJ(eSurf.points(), eSurf.edges(), fEdges, str);
		//}

		faces.append
		(
			walkFace
			(
				eSurf,
				faceI,
				surf.faceNormals()[faceI],
				facePointEdges,

				startEdgeI,
				startVertI,

				visited
			)
		);
	}


	// Check if any unvisited edges left.
	forAll(fEdges, i)
	{
		label edgeI = fEdges[i];

		label stat = visited[edgeI];

		if (eSurf.isSurfaceEdge(edgeI) && stat != BOTH)
		{
			SeriousErrorIn("Foam::intersectedSurface::splitFace")
				<< "Dumping face edges to faceEdges.obj" << endl;

			writeLocalOBJ(points, edges, fEdges, "faceEdges.obj");

			FatalErrorIn("intersectedSurface::splitFace")
			   << "Problem: edge " << edgeI << " vertices "
				<< edges[edgeI] << " on face " << faceI
				<< " has visited status " << stat << " from a "
				 << "righthanded walk along all"
				<< " of the triangle edges. Are the original surfaces"
				<< " closed and non-intersecting?"
				<< abort(FatalError);
		}
		else if (stat != BOTH)
		{
			//{
			//    Pout<< "Dumping faces so far to faces.obj" << nl
			//        << faces << endl;
			//
			//    OFstream str("faces.obj");
			//
			//    forAll(faces, i)
			//    {
			//        writeOBJ(points, faces[i], str);
			//    }
			//}

			Pout<< "** Resplitting **" << endl;

			// Redo face after introducing extra edge. Edge introduced
			// should be one nearest to any fully visited edge.
			return resplitFace
			(
				surf,
				faceI,
				facePointEdges,
				visited,
				eSurf
			);
		}
	}


	// See if normal needs flipping.
	faces.shrink();

	vector n = faces[0].normal(eSurf.points());

	if ((n & surf.faceNormals()[faceI]) < 0)
	{
		forAll(faces, i)
		{
			reverse(faces[i]);
		}
	}

	return faces;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Null constructor
Foam::intersectedSurface::intersectedSurface()
:
	triSurface(),
	intersectionEdges_(0),
	faceMap_(0),
	nSurfacePoints_(0)
{}


// Construct from components
Foam::intersectedSurface::intersectedSurface(const triSurface& surf)
:
	triSurface(surf),
	intersectionEdges_(0),
	faceMap_(0),
	nSurfacePoints_(0)
{}


// Construct from surface and intersection
Foam::intersectedSurface::intersectedSurface
(
	const triSurface& surf,
	const bool isFirstSurface,
	const surfaceIntersection& inter
)
:
	triSurface(),
	intersectionEdges_(0),
	faceMap_(0),
	nSurfacePoints_(surf.nPoints())
{
	if (inter.cutPoints().empty() && inter.cutEdges().empty())
	{
		// No intersection. Make straight copy.
		triSurface::operator=(surf);

		// Identity for face map
		faceMap_.setSize(size());

		forAll(faceMap_, faceI)
		{
			faceMap_[faceI] = faceI;
		}
		return;
	}


	// Calculate face-edge addressing for surface + intersection.
	edgeSurface eSurf(surf, isFirstSurface, inter);

	// Update point information for any removed points. (when are they removed?
	// - but better make sure)
	nSurfacePoints_ = eSurf.nSurfacePoints();

	// Now we have full points, edges and edge addressing for surf. Start
	// extracting faces and triangulate them.

	// Storage for new triangles of surface.
	DynamicList<labelledTri> newTris(eSurf.edges().size()/2);

	// Start in newTris for decomposed face.
	labelList startTriI(surf.size(), 0);

	forAll(surf, faceI)
	{
		startTriI[faceI] = newTris.size();

		if (eSurf.faceEdges()[faceI].size() != surf.faceEdges()[faceI].size())
		{
			// Face has been cut by intersection.
			// Cut face into multiple subfaces. Use faceEdge information
			// from edgeSurface only. (original triSurface 'surf' is used for
			// faceNormals and regionnumber only)
			faceList newFaces
			(
				splitFace
				(
					surf,
					faceI,              // current triangle
					eSurf               // face-edge description of surface
					                    // + intersection
				)
			);
			forAll(newFaces, newFaceI)
			{
				const face& newF = newFaces[newFaceI];

//                {
//                    fileName fName
//                    (
//                        "face_"
//                      + Foam::name(faceI)
//                      + "_subFace_"
//                      + Foam::name(newFaceI)
//                      + ".obj"
//                    );
//                    Pout<< "Writing original face:" << faceI << " subFace:"
//                        << newFaceI << " to " << fName << endl;
//
//                    OFstream str(fName);
//
//                    forAll(newF, fp)
//                    {
//                        meshTools::writeOBJ(str, eSurf.points()[newF[fp]]);
//                    }
//                    str << 'l';
//                    forAll(newF, fp)
//                    {
//                        str << ' ' << fp+1;
//                    }
//                    str<< " 1" << nl;
//                }


				const vector& n = surf.faceNormals()[faceI];
				const label region = surf[faceI].region();

				faceTriangulation tris(eSurf.points(), newF, n);

				forAll(tris, triI)
				{
					const triFace& t = tris[triI];

					forAll(t, i)
					{
					    if (t[i] < 0 || t[i] >= eSurf.points().size())
					    {
					        FatalErrorIn
					        (
					            "intersectedSurface::intersectedSurface"
					        )   << "Face triangulation of face " << faceI
					            << " uses points outside range 0.."
					            << eSurf.points().size()-1 << endl
					            << "Triangulation:"
					            << tris << abort(FatalError);
					    }
					}

					newTris.append(labelledTri(t[0], t[1], t[2], region));
				}
			}
		}
		else
		{
			// Face has not been cut at all. No need to renumber vertices since
			// eSurf keeps surface vertices first.
			newTris.append(surf.localFaces()[faceI]);
		}
	}

	newTris.shrink();


	// Construct triSurface. Note that addressing will be compact since
	// edgeSurface is compact.
	triSurface::operator=
	(
		triSurface
		(
			newTris,
			surf.patches(),
			eSurf.points()
		)
	);

	// Construct mapping back into original surface
	faceMap_.setSize(size());

	for (label faceI = 0; faceI < surf.size()-1; faceI++)
	{
		for (label triI = startTriI[faceI]; triI < startTriI[faceI+1]; triI++)
		{
			faceMap_[triI] = faceI;
		}
	}
	for (label triI = startTriI[surf.size()-1]; triI < size(); triI++)
	{
		faceMap_[triI] = surf.size()-1;
	}


	// Find edges on *this which originate from 'cuts'. (i.e. newEdgeI >=
	// nSurfaceEdges) Renumber edges into local triSurface numbering.

	intersectionEdges_.setSize(eSurf.edges().size() - eSurf.nSurfaceEdges());

	label intersectionEdgeI = 0;

	for
	(
		label edgeI = eSurf.nSurfaceEdges();
		edgeI < eSurf.edges().size();
		edgeI++
	)
	{
		// Get edge vertices in triSurface local numbering
		const edge& e = eSurf.edges()[edgeI];
		label surfStartI = meshPointMap()[e.start()];
		label surfEndI = meshPointMap()[e.end()];

		// Find edge connected to surfStartI which also uses surfEndI.
		label surfEdgeI = -1;

		const labelList& pEdges = pointEdges()[surfStartI];

		forAll(pEdges, i)
		{
			const edge& surfE = edges()[pEdges[i]];

			// Edge already connected to surfStart for sure. See if also
			// connects to surfEnd
			if (surfE.start() == surfEndI || surfE.end() == surfEndI)
			{
				surfEdgeI = pEdges[i];

				break;
			}
		}

		if (surfEdgeI != -1)
		{
			intersectionEdges_[intersectionEdgeI++] = surfEdgeI;
		}
		else
		{
			FatalErrorIn("intersectedSurface::intersectedSurface")
				<< "Cannot find edge among candidates " << pEdges
				<< " which uses points " << surfStartI
				<< " and " << surfEndI
				<< abort(FatalError);
		}
	}
}


// ************************************************************************* //
