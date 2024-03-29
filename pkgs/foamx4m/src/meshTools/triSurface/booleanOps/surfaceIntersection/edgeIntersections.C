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

\*---------------------------------------------------------------------------*/

#include "edgeIntersections.H"
#include "triSurfaceSearch.H"
#include "labelPairLookup.H"
#include "OFstream.H"
#include "HashSet.H"
#include "triSurface.H"
#include "pointIndexHit.H"
#include "treeDataTriSurface.H"
#include "indexedOctree.H"
#include "meshTools.H"
#include "plane.H"
#include "Random.H"
#include "mathematicalConstants.H"
#include "treeBoundBox.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::edgeIntersections, 0);

Foam::scalar Foam::edgeIntersections::alignedCos_ =
	Foam::cos(89.0 * Foam::mathematicalConstant::pi/180.0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::edgeIntersections::checkEdges(const triSurface& surf)
{
	const pointField& localPoints = surf.localPoints();
	const edgeList& edges = surf.edges();
	const labelListList& edgeFaces = surf.edgeFaces();

	treeBoundBox bb(localPoints);

	scalar minSize = SMALL * bb.minDim();

	forAll(edges, edgeI)
	{
		const edge& e = edges[edgeI];

		scalar eMag = e.mag(localPoints);

		if (eMag < minSize)
		{
			WarningIn
			(
				"Foam::edgeIntersections::checkEdges(const triSurface& surf)"
			)   << "Edge " << edgeI << " vertices " << e
				<< " coords:" << localPoints[e[0]] << ' '
				<< localPoints[e[1]] << " is very small compared to bounding"
				<< " box dimensions " << bb << endl
				<< "This might lead to problems in intersection"
				<< endl;
		}

		if (edgeFaces[edgeI].size() == 1)
		{
			WarningIn
			(
				"Foam::edgeIntersections::checkEdges(const triSurface& surf)"
			)   << "Edge " << edgeI << " vertices " << e
				<< " coords:" << localPoints[e[0]] << ' '
				<< localPoints[e[1]] << " has only one face connected to it:"
				<< edgeFaces[edgeI] << endl
				<< "This might lead to problems in intersection"
				<< endl;
		}
	}
}


// Update intersections for selected edges.
void Foam::edgeIntersections::intersectEdges
(
	const triSurface& surf1,
	const pointField& points1,          // surf1 meshPoints (not localPoints!)
	const triSurfaceSearch& querySurf2,
	const scalarField& surf1PointTol,   // surf1 tolerance per point
	const labelList& edgeLabels
)
{
	const triSurface& surf2 = querySurf2.surface();
	const vectorField& normals2 = surf2.faceNormals();

	const labelList& meshPoints = surf1.meshPoints();

	if (debug)
	{
		Pout<< "Calculating intersection of " << edgeLabels.size() << " edges"
			<< " out of " << surf1.nEdges() << " with " << surf2.size()
			<< " triangles ..." << endl;
	}

	// Construct octree.
	const indexedOctree<treeDataTriSurface>& tree = querySurf2.tree();


	label nHits = 0;


	// Go through all edges, calculate intersections
	forAll(edgeLabels, i)
	{
		label edgeI = edgeLabels[i];

		if (debug && (i % 1000 == 0))
		{
			Pout<< "Intersecting edge " << edgeI << " with surface" << endl;
		}

		const edge& e = surf1.edges()[edgeI];

		const point& pStart = points1[meshPoints[e.start()]];
		const point& pEnd = points1[meshPoints[e.end()]];

		const vector eVec(pEnd - pStart);
		const scalar eMag = mag(eVec);
		const vector n(eVec/(eMag + VSMALL));

		// Smallish length for intersection calculation errors.
		const point tolVec = 1e-6*eVec;

		// Start tracking somewhat before pStart and upto somewhat after p1.
		// Note that tolerances here are smaller than those used to classify
		// hit below.
		// This will cause this hit to be marked as degenerate and resolved
		// later on.
		point p0 = pStart - 0.5*surf1PointTol[e[0]]*n;
		const point p1 = pEnd + 0.5*surf1PointTol[e[1]]*n;
		const scalar maxS = mag(p1 - pStart);

		// Get all intersections of the edge with the surface

		DynamicList<pointIndexHit> currentIntersections(100);
		dynamicLabelList currentIntersectionTypes(100);

		while (true)
		{
			pointIndexHit pHit = tree.findLine(p0, p1);

			if (pHit.hit())
			{
				nHits++;

				currentIntersections.append(pHit);

				// Classify point on surface1 edge.
				label edgeEnd = -1;

				if (mag(pHit.hitPoint() - pStart) < surf1PointTol[e[0]])
				{
					edgeEnd = 0;
				}
				else if (mag(pHit.hitPoint() - pEnd) < surf1PointTol[e[1]])
				{
					edgeEnd = 1;
				}
				else if (mag(n & normals2[pHit.index()]) < alignedCos_)
				{
					Pout<< "Flat angle edge:" << edgeI
					    << " face:" << pHit.index()
					    << " cos:" << mag(n & normals2[pHit.index()])
					    << endl;
					edgeEnd = 2;
				}

				currentIntersectionTypes.append(edgeEnd);

				if (edgeEnd == 1)
				{
					// Close to end
					break;
				}
				else
				{
					// Continue tracking. Shift by a small amount.
					p0 = pHit.hitPoint() + tolVec;

					if (((p0-pStart) & n) >= maxS)
					{
					    break;
					}
				}
			}
			else
			{
				// No hit.
				break;
			}
		}


		// Done current edge. Transfer all data into *this
		operator[](edgeI).transfer(currentIntersections);
		classification_[edgeI].transfer(currentIntersectionTypes);
	}

	if (debug)
	{
		Pout<< "Found " << nHits << " intersections of edges with surface ..."
			<< endl;
	}

}


// If edgeI intersections are close to endpoint of edge shift endpoints
// slightly along edge
// Updates
// - points1 with new endpoint position
// - affectedEdges with all edges affected by moving the point
// Returns true if changed anything.
bool Foam::edgeIntersections::inlinePerturb
(
	const triSurface& surf1,
	const scalarField& surf1PointTol,   // surf1 tolerance per point
	const label edgeI,
	Random& rndGen,
	pointField& points1,
	boolList& affectedEdges
) const
{
	bool hasPerturbed = false;

	// Check if edge close to endpoint. Note that we only have to check
	// the intersection closest to the edge endpoints (i.e. first and last in
	// edgeEdgs)

	const labelList& edgeEnds = classification_[edgeI];

	if (edgeEnds.size())
	{
		bool perturbStart = false;
		bool perturbEnd = false;

		// Check first intersection.
		if (edgeEnds[0] == 0)
		{
			perturbStart = true;
		}

		if (edgeEnds[edgeEnds.size()-1] == 1)
		{
			perturbEnd = true;
		}


		if (perturbStart || perturbEnd)
		{
			const edge& e = surf1.edges()[edgeI];

			label v0 = surf1.meshPoints()[e[0]];
			label v1 = surf1.meshPoints()[e[1]];

			vector eVec(points1[v1] - points1[v0]);
			vector n = eVec/mag(eVec);

			if (perturbStart)
			{
				// Perturb with something (hopefully) larger than tolerance.
				scalar t = 4.0*(rndGen.scalar01() - 0.5);
				points1[v0] += t*surf1PointTol[e[0]]*n;

				const labelList& pEdges = surf1.pointEdges()[e[0]];

				forAll(pEdges, i)
				{
					affectedEdges[pEdges[i]] = true;
				}
			}
			if (perturbEnd)
			{
				// Perturb with something larger than tolerance.
				scalar t = 4.0*(rndGen.scalar01() - 0.5);
				points1[v1] += t*surf1PointTol[e[1]]*n;

				const labelList& pEdges = surf1.pointEdges()[e[1]];

				forAll(pEdges, i)
				{
					affectedEdges[pEdges[i]] = true;
				}
			}
			hasPerturbed = true;
		}
	}

	return hasPerturbed;
}


// Perturb single edge endpoint when perpendicular to face
bool Foam::edgeIntersections::rotatePerturb
(
	const triSurface& surf1,
	const scalarField& surf1PointTol,   // surf1 tolerance per point
	const label edgeI,

	Random& rndGen,
	pointField& points1,
	boolList& affectedEdges
) const
{
	const labelList& meshPoints = surf1.meshPoints();

	const labelList& edgeEnds = classification_[edgeI];

	bool hasPerturbed = false;

	forAll(edgeEnds, i)
	{
		if (edgeEnds[i] == 2)
		{
			const edge& e = surf1.edges()[edgeI];

			// Endpoint to modify. Choose either start or end.
			label pointI = e[rndGen.bit()];
			//label pointI = e[0];

			// Generate random vector slightly larger than tolerance.
			vector rndVec = rndGen.vector01() - vector(0.5, 0.5, 0.5);

			// Make sure rndVec only perp to edge
			vector n(points1[meshPoints[e[1]]] - points1[meshPoints[e[0]]]);
			scalar magN = mag(n) + VSMALL;
			n /= magN;

			rndVec -= n*(n & rndVec);

			// Normalize
			rndVec /= mag(rndVec) + VSMALL;

			// Scale to be moved by tolerance.
			rndVec *= 0.01*magN;

			Pout<< "rotating: shifting endpoint " << meshPoints[pointI]
				<< " of edge:" << edgeI << " verts:"
				<< points1[meshPoints[e[0]]] << ' '
				<< points1[meshPoints[e[1]]]
				<< " by " << rndVec
				<< " tol:" << surf1PointTol[pointI] << endl;

			points1[meshPoints[pointI]] += rndVec;

			// Mark edges affected by change to point
			const labelList& pEdges = surf1.pointEdges()[pointI];

			forAll(pEdges, i)
			{
				affectedEdges[pEdges[i]] = true;
			}

			hasPerturbed = true;

			// Enough done for current edge; no need to test other intersections
			// of this edge.
			break;
		}
	}

	return hasPerturbed;
}


// Perturb edge when close to face
bool Foam::edgeIntersections::offsetPerturb
(
	const triSurface& surf1,
	const triSurface& surf2,
	const label edgeI,

	Random& rndGen,
	pointField& points1,
	boolList& affectedEdges
) const
{
	const labelList& meshPoints = surf1.meshPoints();

	const edge& e = surf1.edges()[edgeI];

	const List<pointIndexHit>& hits = operator[](edgeI);


	bool hasPerturbed = false;

	// For all hits on edge
	forAll(hits, i)
	{
		const pointIndexHit& pHit = hits[i];

		// Classify point on face of surface2
		label surf2FaceI = pHit.index();

		const labelledTri& f2 = surf2.localFaces()[surf2FaceI];

		const pointField& surf2Pts = surf2.localPoints();

		label nearType;
		label nearLabel;

		triPointRef tri
		(
			surf2Pts[f2[0]],
			surf2Pts[f2[1]],
			surf2Pts[f2[2]]
		);

		point ctr = tri.centre();

		// Get measure for tolerance.
		scalar tolDim = 0.001*mag(tri.a() - ctr);

		tri.classify(pHit.hitPoint(), tolDim, nearType, nearLabel);

		if (nearType == triPointRef::POINT || nearType == triPointRef::EDGE)
		{
			// Shift edge towards tri centre
			vector offset = 0.01*rndGen.scalar01()*(ctr - pHit.hitPoint());

			// shift e[0]
			points1[meshPoints[e[0]]] += offset;

			// Mark edges affected by change to e0
			const labelList& pEdges0 = surf1.pointEdges()[e[0]];

			forAll(pEdges0, i)
			{
				affectedEdges[pEdges0[i]] = true;
			}

			// shift e[1]
			points1[meshPoints[e[1]]] += offset;

			// Mark edges affected by change to e1
			const labelList& pEdges1 = surf1.pointEdges()[e[1]];

			forAll(pEdges1, i)
			{
				affectedEdges[pEdges1[i]] = true;
			}

			hasPerturbed = true;

			// No need to test any other hits on this edge
			break;
		}
	}

	return hasPerturbed;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct null
Foam::edgeIntersections::edgeIntersections()
:
	List<List<pointIndexHit> >(),
	classification_()
{}


// Construct from surface and tolerance
Foam::edgeIntersections::edgeIntersections
(
	const triSurface& surf1,
	const triSurfaceSearch& query2,
	const scalarField& surf1PointTol
)
:
	List<List<pointIndexHit> >(surf1.nEdges()),
	classification_(surf1.nEdges())
{
	checkEdges(surf1);
	checkEdges(query2.surface());

	// Current set of edges to test
	labelList edgesToTest(surf1.nEdges());

	// Start off with all edges
	forAll(edgesToTest, i)
	{
		edgesToTest[i] = i;
	}

	// Determine intersections for edgesToTest
	intersectEdges
	(
		surf1,
		surf1.points(), // surf1 meshPoints (not localPoints!)
		query2,
		surf1PointTol,
		edgesToTest
	);
}


// Construct from components
Foam::edgeIntersections::edgeIntersections
(
	const List<List<pointIndexHit> >& intersections,
	const labelListList& classification
)
:
	List<List<pointIndexHit> >(intersections),
	classification_(classification)
{}


// * * * * * * * * * * * * * * * Static Functions  * * * * * * * * * * * * * //

Foam::scalarField Foam::edgeIntersections::minEdgeLength(const triSurface& surf)
{
	const pointField& localPoints = surf.localPoints();
	const labelListList& pointEdges = surf.pointEdges();
	const edgeList& edges = surf.edges();

	scalarField minLen(localPoints.size());

	forAll(minLen, pointI)
	{
		const labelList& pEdges = pointEdges[pointI];

		scalar minDist = GREAT;

		forAll(pEdges, i)
		{
			minDist = min(minDist, edges[pEdges[i]].mag(localPoints));
		}

		minLen[pointI] = minDist;
	}
	return minLen;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::label Foam::edgeIntersections::removeDegenerates
(
	const label nIters,
	const triSurface& surf1,
	const triSurfaceSearch& query2,
	const scalarField& surf1PointTol,
	pointField& points1
)
{
	const triSurface& surf2 = query2.surface();

	Random rndGen(356574);

	// Current set of edges to (re)test
	labelList edgesToTest(surf1.nEdges());

	// Start off with all edges
	forAll(edgesToTest, i)
	{
		edgesToTest[i] = i;
	}


	label iter = 0;

	for (; iter < nIters; iter++)
	{
		// Go through all edges to (re)test and perturb points if they are
		// degenerate hits. Mark off edges that need to be recalculated.

		boolList affectedEdges(surf1.nEdges(), false);
		label nShifted = 0;
		label nRotated = 0;
		label nOffset = 0;

		forAll(edgesToTest, i)
		{
			label edgeI = edgesToTest[i];

			// If edge not already marked for retesting
			if (!affectedEdges[edgeI])
			{
				// 1. Check edges close to endpoint and perturb if nessecary.

				bool shiftedEdgeEndPoints =
					inlinePerturb
					(
					    surf1,
					    surf1PointTol,
					    edgeI,
					    rndGen,
					    points1,
					    affectedEdges
					);

				nShifted += (shiftedEdgeEndPoints ? 1 : 0);

				if (!shiftedEdgeEndPoints)
				{
					bool rotatedEdge =
					    rotatePerturb
					    (
					        surf1,
					        surf1PointTol,
					        edgeI,
					        rndGen,
					        points1,
					        affectedEdges
					    );

					nRotated += (rotatedEdge ? 1 : 0);

					if (!rotatedEdge)
					{
					    // 2. we're sure now that the edge actually pierces the
					    // face. Now check the face for intersections close its
					    // points/edges

					    bool offsetEdgePoints =
					        offsetPerturb
					        (
					            surf1,
					            surf2,
					            edgeI,
					            rndGen,
					            points1,
					            affectedEdges
					        );

					    nOffset += (offsetEdgePoints ? 1 : 0);
					}
				}
			}
		}

		if (debug)
		{
			Pout<< "Edges to test : " << nl
				<< "    total:" << edgesToTest.size() << nl
				<< "    resolved by:" << nl
				<< "        shifting   : " << nShifted << nl
				<< "        rotating   : " << nRotated << nl
				<< "        offsetting : " << nOffset << nl
				<< endl;
		}


		if (nShifted == 0 && nRotated == 0 && nOffset == 0)
		{
			// Nothing changed in current iteration. Current hit pattern is
			// without any degenerates.
			break;
		}

		// Repack affected edges
		labelList newEdgesToTest(surf1.nEdges());
		label newEdgeI = 0;

		forAll(affectedEdges, edgeI)
		{
			if (affectedEdges[edgeI])
			{
				newEdgesToTest[newEdgeI++] = edgeI;
			}
		}
		newEdgesToTest.setSize(newEdgeI);

		if (debug)
		{
			Pout<< "Edges to test:" << nl
				<< "    was : " << edgesToTest.size() << nl
				<< "    is  : " << newEdgesToTest.size() << nl
				<< endl;
		}

		// Transfer and test.
		edgesToTest.transfer(newEdgesToTest);

		if (edgesToTest.empty())
		{
			FatalErrorIn("perturb") << "oops" << abort(FatalError);
		}

		// Re intersect moved edges.
		intersectEdges
		(
			surf1,
			points1,          // surf1 meshPoints (not localPoints!)
			query2,
			surf1PointTol,
			edgesToTest
		);
	}

	return iter;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //


// ************************************************************************* //
