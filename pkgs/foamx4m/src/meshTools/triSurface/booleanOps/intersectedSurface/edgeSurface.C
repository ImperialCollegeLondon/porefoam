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

#include "edgeSurface.H"
#include "triSurface.H"
#include "surfaceIntersection.H"
#include "meshTools.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::edgeSurface, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Write whole pointField and edges to stream
void Foam::edgeSurface::writeOBJ
(
	const pointField& points,
	const edgeList& edges,
	Ostream& os
)
{
	forAll(points, pointI)
	{
		const point& pt = points[pointI];

		os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
	}
	forAll(edges, edgeI)
	{
		const edge& e = edges[edgeI];

		os << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
	}
}


// Write whole pointField and selected edges to stream
void Foam::edgeSurface::writeOBJ
(
	const pointField& points,
	const edgeList& edges,
	const labelList& edgeLabels,
	Ostream& os
)
{
	forAll(points, pointI)
	{
		const point& pt = points[pointI];

		os << "v " << pt.x() << ' ' << pt.y() << ' ' << pt.z() << endl;
	}
	forAll(edgeLabels, i)
	{
		const edge& e = edges[edgeLabels[i]];

		os << "l " << e.start()+1 << ' ' << e.end()+1 << endl;
	}
}


// Pointedges in edgeSurface indices only.
void Foam::edgeSurface::calcPointEdges()
{
	pointEdges_.setSize(points_.size());

	labelList pointNEdges(points_.size(), 0);

	forAll(edges_, edgeI)
	{
		const edge& e = edges_[edgeI];

		pointNEdges[e[0]]++;
		pointNEdges[e[1]]++;
	}

	forAll(pointEdges_, pointI)
	{
		pointEdges_[pointI].setSize(pointNEdges[pointI]);
	}

	pointNEdges = 0;

	forAll(edges_, edgeI)
	{
		const edge& e = edges_[edgeI];

		labelList& pEdges0 = pointEdges_[e[0]];
		pEdges0[pointNEdges[e[0]]++] = edgeI;

		labelList& pEdges1 = pointEdges_[e[1]];
		pEdges1[pointNEdges[e[1]]++] = edgeI;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from surface and intersection description
Foam::edgeSurface::edgeSurface
(
	const triSurface& surf,
	const bool isFirstSurface,
	const surfaceIntersection& inter
)
:
	points_(surf.nPoints() + inter.cutPoints().size()),
	nSurfacePoints_(surf.nPoints()),
	edges_(),
	nSurfaceEdges_(surf.nEdges()),
	parentEdges_(0),
	faceEdges_(surf.size()),
	pointEdges_(points_.size())
{
	// Copy points (surface ones first)
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	label pointI = 0;

	const pointField& surfPoints = surf.localPoints();

	forAll(surfPoints, i)
	{
		points_[pointI++] = surfPoints[i];
	}

	const pointField& cutPoints = inter.cutPoints();

	forAll(cutPoints, i)
	{
		points_[pointI++] = cutPoints[i];
	}


	// Copy edges (surface ones first)
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	DynamicList<edge> allEdges(surf.nEdges() + inter.cutEdges().size());
	dynamicLabelList allParentEdges(surf.nEdges());
	List<dynamicLabelList > allFaceEdges(surf.size());


	// Copy surface edges (can be split!)

	const edgeList& surfEdges = surf.edges();

	forAll(surfEdges, edgeI)
	{
		const edge& e = surfEdges[edgeI];

		// Get additional vertices for this edge.
		const labelList& extraVerts = inter.edgeCuts(isFirstSurface)[edgeI];

		// Store current top of allEdges.
		label freeNewEdgeI = allEdges.size();

		if (extraVerts.empty())
		{
			// No cuts across this edge. Note that vertices do not need to be
			// renumbered.
			allEdges.append(e);
		}
		else
		{
			// Edge is cut. From e.start() to extraVerts[0],
			// from extraVerts[i] to i+1 and finally to e.end().
			allEdges.append
			(
				edge
				(
					e.start(),
					extraVerts[0] + nSurfacePoints_
				)
			);

			for (label extraI = 1; extraI < extraVerts.size(); extraI++)
			{
				allEdges.append
				(
					edge
					(
					    extraVerts[extraI-1] + nSurfacePoints_,
					    extraVerts[extraI] + nSurfacePoints_
					)
				);
			}
			allEdges.append
			(
				edge
				(
					extraVerts[extraVerts.size()-1] + nSurfacePoints_,
					e.end()
				)
			);
		}

		// Update allFaceEdges, parentEdges_ for the newly added edges.

		// Add each edge label to all face neighbours of edgeI
		const labelList& myFaces = surf.edgeFaces()[edgeI];

		for (label eI = freeNewEdgeI; eI < allEdges.size(); eI++)
		{
			allParentEdges.append(edgeI);

			forAll(myFaces, myFaceI)
			{
				allFaceEdges[myFaces[myFaceI]].append(eI);
			}
		}
	}

	// Done all (possibly split) surface edges by now.
	nSurfaceEdges_ = allEdges.size();


	// Copy intersection edges
	// (note no parentEdges)
	const edgeList& cutEdges = inter.cutEdges();

	forAll(cutEdges, i)
	{
		const edge& e = cutEdges[i];

		allEdges.append(edge(e[0] + nSurfacePoints_, e[1] + nSurfacePoints_));
	}




	// Add intersection edges to faceEdges
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	for
	(
		labelPairLookup::const_iterator iter = inter.facePairToEdge().begin();
		iter != inter.facePairToEdge().end();
		++iter
	)
	{
		// Edge label in intersection
		label edgeI = iter();

		// Get the face from the correct surface
		const FixedList<label, 2>& twoFaces = iter.key();

		label faceI;

		if (isFirstSurface)
		{
			faceI = twoFaces[0];
		}
		else
		{
			faceI = twoFaces[1];
		}

		// Store on face-edge addressing. (note: offset edge)
		allFaceEdges[faceI].append(edgeI + nSurfaceEdges_);
	}

	// Transfer.
	edges_.transfer(allEdges);
	parentEdges_.transfer(allParentEdges);

	forAll(allFaceEdges, faceI)
	{
		faceEdges_[faceI].transfer(allFaceEdges[faceI]);
	}


	// Additional addressing
	// ~~~~~~~~~~~~~~~~~~~~~

	calcPointEdges();


	if (debug & 4)
	{
		Pout<< "edgeSurface : Dumping faceEdges to files" << endl;

		forAll(faceEdges_, faceI)
		{
			const labelList& fEdges = faceEdges_[faceI];

			if (fEdges.size() != 3)
			{
				fileName faceFName("face_" + name(faceI) + ".obj");
				Pout<< "edgeSurface : Dumping faceEdges for face " << faceI
					<< " to " << faceFName << endl;

				OFstream fStream(faceFName);
				writeOBJ(points_, edges_, fEdges, fStream);
			}
		}

		Pout<< "edgeSurface : Dumping edges to edges.obj" << endl;
		OFstream eStream("edges.obj");
		writeOBJ(points_, edges_, eStream);

		Pout<< "edgeSurface : Dumping intersectionEdges to"
			<< " intersectionEdges.obj" << endl;
		OFstream intEdgesStream("intersectionEdges.obj");

		labelList edgeLabels(edges_.size() - nSurfaceEdges_);

		label i = 0;
		for(label edgeI = nSurfaceEdges_; edgeI < edges_.size(); edgeI++)
		{
			edgeLabels[i++] = edgeI;
		}

		writeOBJ(points_, edges_, edgeLabels, intEdgesStream);
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::edgeSurface::addIntersectionEdges
(
	const label faceI,
	const edgeList& additionalEdges
)
{
	if (debug & 2)
	{
		Pout<< "Old face consisted of edges:" << endl;

		const labelList& fEdges = faceEdges_[faceI];
		forAll(fEdges, i)
		{
			const edge& e = edges_[fEdges[i]];

			Pout<< "    " << fEdges[i] << ' ' << e
				<< points_[e.start()] << ' ' << points_[e.end()] << endl;
		}
	}

	// Make space for additional intersection edges (copies old ones)
	const label oldNEdges = edges_.size();

	edges_.setSize(oldNEdges + additionalEdges.size());

	// Append new intersection edges
	label newEdgeI = oldNEdges;

	forAll(additionalEdges, i)
	{
		edges_[newEdgeI] = additionalEdges[i];  // Vertices already in eSurf
					                            // indices.
		newEdgeI++;
	}

	// Append to faceEdges.
	labelList& fEdges = faceEdges_[faceI];

	label nFEdges = fEdges.size();

	fEdges.setSize(nFEdges + additionalEdges.size());

	forAll(additionalEdges, i)
	{
		fEdges[nFEdges++] = oldNEdges + i;
	}


	// Update pointEdge addressing
	calcPointEdges();


	if (debug & 2)
	{
		const labelList& fEdges = faceEdges_[faceI];

		Pout<< "New face consists of edges:" << endl;
		forAll(fEdges, i)
		{
			const edge& e = edges_[fEdges[i]];

			Pout<< "    " << fEdges[i] << ' ' << e
				<< points_[e.start()] << ' ' << points_[e.end()] << endl;
		}
	}
}


// ************************************************************************* //
