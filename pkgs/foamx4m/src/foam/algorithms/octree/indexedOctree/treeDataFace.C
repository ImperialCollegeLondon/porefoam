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

#include "treeDataFace.H"
#include "polyMesh.H"
#include "triangleFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::treeDataFace, 0);

Foam::scalar Foam::treeDataFace::tolSqr = sqr(1E-6);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::treeDataFace::calcBb(const label faceI) const
{
	const pointField& points = mesh_.points();

	const face& f = mesh_.faces()[faceI];

	treeBoundBox bb(points[f[0]], points[f[0]]);

	for (label fp = 1; fp < f.size(); fp++)
	{
		const point& p = points[f[fp]];

		bb.min() = min(bb.min(), p);
		bb.max() = max(bb.max(), p);
	}
	return bb;
}


void Foam::treeDataFace::update()
{
	forAll(faceLabels_, i)
	{
		isTreeFace_.set(faceLabels_[i], 1);
	}

	if (cacheBb_)
	{
		bbs_.setSize(faceLabels_.size());

		forAll(faceLabels_, i)
		{
			bbs_[i] = calcBb(faceLabels_[i]);
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::treeDataFace::treeDataFace
(
	const bool cacheBb,
	const primitiveMesh& mesh,
	const labelList& faceLabels
)
:
	mesh_(mesh),
	faceLabels_(faceLabels),
	isTreeFace_(mesh.nFaces(), 0),
	cacheBb_(cacheBb)
{
	update();
}


Foam::treeDataFace::treeDataFace
(
	const bool cacheBb,
	const primitiveMesh& mesh
)
:
	mesh_(mesh),
	faceLabels_(identity(mesh_.nFaces())),
	isTreeFace_(mesh.nFaces(), 0),
	cacheBb_(cacheBb)
{
	update();
}


Foam::treeDataFace::treeDataFace
(
	const bool cacheBb,
	const polyPatch& patch
)
:
	mesh_(patch.boundaryMesh().mesh()),
	faceLabels_
	(
		identity(patch.size())
	  + patch.start()
	),
	isTreeFace_(mesh_.nFaces(), 0),
	cacheBb_(cacheBb)
{
	update();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataFace::points() const
{
	pointField cc(faceLabels_.size());

	forAll(faceLabels_, i)
	{
		cc[i] = mesh_.faceCentres()[faceLabels_[i]];
	}

	return cc;
}


//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
//  Only makes sense for closed surfaces.
Foam::label Foam::treeDataFace::getVolumeType
(
	const indexedOctree<treeDataFace>& oc,
	const point& sample
) const
{
	// Need to determine whether sample is 'inside' or 'outside'
	// Done by finding nearest face. This gives back a face which is
	// guaranteed to contain nearest point. This point can be
	// - in interior of face: compare to face normal
	// - on edge of face: compare to edge normal
	// - on point of face: compare to point normal
	// Unfortunately the octree does not give us back the intersection point
	// or where on the face it has hit so we have to recreate all that
	// information.


	// Find nearest face to sample
	pointIndexHit info = oc.findNearest(sample, sqr(GREAT));

	if (info.index() == -1)
	{
		FatalErrorIn
		(
			"treeDataFace::getSampleType"
			"(indexedOctree<treeDataFace>&, const point&)"
		)   << "Could not find " << sample << " in octree."
			<< abort(FatalError);
	}


	// Get actual intersection point on face
	label faceI = faceLabels_[info.index()];

	if (debug & 2)
	{
		Pout<< "getSampleType : sample:" << sample
			<< " nearest face:" << faceI;
	}

	const pointField& points = mesh_.points();

	// Retest to classify where on face info is. Note: could be improved. We
	// already have point.

	const face& f = mesh_.faces()[faceI];
	const vector& area = mesh_.faceAreas()[faceI];
	const point& fc = mesh_.faceCentres()[faceI];

	pointHit curHit = f.nearestPoint(sample, points);
	const point& curPt = curHit.rawPoint();

	//
	// 1] Check whether sample is above face
	//

	if (curHit.hit())
	{
		// Nearest point inside face. Compare to face normal.

		if (debug & 2)
		{
			Pout<< " -> face hit:" << curPt
				<< " comparing to face normal " << area << endl;
		}
		return indexedOctree<treeDataFace>::getSide(area, sample - curPt);
	}

	if (debug & 2)
	{
		Pout<< " -> face miss:" << curPt;
	}

	//
	// 2] Check whether intersection is on one of the face vertices or
	//    face centre
	//

	const scalar typDimSqr = mag(area) + VSMALL;

	forAll(f, fp)
	{
		if ((magSqr(points[f[fp]] - curPt)/typDimSqr) < tolSqr)
		{
			// Face intersection point equals face vertex fp

			// Calculate point normal (wrong: uses face normals instead of
			// triangle normals)
			const labelList& pFaces = mesh_.pointFaces()[f[fp]];

			vector pointNormal(vector::zero);

			forAll(pFaces, i)
			{
				if (isTreeFace_.get(pFaces[i]) == 1)
				{
					vector n = mesh_.faceAreas()[pFaces[i]];
					n /= mag(n) + VSMALL;

					pointNormal += n;
				}
			}

			if (debug & 2)
			{
					Pout<< " -> face point hit :" << points[f[fp]]
					    << " point normal:" << pointNormal
					    << " distance:"
					    << magSqr(points[f[fp]] - curPt)/typDimSqr << endl;
			}
			return indexedOctree<treeDataFace>::getSide
			(
				pointNormal,
				sample - curPt
			);
		}
	}
	if ((magSqr(fc - curPt)/typDimSqr) < tolSqr)
	{
		// Face intersection point equals face centre. Normal at face centre
		// is already average of face normals

		if (debug & 2)
		{
			Pout<< " -> centre hit:" << fc
				<< " distance:" << magSqr(fc - curPt)/typDimSqr << endl;
		}

		return indexedOctree<treeDataFace>::getSide(area,  sample - curPt);
	}



	//
	// 3] Get the 'real' edge the face intersection is on
	//

	const labelList& myEdges = mesh_.faceEdges()[faceI];

	forAll(myEdges, myEdgeI)
	{
		const edge& e = mesh_.edges()[myEdges[myEdgeI]];

		pointHit edgeHit =
			line<point, const point&>
			(
				points[e.start()],
				points[e.end()]
			).nearestDist(sample);


		if ((magSqr(edgeHit.rawPoint() - curPt)/typDimSqr) < tolSqr)
		{
			// Face intersection point lies on edge e

			// Calculate edge normal (wrong: uses face normals instead of
			// triangle normals)
			const labelList& eFaces = mesh_.edgeFaces()[myEdges[myEdgeI]];

			vector edgeNormal(vector::zero);

			forAll(eFaces, i)
			{
				if (isTreeFace_.get(eFaces[i]) == 1)
				{
					vector n = mesh_.faceAreas()[eFaces[i]];
					n /= mag(n) + VSMALL;

					edgeNormal += n;
				}
			}

			if (debug & 2)
			{
				Pout<< " -> real edge hit point:" << edgeHit.rawPoint()
					<< " comparing to edge normal:" << edgeNormal
					<< endl;
			}

			// Found face intersection point on this edge. Compare to edge
			// normal
			return indexedOctree<treeDataFace>::getSide
			(
				edgeNormal,
				sample - curPt
			);
		}
	}


	//
	// 4] Get the internal edge the face intersection is on
	//

	forAll(f, fp)
	{
		pointHit edgeHit = line<point, const point&>
		(
			points[f[fp]],
			fc
		).nearestDist(sample);

		if ((magSqr(edgeHit.rawPoint() - curPt)/typDimSqr) < tolSqr)
		{
			// Face intersection point lies on edge between two face triangles

			// Calculate edge normal as average of the two triangle normals
			vector e = points[f[fp]] - fc;
			vector ePrev = points[f[f.rcIndex(fp)]] - fc;
			vector eNext = points[f[f.fcIndex(fp)]] - fc;

			vector nLeft = ePrev ^ e;
			nLeft /= mag(nLeft) + VSMALL;

			vector nRight = e ^ eNext;
			nRight /= mag(nRight) + VSMALL;

			if (debug & 2)
			{
				Pout<< " -> internal edge hit point:" << edgeHit.rawPoint()
					<< " comparing to edge normal "
					<< 0.5*(nLeft + nRight)
					<< endl;
			}

			// Found face intersection point on this edge. Compare to edge
			// normal
			return indexedOctree<treeDataFace>::getSide
			(
				0.5*(nLeft + nRight),
				sample - curPt
			);
		}
	}

	if (debug & 2)
	{
		Pout<< "Did not find sample " << sample
			<< " anywhere related to nearest face " << faceI << endl
			<< "Face:";

		forAll(f, fp)
		{
			Pout<< "    vertex:" << f[fp] << "  coord:" << points[f[fp]]
				<< endl;
		}
	}

	// Can't determine status of sample with respect to nearest face.
	// Either
	// - tolerances are wrong. (if e.g. face has zero area)
	// - or (more likely) surface is not closed.

	return indexedOctree<treeDataFace>::UNKNOWN;
}


// Check if any point on shape is inside cubeBb.
bool Foam::treeDataFace::overlaps
(
	const label index,
	const treeBoundBox& cubeBb
) const
{
	// 1. Quick rejection: bb does not intersect face bb at all
	if (cacheBb_)
	{
		if (!cubeBb.overlaps(bbs_[index]))
		{
			return false;
		}
	}
	else
	{
		if (!cubeBb.overlaps(calcBb(faceLabels_[index])))
		{
			return false;
		}
	}

	const pointField& points = mesh_.points();


	// 2. Check if one or more face points inside
	label faceI = faceLabels_[index];

	const face& f = mesh_.faces()[faceI];

	forAll(f, fp)
	{
		if (cubeBb.contains(points[f[fp]]))
		{
			return true;
		}
	}

	// 3. Difficult case: all points are outside but connecting edges might
	// go through cube. Use triangle-bounding box intersection.
	const point& fc = mesh_.faceCentres()[faceI];

	forAll(f, fp)
	{
		bool triIntersects = triangleFuncs::intersectBb
		(
			points[f[fp]],
			points[f[f.fcIndex(fp)]],
			fc,
			cubeBb
		);

		if (triIntersects)
		{
			return true;
		}
	}
	return false;
}


// Calculate nearest point to sample. Updates (if any) nearestDistSqr, minIndex,
// nearestPoint.
void Foam::treeDataFace::findNearest
(
	const labelList& indices,
	const point& sample,

	scalar& nearestDistSqr,
	label& minIndex,
	point& nearestPoint
) const
{
	forAll(indices, i)
	{
		label index = indices[i];

		const face& f = mesh_.faces()[faceLabels_[index]];

		pointHit nearHit = f.nearestPoint(sample, mesh_.points());
		scalar distSqr = sqr(nearHit.distance());

		if (distSqr < nearestDistSqr)
		{
			nearestDistSqr = distSqr;
			minIndex = index;
			nearestPoint = nearHit.rawPoint();
		}
	}
}


bool Foam::treeDataFace::intersects
(
	const label index,
	const point& start,
	const point& end,
	point& intersectionPoint
) const
{
	// Do quick rejection test
	if (cacheBb_)
	{
		const treeBoundBox& faceBb = bbs_[index];

		if ((faceBb.posBits(start) & faceBb.posBits(end)) != 0)
		{
			// start and end in same block outside of faceBb.
			return false;
		}
	}

	label faceI = faceLabels_[index];

	const vector dir(end - start);

	pointHit inter = mesh_.faces()[faceI].fastIntersection
	(
		start,
		dir,
		mesh_.faceCentres()[faceI],
		mesh_.points(),
		intersection::HALF_RAY
	);

	if (inter.hit() && inter.distance() <= 1)
	{
		// Note: no extra test on whether intersection is in front of us
		// since using half_ray
		intersectionPoint = inter.hitPoint();
		return true;
	}
	else
	{
		return false;
	}
}


// ************************************************************************* //
