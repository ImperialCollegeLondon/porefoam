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

#include "treeDataPoint.H"
#include "treeBoundBox.H"
#include "indexedOctree.H"
#include "polyMesh.H"
#include "triangleFuncs.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::treeDataPoint, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::treeDataPoint::treeDataPoint(const pointField& points)
:
	points_(points)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataPoint::points() const
{
	return points_;
}


//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
//  Only makes sense for closed surfaces.
Foam::label Foam::treeDataPoint::getVolumeType
(
	const indexedOctree<treeDataPoint>& oc,
	const point& sample
) const
{
	return indexedOctree<treeDataPoint>::UNKNOWN;
}


// Check if any point on shape is inside cubeBb.
bool Foam::treeDataPoint::overlaps
(
	const label index,
	const treeBoundBox& cubeBb
) const
{
	return cubeBb.contains(points_[index]);
}


// Calculate nearest point to sample. Updates (if any) nearestDistSqr,
//  minIndex, nearestPoint.
void Foam::treeDataPoint::findNearest
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

		const point& pt = points_[index];

		scalar distSqr = magSqr(pt - sample);

		if (distSqr < nearestDistSqr)
		{
			nearestDistSqr = distSqr;
			minIndex = index;
			nearestPoint = pt;
		}
	}
}


//- Calculates nearest (to line) point in shape.
//  Returns point and distance (squared)
void Foam::treeDataPoint::findNearest
(
	const labelList& indices,
	const linePointRef& ln,

	treeBoundBox& tightest,
	label& minIndex,
	point& linePoint,
	point& nearestPoint
) const
{
	// Best so far
	scalar nearestDistSqr = magSqr(linePoint - nearestPoint);

	forAll(indices, i)
	{
		label index = indices[i];

		const point& shapePt = points_[index];

		if (tightest.contains(shapePt))
		{
			// Nearest point on line
			pointHit pHit = ln.nearestDist(shapePt);
			scalar distSqr = sqr(pHit.distance());

			if (distSqr < nearestDistSqr)
			{
				nearestDistSqr = distSqr;
				minIndex = index;
				linePoint = pHit.rawPoint();
				nearestPoint = shapePt;

				{
					point& minPt = tightest.min();
					minPt = min(ln.start(), ln.end());
					minPt.x() -= pHit.distance();
					minPt.y() -= pHit.distance();
					minPt.z() -= pHit.distance();
				}
				{
					point& maxPt = tightest.max();
					maxPt = max(ln.start(), ln.end());
					maxPt.x() += pHit.distance();
					maxPt.y() += pHit.distance();
					maxPt.z() += pHit.distance();
				}
			}
		}
	}
}


// ************************************************************************* //
