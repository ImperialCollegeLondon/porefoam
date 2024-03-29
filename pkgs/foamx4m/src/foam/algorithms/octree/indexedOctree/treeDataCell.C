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

#include "treeDataCell.H"
#include "indexedOctree.H"
#include "primitiveMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::treeDataCell, 0);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::treeBoundBox Foam::treeDataCell::calcCellBb(const label cellI) const
{
	const cellList& cells = mesh_.cells();
	const faceList& faces = mesh_.faces();
	const pointField& points = mesh_.points();

	treeBoundBox cellBb
	(
		vector(GREAT, GREAT, GREAT),
		vector(-GREAT, -GREAT, -GREAT)
	);

	const cell& cFaces = cells[cellI];

	forAll(cFaces, cFaceI)
	{
		const face& f = faces[cFaces[cFaceI]];

		forAll(f, fp)
		{
			const point& p = points[f[fp]];

			cellBb.min() = min(cellBb.min(), p);
			cellBb.max() = max(cellBb.max(), p);
		}
	}
	return cellBb;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::treeDataCell::treeDataCell
(
	const bool cacheBb,
	const primitiveMesh& mesh,
	const labelList& cellLabels
)
:
	mesh_(mesh),
	cellLabels_(cellLabels),
	cacheBb_(cacheBb)
{
	if (cacheBb_)
	{
		bbs_.setSize(cellLabels_.size());

		forAll(cellLabels_, i)
		{
			bbs_[i] = calcCellBb(cellLabels_[i]);
		}
	}
}


Foam::treeDataCell::treeDataCell
(
	const bool cacheBb,
	const primitiveMesh& mesh
)
:
	mesh_(mesh),
	cellLabels_(identity(mesh_.nCells())),
	cacheBb_(cacheBb)
{
	if (cacheBb_)
	{
		bbs_.setSize(cellLabels_.size());

		forAll(cellLabels_, i)
		{
			bbs_[i] = calcCellBb(cellLabels_[i]);
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::treeDataCell::points() const
{
	pointField cc(cellLabels_.size());

	forAll(cellLabels_, i)
	{
		cc[i] = mesh_.cellCentres()[cellLabels_[i]];
	}

	return cc;
}


// Check if any point on shape is inside cubeBb.
bool Foam::treeDataCell::overlaps
(
	const label index,
	const treeBoundBox& cubeBb
) const
{
	if (cacheBb_)
	{
		return cubeBb.overlaps(bbs_[index]);
	}
	else
	{
		return cubeBb.overlaps(calcCellBb(cellLabels_[index]));
	}
}


// Calculate nearest point to sample. Updates (if any) nearestDistSqr, minIndex,
// nearestPoint.
void Foam::treeDataCell::findNearest
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
		label cellI = cellLabels_[index];
		scalar distSqr = magSqr(sample - mesh_.cellCentres()[cellI]);

		if (distSqr < nearestDistSqr)
		{
			nearestDistSqr = distSqr;
			minIndex = index;
			nearestPoint = mesh_.cellCentres()[cellI];
		}
	}
}


bool Foam::treeDataCell::intersects
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
		const treeBoundBox& cellBb = bbs_[index];

		if ((cellBb.posBits(start) & cellBb.posBits(end)) != 0)
		{
			// start and end in same block outside of cellBb.
			return false;
		}
	}
	else
	{
		const treeBoundBox cellBb = calcCellBb(cellLabels_[index]);

		if ((cellBb.posBits(start) & cellBb.posBits(end)) != 0)
		{
			// start and end in same block outside of cellBb.
			return false;
		}
	}


	// Do intersection with all faces of cell
	// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	// Disable picking up intersections behind us.
	scalar oldTol = intersection::setPlanarTol(0.0);

	const cell& cFaces = mesh_.cells()[cellLabels_[index]];

	const vector dir(end - start);
	scalar minDistSqr = magSqr(dir);
	bool hasMin = false;

	forAll(cFaces, i)
	{
		const face& f = mesh_.faces()[cFaces[i]];

		pointHit inter = f.ray
		(
			start,
			dir,
			mesh_.points(),
			intersection::HALF_RAY
		);

		if (inter.hit() && sqr(inter.distance()) <= minDistSqr)
		{
			// Note: no extra test on whether intersection is in front of us
			// since using half_ray AND zero tolerance. (note that tolerance
			// is used to look behind us)
			minDistSqr = sqr(inter.distance());
			intersectionPoint = inter.hitPoint();
			hasMin = true;
		}
	}

	// Restore picking tolerance
	intersection::setPlanarTol(oldTol);

	return hasMin;
}


// ************************************************************************* //
