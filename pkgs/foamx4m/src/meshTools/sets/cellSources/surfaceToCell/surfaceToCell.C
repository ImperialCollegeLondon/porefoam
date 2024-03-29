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

#include "surfaceToCell.H"
#include "polyMesh.H"
#include "meshSearch.H"
#include "triSurface.H"
#include "triSurfaceSearch.H"
#include "cellClassification.H"

#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(surfaceToCell, 0);

addToRunTimeSelectionTable(topoSetSource, surfaceToCell, word);

addToRunTimeSelectionTable(topoSetSource, surfaceToCell, istream);

}


Foam::topoSetSource::addToUsageTable Foam::surfaceToCell::usage_
(
	surfaceToCell::typeName,
	"\n    Usage: surfaceToCell"
	"<surface> <outsidePoints> <cut> <inside> <outside> <near> <curvature>\n\n"
	"    <surface> name of triSurface\n"
	"    <outsidePoints> list of points that define outside\n"
	"    <cut> boolean whether to include cells cut by surface\n"
	"    <inside>   ,,                 ,,       inside surface\n"
	"    <outside>  ,,                 ,,       outside surface\n"
	"    <near> scalar; include cells with centre <= near to surface\n"
	"    <curvature> scalar; include cells close to strong curvature"
	" on surface\n"
	"    (curvature defined as difference in surface normal at nearest"
	" point on surface for each vertex of cell)\n\n"
);


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::label Foam::surfaceToCell::getNearest
(
	const triSurfaceSearch& querySurf,
	const label pointI,
	const point& pt,
	const vector& span,
	Map<label>& cache
)
{
	Map<label>::const_iterator iter = cache.find(pointI);

	if (iter != cache.end())
	{
		// Found cached answer
		return iter();
	}
	else
	{
		pointIndexHit inter = querySurf.nearest(pt, span);

		// Triangle label (can be -1)
		label triI = inter.index();

		// Store triangle on point
		cache.insert(pointI, triI);

		return triI;
	}
}


// Return true if nearest surface to points on cell makes largish angle
// with nearest surface to cell centre. Returns false otherwise. Points visited
// are cached in pointToNearest
bool Foam::surfaceToCell::differingPointNormals
(
	const triSurfaceSearch& querySurf,

	const vector& span,         // current search span
	const label cellI,
	const label cellTriI,       // nearest (to cell centre) surface triangle

	Map<label>& pointToNearest  // cache for nearest triangle to point
) const
{
	const triSurface& surf = querySurf.surface();
	const vectorField& normals = surf.faceNormals();

	const faceList& faces = mesh().faces();
	const pointField& points = mesh().points();

	const labelList& cFaces = mesh().cells()[cellI];

	forAll(cFaces, cFaceI)
	{
		const face& f = faces[cFaces[cFaceI]];

		forAll(f, fp)
		{
			label pointI = f[fp];

			label pointTriI =
				getNearest
				(
					querySurf,
					pointI,
					points[pointI],
					span,
					pointToNearest
				);

			if (pointTriI != -1 && pointTriI != cellTriI)
			{
				scalar cosAngle = normals[pointTriI] & normals[cellTriI];

				if (cosAngle < 0.9)
				{
					return true;
				}
			}
		}
	}
	return false;
}


void Foam::surfaceToCell::combine(topoSet& set, const bool add) const
{
	cpuTime timer;

	if (includeCut_ || includeInside_ || includeOutside_)
	{
		//
		// Cut cells with surface and classify cells
		//


		// Construct search engine on mesh

		meshSearch queryMesh(mesh_, true);


		// Check all 'outside' points
		forAll(outsidePoints_, outsideI)
		{
			const point& outsidePoint = outsidePoints_[outsideI];

			// Find cell point is in. Linear search.
			if (queryMesh.findCell(outsidePoint, -1, false) == -1)
			{
				FatalErrorIn("surfaceToCell::combine(topoSet&, const bool)")
					<< "outsidePoint " << outsidePoint
					<< " is not inside any cell"
					<< exit(FatalError);
			}
		}

		// Cut faces with surface and classify cells

		cellClassification cellType
		(
			mesh_,
			queryMesh,
			querySurf(),
			outsidePoints_
		);


		Info<< "    Marked inside/outside in = "
			<< timer.cpuTimeIncrement() << " s" << endl << endl;


		forAll(cellType, cellI)
		{
			label cType = cellType[cellI];

			if
			(
				(
					includeCut_
				 && (cType == cellClassification::CUT)
				)
			 || (
					includeInside_
				 && (cType == cellClassification::INSIDE)
				)
			 || (
					includeOutside_
				 && (cType == cellClassification::OUTSIDE)
				)
			)
			{
				addOrDelete(set, cellI, add);
			}
		}
	}


	if (nearDist_ > 0)
	{
		//
		// Determine distance to surface
		//

		const pointField& ctrs = mesh_.cellCentres();

		// Box dimensions to search in octree.
		const vector span(nearDist_, nearDist_, nearDist_);


		if (curvature_ < -1)
		{
			Info<< "    Selecting cells with cellCentre closer than "
				<< nearDist_ << " to surface" << endl;

			// No need to test curvature. Insert near cells into set.

			forAll(ctrs, cellI)
			{
				const point& c = ctrs[cellI];

				pointIndexHit inter = querySurf().nearest(c, span);

				if (inter.hit() && (mag(inter.hitPoint() - c) < nearDist_))
				{
					addOrDelete(set, cellI, add);
				}
			}

			Info<< "    Determined nearest surface point in = "
				<< timer.cpuTimeIncrement() << " s" << endl << endl;

		}
		else
		{
			// Test near cells for curvature

			Info<< "    Selecting cells with cellCentre closer than "
				<< nearDist_ << " to surface and curvature factor"
				<< " less than " << curvature_ << endl;

			// Cache for nearest surface triangle for a point
			Map<label> pointToNearest(mesh_.nCells()/10);

			forAll(ctrs, cellI)
			{
				const point& c = ctrs[cellI];

				pointIndexHit inter = querySurf().nearest(c, span);

				if (inter.hit() && (mag(inter.hitPoint() - c) < nearDist_))
				{
					if
					(
					    differingPointNormals
					    (
					        querySurf(),
					        span,
					        cellI,
					        inter.index(),      // nearest surface triangle
					        pointToNearest
					    )
					)
					{
					    addOrDelete(set, cellI, add);
					}
				}
			}

			Info<< "    Determined nearest surface point in = "
				<< timer.cpuTimeIncrement() << " s" << endl << endl;
		}
	}
}


void Foam::surfaceToCell::checkSettings() const
{
	if
	(
		(nearDist_ < 0)
	 && (curvature_ < -1)
	 && (
			(includeCut_ && includeInside_ && includeOutside_)
		 || (!includeCut_ && !includeInside_ && !includeOutside_)
		)
	)
	{
		FatalErrorIn
		(
			"surfaceToCell:checkSettings()"
		)   << "Illegal include cell specification."
			<< " Result would be either all or no cells." << endl
			<< "Please set one of includeCut, includeInside, includeOutside"
			<< " to true, set nearDistance to a value > 0"
			<< " or set curvature to a value -1 .. 1."
			<< exit(FatalError);
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::surfaceToCell::surfaceToCell
(
	const polyMesh& mesh,
	const fileName& surfName,
	const pointField& outsidePoints,
	const bool includeCut,
	const bool includeInside,
	const bool includeOutside,
	const scalar nearDist,
	const scalar curvature
)
:
	topoSetSource(mesh),
	surfName_(surfName),
	outsidePoints_(outsidePoints),
	includeCut_(includeCut),
	includeInside_(includeInside),
	includeOutside_(includeOutside),
	nearDist_(nearDist),
	curvature_(curvature),
	surfPtr_(new triSurface(surfName_)),
	querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
	IOwnPtrs_(true)
{
	checkSettings();
}


// Construct from components. Externally supplied surface.
Foam::surfaceToCell::surfaceToCell
(
	const polyMesh& mesh,
	const fileName& surfName,
	const triSurface& surf,
	const triSurfaceSearch& querySurf,
	const pointField& outsidePoints,
	const bool includeCut,
	const bool includeInside,
	const bool includeOutside,
	const scalar nearDist,
	const scalar curvature
)
:
	topoSetSource(mesh),
	surfName_(surfName),
	outsidePoints_(outsidePoints),
	includeCut_(includeCut),
	includeInside_(includeInside),
	includeOutside_(includeOutside),
	nearDist_(nearDist),
	curvature_(curvature),
	surfPtr_(&surf),
	querySurfPtr_(&querySurf),
	IOwnPtrs_(false)
{
	checkSettings();
}


// Construct from dictionary
Foam::surfaceToCell::surfaceToCell
(
	const polyMesh& mesh,
	const dictionary& dict
)
:
	topoSetSource(mesh),
	surfName_(dict.lookup("file")),
	outsidePoints_(dict.lookup("outsidePoints")),
	includeCut_(readBool(dict.lookup("includeCut"))),
	includeInside_(readBool(dict.lookup("includeInside"))),
	includeOutside_(readBool(dict.lookup("includeOutside"))),
	nearDist_(readScalar(dict.lookup("nearDistance"))),
	curvature_(readScalar(dict.lookup("curvature"))),
	surfPtr_(new triSurface(surfName_)),
	querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
	IOwnPtrs_(true)
{
	checkSettings();
}


// Construct from Istream
Foam::surfaceToCell::surfaceToCell
(
	const polyMesh& mesh,
	Istream& is
)
:
	topoSetSource(mesh),
	surfName_(checkIs(is)),
	outsidePoints_(checkIs(is)),
	includeCut_(readBool(checkIs(is))),
	includeInside_(readBool(checkIs(is))),
	includeOutside_(readBool(checkIs(is))),
	nearDist_(readScalar(checkIs(is))),
	curvature_(readScalar(checkIs(is))),
	surfPtr_(new triSurface(surfName_)),
	querySurfPtr_(new triSurfaceSearch(*surfPtr_)),
	IOwnPtrs_(true)
{
	checkSettings();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfaceToCell::~surfaceToCell()
{
	if (IOwnPtrs_)
	{
		if (surfPtr_)
		{
			delete surfPtr_;
		}
		if (querySurfPtr_)
		{
			delete querySurfPtr_;
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::surfaceToCell::applyToSet
(
	const topoSetSource::setAction action,
	topoSet& set
) const
{
	if ( (action == topoSetSource::NEW) || (action == topoSetSource::ADD))
	{
		Info<< "    Adding cells in relation to surface " << surfName_
			<< " ..." << endl;

		combine(set, true);
	}
	else if (action == topoSetSource::DELETE)
	{
		Info<< "    Removing cells in relation to surface " << surfName_
			<< " ..." << endl;

		combine(set, false);
	}
}


// ************************************************************************* //
