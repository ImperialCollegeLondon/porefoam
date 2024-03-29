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

#include "triSurfaceMesh.H"
#include "Random.H"
#include "addToRunTimeSelectionTable.H"
#include "EdgeMap.H"
#include "triSurfaceFields.H"
#include "foamTime.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{

defineTypeNameAndDebug(triSurfaceMesh, 0);
addToRunTimeSelectionTable(searchableSurface, triSurfaceMesh, dict);

}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Check file existence
const Foam::fileName& Foam::triSurfaceMesh::checkFile
(
	const fileName& fName,
	const fileName& objectName
)
{
	if (fName.empty())
	{
		FatalErrorIn
		(
			"triSurfaceMesh::checkFile(const fileName&, const fileName&)"
		)   << "Cannot find triSurfaceMesh starting from "
			<< objectName << exit(FatalError);
	}

	return fName;
}


bool Foam::triSurfaceMesh::addFaceToEdge
(
	const edge& e,
	EdgeMap<label>& facesPerEdge
)
{
	EdgeMap<label>::iterator eFnd = facesPerEdge.find(e);
	if (eFnd != facesPerEdge.end())
	{
		if (eFnd() == 2)
		{
			return false;
		}
		eFnd()++;
	}
	else
	{
		facesPerEdge.insert(e, 1);
	}
	return true;
}


bool Foam::triSurfaceMesh::isSurfaceClosed() const
{
	// Construct pointFaces. Let's hope surface has compact point
	// numbering ...
	labelListList pointFaces;
	invertManyToMany(points().size(), *this, pointFaces);

	// Loop over all faces surrounding point. Count edges emanating from point.
	// Every edge should be used by two faces exactly.
	// To prevent doing work twice per edge only look at edges to higher
	// point
	EdgeMap<label> facesPerEdge(100);
	forAll(pointFaces, pointI)
	{
		const labelList& pFaces = pointFaces[pointI];

		facesPerEdge.clear();
		forAll(pFaces, i)
		{
			const labelledTri& f = triSurface::operator[](pFaces[i]);
			label fp = findIndex(f, pointI);

			// Something weird: if I expand the code of addFaceToEdge in both
			// below instances it gives a segmentation violation on some
			// surfaces. Compiler (4.3.2) problem?


			// Forward edge
			label nextPointI = f[f.fcIndex(fp)];

			if (nextPointI > pointI)
			{
				bool okFace = addFaceToEdge
				(
					edge(pointI, nextPointI),
					facesPerEdge
				);

				if (!okFace)
				{
					return false;
				}
			}
			// Reverse edge
			label prevPointI = f[f.rcIndex(fp)];

			if (prevPointI > pointI)
			{
				bool okFace = addFaceToEdge
				(
					edge(pointI, prevPointI),
					facesPerEdge
				);

				if (!okFace)
				{
					return false;
				}
			}
		}

		// Check for any edges used only once.
		forAllConstIter(EdgeMap<label>, facesPerEdge, iter)
		{
			if (iter() != 2)
			{
				return false;
			}
		}
	}

	return true;
}


// Gets all intersections after initial one. Adds smallVec and starts tracking
// from there.
void Foam::triSurfaceMesh::getNextIntersections
(
	const indexedOctree<treeDataTriSurface>& octree,
	const point& start,
	const point& end,
	const vector& smallVec,
	DynamicList<pointIndexHit, 1, 1>& hits
)
{
	const vector dirVec(end-start);
	const scalar magSqrDirVec(magSqr(dirVec));

	// Initial perturbation amount
	vector perturbVec(smallVec);

	while (true)
	{
		// Start tracking from last hit.
		point pt = hits[hits.size()-1].hitPoint() + perturbVec;

		if (((pt-start)&dirVec) > magSqrDirVec)
		{
			return;
		}

		// See if any intersection between pt and end
		pointIndexHit inter = octree.findLine(pt, end);

		if (!inter.hit())
		{
			return;
		}

		// Check if already found this intersection
		bool duplicateHit = false;
		forAllReverse(hits, i)
		{
			if (hits[i].index() == inter.index())
			{
				duplicateHit = true;
				break;
			}
		}


		if (duplicateHit)
		{
			// Hit same triangle again. Increase perturbVec and try again.
			perturbVec *= 2;
		}
		else
		{
			// Proper hit
			hits.append(inter);
			// Restore perturbVec
			perturbVec = smallVec;
		}
	}
}


void Foam::triSurfaceMesh::calcBounds(boundBox& bb, label& nPoints) const
{
	// Unfortunately nPoints constructs meshPoints() so do compact version
	// ourselves

	const triSurface& s = static_cast<const triSurface&>(*this);

	PackedBoolList pointIsUsed(points().size());

	nPoints = 0;
	bb = boundBox::invertedBox;

	forAll(s, triI)
	{
		const labelledTri& f = s[triI];

		forAll(f, fp)
		{
			label pointI = f[fp];
			if (pointIsUsed.set(pointI, 1u))
			{
				bb.min() = ::Foam::min(bb.min(), points()[pointI]);
				bb.max() = ::Foam::max(bb.max(), points()[pointI]);
				nPoints++;
			}
		}
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io, const triSurface& s)
:
	searchableSurface(io),
	objectRegistry
	(
		IOobject
		(
			io.name(),
			io.instance(),
			io.local(),
			io.db(),
			io.readOpt(),
			io.writeOpt(),
			false       // searchableSurface already registered under name
		)
	),
	triSurface(s),
	tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
	maxTreeDepth_(10),
	surfaceClosed_(-1)
{}


Foam::triSurfaceMesh::triSurfaceMesh(const IOobject& io)
:
	// Find instance for triSurfaceMesh
	searchableSurface(io),
	objectRegistry
	(
		IOobject
		(
			io.name(),
			static_cast<const searchableSurface&>(*this).instance(),
			io.local(),
			io.db(),
			io.readOpt(),
			io.writeOpt(),
			false       // searchableSurface already registered under name
		)
	),
	triSurface
	(
		checkFile
		(
			searchableSurface::filePath(),
			searchableSurface::objectPath()
		)
	),
	tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
	maxTreeDepth_(10),
	surfaceClosed_(-1)
{}


Foam::triSurfaceMesh::triSurfaceMesh
(
	const IOobject& io,
	const dictionary& dict
)
:
	searchableSurface(io),
	objectRegistry
	(
		IOobject
		(
			io.name(),
			static_cast<const searchableSurface&>(*this).instance(),
			io.local(),
			io.db(),
			io.readOpt(),
			io.writeOpt(),
			false       // searchableSurface already registered under name
		)
	),
	triSurface
	(
		checkFile
		(
			searchableSurface::filePath(),
			searchableSurface::objectPath()
		)
	),
	tolerance_(indexedOctree<treeDataTriSurface>::perturbTol()),
	maxTreeDepth_(10),
	surfaceClosed_(-1)
{
	scalar scaleFactor = 0;

	// Allow rescaling of the surface points
	// eg, CAD geometries are often done in millimeters
	if (dict.readIfPresent("scale", scaleFactor) && scaleFactor > 0)
	{
		Info<< searchableSurface::name() << " : using scale " << scaleFactor
			<< endl;
		triSurface::scalePoints(scaleFactor);
	}


	// Have optional non-standard search tolerance for gappy surfaces.
	if (dict.readIfPresent("tolerance", tolerance_) && tolerance_ > 0)
	{
		Info<< searchableSurface::name() << " : using intersection tolerance "
			<< tolerance_ << endl;
	}


	// Have optional non-standard tree-depth to limit storage.
	if (dict.readIfPresent("maxTreeDepth", maxTreeDepth_) && maxTreeDepth_ > 0)
	{
		Info<< searchableSurface::name() << " : using maximum tree depth "
			<< maxTreeDepth_ << endl;
	}
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::triSurfaceMesh::~triSurfaceMesh()
{
	clearOut();
}


void Foam::triSurfaceMesh::clearOut()
{
	tree_.clear();
	edgeTree_.clear();
	triSurface::clearOut();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::pointField Foam::triSurfaceMesh::coordinates() const
{
	// Use copy to calculate face centres so they don't get stored
	return PrimitivePatch<labelledTri, SubList, const pointField&>
	(
		SubList<labelledTri>(*this, triSurface::size()),
		triSurface::points()
	).faceCentres();
}


void Foam::triSurfaceMesh::movePoints(const pointField& newPoints)
{
	tree_.clear();
	edgeTree_.clear();
	triSurface::movePoints(newPoints);
}


const Foam::indexedOctree<Foam::treeDataTriSurface>&
Foam::triSurfaceMesh::tree() const
{
	if (tree_.empty())
	{
		// Calculate bb without constructing local point numbering.
		treeBoundBox bb;
		label nPoints;
		calcBounds(bb, nPoints);

		if (nPoints != points().size())
		{
			WarningIn("triSurfaceMesh::tree() const")
				<< "Surface " << searchableSurface::name()
				<< " does not have compact point numbering."
				<< " Of " << points().size() << " only " << nPoints
				<< " are used. This might give problems in some routines."
				<< endl;
		}


		// Random number generator. Bit dodgy since not exactly random ;-)
		Random rndGen(65431);

		// Slightly extended bb. Slightly off-centred just so on symmetric
		// geometry there are less face/edge aligned items.
		bb = bb.extend(rndGen, 1E-4);
		bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
		bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

		scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
		indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

		tree_.reset
		(
			new indexedOctree<treeDataTriSurface>
			(
				treeDataTriSurface(*this),
				bb,
				maxTreeDepth_,  // maxLevel
				10,             // leafsize
				3.0             // duplicity
			)
		);

		indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
	}

	return tree_();
}


const Foam::indexedOctree<Foam::treeDataEdge>&
 Foam::triSurfaceMesh::edgeTree() const
{
	if (edgeTree_.empty())
	{
		// Boundary edges
		labelList bEdges
		(
			identity
			(
				nEdges()
			   -nInternalEdges()
			)
		  + nInternalEdges()
		);

		treeBoundBox bb;
		label nPoints;
		calcBounds(bb, nPoints);

		// Random number generator. Bit dodgy since not exactly random ;-)
		Random rndGen(65431);

		// Slightly extended bb. Slightly off-centred just so on symmetric
		// geometry there are less face/edge aligned items.

		bb = bb.extend(rndGen, 1E-4);
		bb.min() -= point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);
		bb.max() += point(ROOTVSMALL, ROOTVSMALL, ROOTVSMALL);

		edgeTree_.reset
		(
			new indexedOctree<treeDataEdge>
			(
				treeDataEdge
				(
					false,          // cachebb
					edges(),        // edges
					localPoints(),  // points
					bEdges          // selected edges
				),
				bb,                 // bb
				maxTreeDepth_,      // maxLevel
				10,                 // leafsize
				3.0                 // duplicity
			)
		);
	}
	return edgeTree_();
}


const Foam::wordList& Foam::triSurfaceMesh::regions() const
{
	if (regions_.empty())
	{
		regions_.setSize(patches().size());
		forAll(regions_, regionI)
		{
			regions_[regionI] = patches()[regionI].name();
		}
	}
	return regions_;
}


// Find out if surface is closed.
bool Foam::triSurfaceMesh::hasVolumeType() const
{
	if (surfaceClosed_ == -1)
	{
		if (isSurfaceClosed())
		{
			surfaceClosed_ = 1;
		}
		else
		{
			surfaceClosed_ = 0;
		}
	}

	return surfaceClosed_ == 1;
}


void Foam::triSurfaceMesh::findNearest
(
	const pointField& samples,
	const scalarField& nearestDistSqr,
	List<pointIndexHit>& info
) const
{
	const indexedOctree<treeDataTriSurface>& octree = tree();

	info.setSize(samples.size());

	scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
	indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

	forAll(samples, i)
	{
		static_cast<pointIndexHit&>(info[i]) = octree.findNearest
		(
			samples[i],
			nearestDistSqr[i]
		);
	}

	indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::findLine
(
	const pointField& start,
	const pointField& end,
	List<pointIndexHit>& info
) const
{
	const indexedOctree<treeDataTriSurface>& octree = tree();

	info.setSize(start.size());

	scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
	indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

	forAll(start, i)
	{
		static_cast<pointIndexHit&>(info[i]) = octree.findLine
		(
			start[i],
			end[i]
		);
	}

	indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::findLineAny
(
	const pointField& start,
	const pointField& end,
	List<pointIndexHit>& info
) const
{
	const indexedOctree<treeDataTriSurface>& octree = tree();

	info.setSize(start.size());

	scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
	indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

	forAll(start, i)
	{
		static_cast<pointIndexHit&>(info[i]) = octree.findLineAny
		(
			start[i],
			end[i]
		);
	}

	indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::findLineAll
(
	const pointField& start,
	const pointField& end,
	List<List<pointIndexHit> >& info
) const
{
	const indexedOctree<treeDataTriSurface>& octree = tree();

	info.setSize(start.size());

	scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
	indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

	// Work array
	DynamicList<pointIndexHit, 1, 1> hits;

	// Tolerances:
	// To find all intersections we add a small vector to the last intersection
	// This is chosen such that
	// - it is significant (SMALL is smallest representative relative tolerance;
	//   we need something bigger since we're doing calculations)
	// - if the start-end vector is zero we still progress
	const vectorField dirVec(end-start);
	const scalarField magSqrDirVec(magSqr(dirVec));
	const vectorField smallVec
	(
		indexedOctree<treeDataTriSurface>::perturbTol()*dirVec
	  + vector(ROOTVSMALL,ROOTVSMALL,ROOTVSMALL)
	);

	forAll(start, pointI)
	{
		// See if any intersection between pt and end
		pointIndexHit inter = octree.findLine(start[pointI], end[pointI]);

		if (inter.hit())
		{
			hits.clear();
			hits.append(inter);

			getNextIntersections
			(
				octree,
				start[pointI],
				end[pointI],
				smallVec[pointI],
				hits
			);

			info[pointI].transfer(hits);
		}
		else
		{
			info[pointI].clear();
		}
	}

	indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


void Foam::triSurfaceMesh::getRegion
(
	const List<pointIndexHit>& info,
	labelList& region
) const
{
	region.setSize(info.size());
	forAll(info, i)
	{
		if (info[i].hit())
		{
			region[i] = triSurface::operator[](info[i].index()).region();
		}
		else
		{
			region[i] = -1;
		}
	}
}


void Foam::triSurfaceMesh::getNormal
(
	const List<pointIndexHit>& info,
	vectorField& normal
) const
{
	normal.setSize(info.size());

	forAll(info, i)
	{
		if (info[i].hit())
		{
			label triI = info[i].index();
			//- Cached:
			//normal[i] = faceNormals()[triI];

			//- Uncached
			normal[i] = triSurface::operator[](triI).normal(points());
			normal[i] /= mag(normal[i]) + VSMALL;
		}
		else
		{
			// Set to what?
			normal[i] = vector::zero;
		}
	}
}


void Foam::triSurfaceMesh::setField(const labelList& values)
{
	autoPtr<triSurfaceLabelField> fldPtr
	(
		new triSurfaceLabelField
		(
			IOobject
			(
				"values",
				objectRegistry::time().timeName(),  // instance
				"triSurface",                       // local
				*this,
				IOobject::NO_READ,
				IOobject::AUTO_WRITE
			),
			*this,
			dimless,
			labelField(values)
		)
	);

	// Store field on triMesh
	fldPtr.ptr()->store();
}


void Foam::triSurfaceMesh::getField
(
	const List<pointIndexHit>& info,
	labelList& values
) const
{
	if (foundObject<triSurfaceLabelField>("values"))
	{
		values.setSize(info.size());

		const triSurfaceLabelField& fld = lookupObject<triSurfaceLabelField>
		(
			"values"
		);

		forAll(info, i)
		{
			if (info[i].hit())
			{
				values[i] = fld[info[i].index()];
			}
		}
	}
}


void Foam::triSurfaceMesh::getVolumeType
(
	const pointField& points,
	List<volumeType>& volType
) const
{
	volType.setSize(points.size());

	scalar oldTol = indexedOctree<treeDataTriSurface>::perturbTol();
	indexedOctree<treeDataTriSurface>::perturbTol() = tolerance_;

	forAll(points, pointI)
	{
		const point& pt = points[pointI];

		// - use cached volume type per each tree node
		// - cheat conversion since same values
		volType[pointI] = static_cast<volumeType>(tree().getVolumeType(pt));
	}

	indexedOctree<treeDataTriSurface>::perturbTol() = oldTol;
}


//- Write using given format, version and compression
bool Foam::triSurfaceMesh::writeObject
(
	IOstream::streamFormat fmt,
	IOstream::versionNumber ver,
	IOstream::compressionType cmp
) const
{
	fileName fullPath(searchableSurface::objectPath());

	if (!mkDir(fullPath.path()))
	{
		return false;
	}

	triSurface::write(fullPath);

	if (!isFile(fullPath))
	{
		return false;
	}

	//return objectRegistry::writeObject(fmt, ver, cmp);
	return true;
}


// ************************************************************************* //
