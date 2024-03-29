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
	Foam::triSurfaceMesh

Description
	IOoject and searching on triSurface

	Note: when constructing from dictionary has optional parameters:
		- scale     : scaling factor.
		- tolerance : relative tolerance for doing intersections
					  (see triangle::intersection)

SourceFiles
	triSurfaceMesh.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceMesh_H
#define triSurfaceMesh_H

#include "treeBoundBox.H"
#include "searchableSurface.H"
#include "objectRegistry.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"
#include "octreeDataTriSurfaceTreeLeaf.H"
#include "treeDataEdge.H"
#include "EdgeMap.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class triSurfaceMesh
:
	public searchableSurface,
	public objectRegistry,      // so we can store fields
	public triSurface
{
private:

	// Private member data

		//- Optional tolerance to use in searches
		scalar tolerance_;

		//- Optional max tree depth of octree
		label maxTreeDepth_;

		//- Search tree (triangles)
		mutable autoPtr<indexedOctree<treeDataTriSurface> > tree_;

		//- Search tree for boundary edges.
		mutable autoPtr<indexedOctree<treeDataEdge> > edgeTree_;

		//- Names of regions
		mutable wordList regions_;

		//- Is surface closed
		mutable label surfaceClosed_;


	// Private Member Functions

		////- Helper: find instance of files without header
		//static word findRawInstance
		//(
		//    const Time&,
		//    const fileName&,
		//    const word&
		//);

		//- Check file existence
		static const fileName& checkFile
		(
			const fileName& fName,
			const fileName& objectName
		);

		//- Helper function for isSurfaceClosed
		static bool addFaceToEdge
		(
			const edge&,
			EdgeMap<label>&
		);

		//- Check whether surface is closed without calculating any permanent
		//  addressing.
		bool isSurfaceClosed() const;

		//- Steps to next intersection. Adds smallVec and starts tracking
		//  from there.
		static void getNextIntersections
		(
			const indexedOctree<treeDataTriSurface>& octree,
			const point& start,
			const point& end,
			const vector& smallVec,
			DynamicList<pointIndexHit, 1, 1>& hits
		);

		//- Disallow default bitwise copy construct
		triSurfaceMesh(const triSurfaceMesh&);

		//- Disallow default bitwise assignment
		void operator=(const triSurfaceMesh&);


protected:

		//- Calculate (number of)used points and their bounding box
		void calcBounds(boundBox& bb, label& nPoints) const;

public:

	//- Runtime type information
	TypeName("triSurfaceMesh");


	// Constructors

		//- Construct from triSurface
		triSurfaceMesh(const IOobject&, const triSurface&);

		//- Construct read.
		triSurfaceMesh(const IOobject& io);

		//- Construct from IO and dictionary (used by searchableSurface).
		//  Dictionary may contain a 'scale' entry (eg, 0.001: mm -> m)
		triSurfaceMesh
		(
			const IOobject& io,
			const dictionary& dict
		);


	// Destructor

		virtual ~triSurfaceMesh();

		//- Clear storage
		void clearOut();


	// Member Functions

		//- Move points
		virtual void movePoints(const pointField&);

		//- Demand driven construction of octree
		const indexedOctree<treeDataTriSurface>& tree() const;

		//- Demand driven contruction of octree for boundary edges
		const indexedOctree<treeDataEdge>& edgeTree() const;


		// searchableSurface implementation

			virtual const wordList& regions() const;

			//- Whether supports volume type below. I.e. whether is closed.
			virtual bool hasVolumeType() const;

			//- Range of local indices that can be returned.
			virtual label size() const
			{
				return triSurface::size();
			}

			//- Get representative set of element coordinates
			//  Usually the element centres (should be of length size()).
			virtual pointField coordinates() const;

			virtual void findNearest
			(
				const pointField& sample,
				const scalarField& nearestDistSqr,
				List<pointIndexHit>&
			) const;

			virtual void findLine
			(
				const pointField& start,
				const pointField& end,
				List<pointIndexHit>&
			) const;

			virtual void findLineAny
			(
				const pointField& start,
				const pointField& end,
				List<pointIndexHit>&
			) const;

			//- Get all intersections in order from start to end.
			virtual void findLineAll
			(
				const pointField& start,
				const pointField& end,
				List<List<pointIndexHit> >&
			) const;

			//- From a set of points and indices get the region
			virtual void getRegion
			(
				const List<pointIndexHit>&,
				labelList& region
			) const;

			//- From a set of points and indices get the normal
			virtual void getNormal
			(
				const List<pointIndexHit>&,
				vectorField& normal
			) const;

			//- Determine type (inside/outside/mixed) for point. unknown if
			//  cannot be determined (e.g. non-manifold surface)
			virtual void getVolumeType
			(
				const pointField&,
				List<volumeType>&
			) const;


		// Other

			//- Set bounds of surface. Bounds currently set as list of
			//  bounding boxes. The bounds are hints to the surface as for
			//  the range of queries it can expect. faceMap/pointMap can be
			//  set if the surface has done any redistribution.
			virtual void distribute
			(
				const List<treeBoundBox>&,
				const bool keepNonLocal,
				autoPtr<mapDistribute>& faceMap,
				autoPtr<mapDistribute>& pointMap
			)
			{}

			//- WIP. Store element-wise field.
			virtual void setField(const labelList& values);

			//- WIP. From a set of hits (points and
			//  indices) get the specified field. Misses do not get set.
			virtual void getField
			(
				const List<pointIndexHit>&,
				labelList&
			) const;


		// regIOobject implementation

			bool writeData(Ostream&) const
			{
				notImplemented("triSurfaceMesh::writeData(Ostream&) const");

				return false;
			}

			//- Write using given format, version and compression
			virtual bool writeObject
			(
				IOstream::streamFormat fmt,
				IOstream::versionNumber ver,
				IOstream::compressionType cmp
			) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
