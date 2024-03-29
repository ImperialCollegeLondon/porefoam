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
	Foam::searchableSurfaceCollection

Description
	Set of transformed searchableSurfaces. Does not do boolean operations.
	So when meshing might find parts 'inside'.

SourceFiles
	searchableSurfaceCollection.C

\*---------------------------------------------------------------------------*/

#ifndef searchableSurfaceCollection_H
#define searchableSurfaceCollection_H

#include "searchableSurface.H"
#include "treeBoundBox.H"
#include "coordinateSystem.H"
#include "UPtrList.H"
#include "Switch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


class searchableSurfaceCollection
:
	public searchableSurface
{
private:

	// Private Member Data

		// Per instance data

			//- instance name
			wordList instance_;

			//- scaling vector
			vectorField scale_;

			//- transformation
			PtrList<coordinateSystem> transform_;

			UPtrList<searchableSurface> subGeom_;

			Switch mergeSubRegions_;

			//- offsets for indices coming from different surfaces
			//  (sized with size() of each surface)
			labelList indexOffset_;

		//- Region names
		mutable wordList regions_;
		//- From individual regions to collection regions
		mutable labelList regionOffset_;


	// Private Member Functions

		//- Find point nearest to sample. Updates minDistSqr. Sets nearestInfo
		//  and surface index
		void findNearest
		(
			const pointField& samples,
			scalarField& minDistSqr,
			List<pointIndexHit>& nearestInfo,
			labelList& nearestSurf
		) const;

		//- Sort hits into per-surface bins. Misses are rejected.
		//  Maintains map back to position
		void sortHits
		(
			const List<pointIndexHit>& info,
			List<List<pointIndexHit> >& surfInfo,
			labelListList& infoMap
		) const;


		//- Disallow default bitwise copy construct
		searchableSurfaceCollection(const searchableSurfaceCollection&);

		//- Disallow default bitwise assignment
		void operator=(const searchableSurfaceCollection&);


public:

	//- Runtime type information
	TypeName("searchableSurfaceCollection");


	// Constructors

		//- Construct from dictionary (used by searchableSurface)
		searchableSurfaceCollection
		(
			const IOobject& io,
			const dictionary& dict
		);

	// Destructor

		virtual ~searchableSurfaceCollection();


	// Member Functions

		virtual const wordList& regions() const;

		//- Whether supports volume type below
		virtual bool hasVolumeType() const
		{
			return false;
		}

		//- Range of local indices that can be returned.
		virtual label size() const;

		//- Get representative set of element coordinates
		//  Usually the element centres (should be of length size()).
		virtual pointField coordinates() const;


		// Multiple point queries.

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
			);

			//- WIP. Store element-wise field.
			virtual void setField(const labelList& values);

			//- WIP. From a set of hits (points and
			//  indices) get the specified field. Misses do not get set. Return
			//  empty field if not supported.
			virtual void getField(const List<pointIndexHit>&, labelList&) const;

		// regIOobject implementation

			bool writeData(Ostream&) const
			{
				notImplemented
				(
					"searchableSurfaceCollection::writeData(Ostream&) const"
				);
				return false;
			}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
