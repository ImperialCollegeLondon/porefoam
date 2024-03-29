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
	Foam::enrichedPatch

Description
	The enriched patch contains a double set of faces from the two
	sides of the sliding interface before the cutting.

	The patch basically consists of two over-lapping sets of faces sitting
	on a common point support, where every edge may be shared by more than
	2 faces.  The patch points are collected in a map.  Additional
	information needed for cutting is the point insertion into every edge
	of master and slave.

	Note:
	If new points are created during master-slave edge cutting, they
	should be registred with the pointMap.

Author
	Hrvoje Jasak, Nabla Ltd. and Wikki Ltd.  All rights reserved.
	Copyright Hrvoje Jasak

SourceFiles
	enrichedPatch.C
	enrichedPatchCutFaces.C
	enrichedPatchFaces.C
	enrichedPatchPointMap.C
	enrichedPatchPointMergeMap.C
	enrichedPatchPointPoints.C

\*---------------------------------------------------------------------------*/

#ifndef enrichedPatch_H
#define enrichedPatch_H

#include "primitiveFacePatch.H"
#include "Map.H"
#include "point.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


class enrichedPatch
{
	// Private data

		//- Reference to master patch
		const primitiveFacePatch& masterPatch_;

		//- Reference to slave patch
		const primitiveFacePatch& slavePatch_;

		//- Map of points supporting patch faces
		mutable Map<point> pointMap_;

		//- Has the point map been completed?
		mutable bool pointMapComplete_;

		//- Map of point merges
		mutable Map<label> pointMergeMap_;

		//- Slave point point hits
		const labelList& slavePointPointHits_;

		//- Slave point edge hits
		const labelList& slavePointEdgeHits_;

		//- Slave point face hits
		const List<objectHit>& slavePointFaceHits_;


	// Demand-driven private data

		//- Enriched patch
		mutable faceList* enrichedFacesPtr_;

		//- Mesh points
		mutable labelList* meshPointsPtr_;

		//- Local faces
		mutable faceList* localFacesPtr_;

		//- Local points
		mutable pointField* localPointsPtr_;

		//- Point-point addressing
		mutable labelListList* pointPointsPtr_;

		// Master point addressing
		mutable Map<labelList>* masterPointFacesPtr_;


		// Cut faces and addressing

			//- Cut faces
			mutable faceList* cutFacesPtr_;

			//- Cut face master
			//  - the face on the master patch for internal faces
			//  - the creator face for boundary face
			mutable labelList* cutFaceMasterPtr_;

			//- Cut face slave
			//  - the face on the slave patch for internal faces
			//  - -1 for boundary face
			mutable labelList* cutFaceSlavePtr_;

			//- Local cut faces
			mutable faceList* localCutFacesPtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		enrichedPatch(const enrichedPatch&);

		//- Disallow default bitwise assignment
		void operator=(const enrichedPatch&);

		// Creation of demand-driven private data

			//- Calculate point merge map
			void calcPointMergeMap() const;

			//- Calculate mesh points
			void calcMeshPoints() const;

			//- Calculate local points
			void calcLocalPoints() const;

			//- Calculate local faces
			void calcLocalFaces() const;

			//- Calculate point-point addressing
			void calcPointPoints() const;

			//- Calculate master point addressing
			void calcMasterPointFaces() const;

			//- Calculate cut faces
			void calcCutFaces() const;

			//- Calculate local cut faces
			void calcLocalCutFaces() const;

			//- Clear cut faces
			void clearCutFaces();

			//- Clear out demand-driven data
			void clearOut();


	// Static data members

		//- Estimated ratio of original-to-enriched face size
		static const label enrichedFaceRatio_;

		//- Estimated number of master face hits by slave points
		static const label nFaceHits_;

		//- Size of face on which the check is forced
		static const label maxFaceSizeDebug_;

		//- Out of plane tolerance used to reject points in right-hand walk
		static const scalar outOfPlaneTol_;


public:

	// Static data members
	ClassName("enrichedPatch");

	// Constructors

		//- Construct from components
		enrichedPatch
		(
			const primitiveFacePatch& masterPatch,
			const primitiveFacePatch& slavePatch,
			const labelList& slavePointPointHits,
			const labelList& slavePointEdgeHits,
			const List<objectHit>& slavePointFaceHits
		);


	// Destructor

		~enrichedPatch();


	// Member Functions

		// Access

			//- Return non-const access to point map to add points
			Map<point>& pointMap();

			//- Return map of points
			const Map<point>& pointMap() const;

			//- Return map of point merges
			Map<label>& pointMergeMap()
			{
				return pointMergeMap_;
			}

			//- Return map of point merges
			const Map<label>& pointMergeMap() const
			{
				return pointMergeMap_;
			}

			//- Complete point map
			void completePointMap() const;


		// Topological data

			//- Calculate enriched faces
			void calcEnrichedFaces
			(
				const labelListList& pointsIntoMasterEdges,
				const labelListList& pointsIntoSlaveEdges,
				const pointField& projectedSlavePoints
			);

			//- Return enriched faces
			const faceList& enrichedFaces() const;

			//- Return mesh points
			const labelList& meshPoints() const;

			//- Return local faces
			const faceList& localFaces() const;

			//- Return local points
			const pointField& localPoints() const;

			//- Return point-point addressing
			const labelListList& pointPoints() const;

			//- Master point face addressing
			const Map<labelList>& masterPointFaces() const;


		// Cut faces

			//- Return list of cut faces
			const faceList& cutFaces() const;

			//- Return cut face master list
			const labelList& cutFaceMaster() const;

			//- Return cut face slave list
			const labelList& cutFaceSlave() const;

			//- Return local cut faces
			const faceList& localCutFaces() const;


		//- Check if the patch is fully supported
		bool checkSupport() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
