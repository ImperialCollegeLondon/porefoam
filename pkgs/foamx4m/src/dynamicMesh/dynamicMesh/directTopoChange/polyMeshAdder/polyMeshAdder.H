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
	Foam::polyMeshAdder

Description
	Adds two meshes without using any polyMesh morphing.

	Gets faces to couple as faceCoupleInfo which is list of faces on both
	meshes. Holds map from last mesh addition.

SourceFiles
	polyMeshAdder.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshAdder_H
#define polyMeshAdder_H

#include "autoPtr.H"
#include "polyMesh.H"
#include "mapAddedPolyMesh.H"
#include "faceCoupleInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class IOobject;
class directTopoChange;


class polyMeshAdder
{

private:
	// Private data


	// Private Member Functions

		//- Append all elements of a list to a DynamicList
		template<class T>
		static void append(const List<T>&, DynamicList<T>&);

		//- Append all mapped elements of a list to a DynamicList
		static void append
		(
			const labelList& map,
			const labelList& lst,
			dynamicLabelList&
		);

		//- Index of patch in allPatches. Add if nonexisting.
		static label patchIndex
		(
			const polyPatch&,
			DynamicList<word>& allPatchNames,
			DynamicList<word>& allPatchTypes
		);

		//- Index of zone in all zones
		static label zoneIndex(const word&, DynamicList<word>&);

		static void mergePatchNames
		(
			const polyBoundaryMesh& patches0,
			const polyBoundaryMesh& patches1,

			DynamicList<word>& allPatchNames,
			DynamicList<word>& allPatchTypes,

			labelList& from1ToAllPatches,
			labelList& fromAllTo1Patches
		);

		//- Get starts of patches
		static labelList getPatchStarts(const polyBoundaryMesh&);
		//- Get sizes of patches
		static labelList getPatchSizes(const polyBoundaryMesh&);

		static List<polyPatch*> combinePatches
		(
			const polyMesh& mesh0,
			const polyMesh& mesh1,
			const polyBoundaryMesh& allBoundaryMesh,
			const label nAllPatches,
			const labelList& fromAllTo1Patches,

			const label nInternalFaces,
			const labelList& nFaces,

			labelList& from0ToAllPatches,
			labelList& from1ToAllPatches
		);

		//- Determine order for internalFaces to be upper-triangular.
		//  Does not change order of external faces.
		static labelList getFaceOrder
		(
			const cellList& cells,
			const label nInternalFaces,
			const labelList& owner,
			const labelList& neighbour
		);

		//- Extends face f with split points. cutEdgeToPoints gives for every
		// edge the points introduced inbetween the endpoints.
		static void insertVertices
		(
			const edgeLookup& cutEdgeToPoints,
			const Map<label>& meshToMaster,
			const labelList& masterToCutPoints,
			const face& masterFace,

			dynamicLabelList& workFace,
			face& f
		);

		//- Merges primitives of two meshes.
		static void mergePrimitives
		(
			const polyMesh& mesh0,
			const polyMesh& mesh1,
			const faceCoupleInfo& coupleInfo,

			const label nAllPatches,
			const labelList& fromAllTo1Patches,
			const labelList& from1ToAllPatches,

			pointField& points,
			labelList& from0ToAllPoints,
			labelList& from1ToAllPoints,

			faceList& allFaces,
			labelList& allOwner,
			labelList& allNeighbour,
			label& nInternalFaces,
			labelList& nFacesPerPatch,
			label& nCells,

			labelList& from0ToAllFaces,
			labelList& from1ToAllFaces,
			labelList& from1ToAllCells
		);

		//- Merge point zones
		static void mergePointZones
		(
			const pointZoneMesh& pz0,
			const pointZoneMesh& pz1,
			const labelList& from0ToAllPoints,
			const labelList& from1ToAllPoints,

			DynamicList<word>& zoneNames,
			labelList& from1ToAll,
			List<dynamicLabelList >& pzPoints
		);

		//- Merge face zones
		static void mergeFaceZones
		(
			const faceZoneMesh& fz0,
			const faceZoneMesh& fz1,
			const labelList& from0ToAllFaces,
			const labelList& from1ToAllFaces,

			DynamicList<word>& zoneNames,
			labelList& from1ToAll,
			List<dynamicLabelList >& fzFaces,
			List<DynamicList<bool> >& fzFlips
		);

		//- Merge cell zones
		static void mergeCellZones
		(
			const cellZoneMesh& cz0,
			const cellZoneMesh& cz1,
			const labelList& from1ToAllCells,

			DynamicList<word>& zoneNames,
			labelList& from1ToAll,
			List<dynamicLabelList >& czCells
		);

		//- Merge point/face/cell zone information
		static void mergeZones
		(
			const polyMesh& mesh0,
			const polyMesh& mesh1,
			const labelList& from0ToAllPoints,
			const labelList& from0ToAllFaces,
			const labelList& from1ToAllPoints,
			const labelList& from1ToAllFaces,
			const labelList& from1ToAllCells,

			DynamicList<word>& pointZoneNames,
			List<dynamicLabelList >& pzPoints,

			DynamicList<word>& faceZoneNames,
			List<dynamicLabelList >& fzFaces,
			List<DynamicList<bool> >& fzFlips,

			DynamicList<word>& cellZoneNames,
			List<dynamicLabelList >& czCells
		);

		//- Create new zones and add to new mesh.
		static void addZones
		(
			const DynamicList<word>& pointZoneNames,
			const List<dynamicLabelList >& pzPoints,

			const DynamicList<word>& faceZoneNames,
			const List<dynamicLabelList >& fzFaces,
			const List<DynamicList<bool> >& fzFlips,

			const DynamicList<word>& cellZoneNames,
			const List<dynamicLabelList >& czCells,

			polyMesh& mesh
		);

public:


	// Member Functions


		//- Add two polyMeshes. Returns new polyMesh and map construct.
		static autoPtr<polyMesh> add
		(
			const IOobject& io,
			const polyMesh& mesh0,
			const polyMesh& mesh1,
			const faceCoupleInfo& coupleInfo,
			autoPtr<mapAddedPolyMesh>& mapPtr
		);

		//- Inplace add mesh to polyMesh. Returns map construct.
		static autoPtr<mapAddedPolyMesh> add
		(
			polyMesh& mesh0,
			const polyMesh& mesh1,
			const faceCoupleInfo& coupleInfo,
			const bool validBoundary = true
		);

		// Point merging

			//- Find topologically&geometrically shared points.
			//  - should only be called for parallel correct mesh
			//    (since uses mesh.globalData)
			//  - returns Map from point to master point (all in mesh point
			//    labels) for any sets of points that need to be merged.
			static Map<label> findSharedPoints
			(
				const polyMesh&,
				const scalar mergeTol
			);

			//- Helper: Merge points.
			//  - Gets map from point to destination point
			//  - Removes all points that don't map to themselves
			//  - Modifies all faces that use the points to be removed.
			static void mergePoints
			(
				const polyMesh&,
				const Map<label>& pointToMaster,
				directTopoChange& meshMod
			);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "polyMeshAdderTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
