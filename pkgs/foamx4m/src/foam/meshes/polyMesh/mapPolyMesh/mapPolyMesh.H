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
	Foam::mapPolyMesh

Description
	Class containing mesh-to-mesh mapping information after a change
	in polyMesh topology.

	General:
		- pointMap/faceMap/cellMap: \n
		  from current mesh back to previous mesh.
		  (so to 'pull' the information onto the current mesh)
		- reversePointMap/faceMap/cellMap: \n
		  from previous mesh to current. (so to 'push' information)

	In the topology change points/faces/cells
	- can be unchanged. (faces might be renumbered though)
	- can be removed (into nothing)
	- can be removed into/merged with existing same entity
	  (so point merged with other point, face with other face, cell with
	   other cell. Note that probably only cell with cell is relevant)
	- can be added from existing same 'master' entity
	  (so point from point, face from face and cell from cell)
	- can be inflated: face out of edge or point,
		cell out of face, edge or point.
	- can be appended: added 'out of nothing'.

	All this information is nessecary to correctly map fields.

	@par points

	- unchanged:
		- pointMap[pointI] contains old point label
		- reversePointMap[oldPointI] contains new point label
	- removed:
		- reversePointMap[oldPointI] contains -1
	- merged into point:
		- reversePointMap[oldPointI] contains <-1 : -newPointI-2
		- pointMap[pointI] contains the old master point label
		- pointsFromPoints gives for pointI all the old point labels
		  (including the old master point!)
	- added-from-same:
		- pointMap[pointI] contains the old master point label
	- appended:
		- pointMap[pointI] contains -1

	@par faces

	- unchanged:
		- faceMap[faceI] contains old face label
		- reverseFaceMap[oldFaceI] contains new face label
	- removed:
		- reverseFaceMap[oldFaceI] contains -1
	- merged into face:
		- reverseFaceMap[oldFaceI] contains <-1 : -newFaceI-2
		- faceMap[faceI] contains the old master face label
		- facesFromFaces gives for faceI all the old face labels
		  (including the old master face!)
	- added-from-same:
		- faceMap[faceI] contains the old master face label
	- inflated-from-edge:
		- faceMap[faceI] contains -1
		- facesFromEdges contains an entry with
			- faceI
			- list of faces(*) on old mesh that connected to the old edge
	- inflated-from-point:
		- faceMap[faceI] contains -1
		- facesFromPoints contains an entry with
			- faceI
			- list of faces(*) on old mesh that connected to the old point
	- appended:
		- faceMap[faceI] contains -1

	Note (*) \n
	   if the newly inflated face is a boundary face the list of faces will
	   only be boundary faces; if the new face is an internal face they
	   will only be internal faces.

	@par cells

	- unchanged:
		- cellMap[cellI] contains old cell label
		- reverseCellMap[oldCellI] contains new cell label
	- removed:
		- reverseCellMap[oldCellI] contains -1
	- merged into cell:
		- reverseCellMap[oldCellI] contains <-1 : -newCellI-2
		- cellMap[cellI] contains the old master cell label
		- cellsFromCells gives for cellI all the old cell labels
		  (including the old master cell!)
	- added-from-same:
		- cellMap[cellI] contains the old master cell label
	- inflated-from-face:
		- cellMap[cellI] contains -1
		- cellsFromFaces contains an entry with
			- cellI
			- list of cells on old mesh that connected to the old face
	- inflated-from-edge:
		- cellMap[cellI] contains -1
		- cellsFromEdges contains an entry with
			- cellI
			- list of cells on old mesh that connected to the old edge
	- inflated-from-point:
		- cellMap[cellI] contains -1
		- cellsFromPoints contains an entry with
			- cellI
			- list of cells on old mesh that connected to the old point
	- appended:
		- cellMap[cellI] contains -1


SourceFiles
	mapPolyMesh.C

\*---------------------------------------------------------------------------*/

#ifndef mapPolyMesh_H
#define mapPolyMesh_H

#include "boolList.H"
#include "labelList.H"
#include "objectMap.H"
#include "pointField.H"
#include "HashSet.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMesh;


class mapPolyMesh
:
	public refCount
{
	// Private data

		//- Reference to polyMesh
		const polyMesh& mesh_;

		//- Is the mesh morphing?
		//  Changed signalling.  HJ, 27/Nov/2009
		const bool morphing_;

		//- Number of old live points
		const label nOldPoints_;

		//- Number of old live faces
		const label nOldFaces_;

		//- Number of old live cells
		const label nOldCells_;

		//- Old point map.
		//  Contains the old point label for all new points.
		//  - for preserved points this is the old point label.
		//  - for added points this is the master point ID
		//  - for points added with no master, this is -1
		//  Size of the list equals the size of new points
		const labelList pointMap_;

		//- Points resulting from merging points
		const List<objectMap> pointsFromPointsMap_;

		//- Old face map.
		//  Contains a list of old face labels for every new face.
		//  Size of the list equals the number of new faces
		//  - for preserved faces this is the old face label.
		//  - for faces added from faces this is the master face ID
		//  - for faces added with no master, this is -1
		//  - for faces added from points or edges, this is -1
		const labelList faceMap_;

		//- Faces inflated from points
		const List<objectMap> facesFromPointsMap_;

		//- Faces inflated from edges
		const List<objectMap> facesFromEdgesMap_;

		//- Faces resulting from merging faces
		const List<objectMap> facesFromFacesMap_;

		//- Old cell map.
		//  Contains old cell label for all preserved cells.
		//  Size of the list equals the number or preserved cells
		const labelList cellMap_;

		//- Cells inflated from points
		const List<objectMap> cellsFromPointsMap_;

		//- Cells inflated from edges
		const List<objectMap> cellsFromEdgesMap_;

		//- Cells inflated from faces
		const List<objectMap> cellsFromFacesMap_;

		//- Cells resulting from merging cells
		const List<objectMap> cellsFromCellsMap_;

		//- Reverse point map
		const labelList reversePointMap_;

		//- Reverse face map
		const labelList reverseFaceMap_;

		//- Reverse cell map
		const labelList reverseCellMap_;

		//- Map of flipped face flux faces
		const labelHashSet flipFaceFlux_;

		//- Patch mesh point renumbering
		const labelListList patchPointMap_;

		//- Point zone renumbering
		//  For every preserved point in zone give the old position.
		//  For added points, the index is set to -1
		const labelListList pointZoneMap_;

		//- Face zone point renumbering
		//  For every preserved point in zone give the old position.
		//  For added points, the index is set to -1
		const labelListList faceZonePointMap_;

		//- Face zone face renumbering
		//  For every preserved face in zone give the old position.
		//  For added faces, the index is set to -1
		const labelListList faceZoneFaceMap_;

		//- Cell zone renumbering
		//  For every preserved cell in zone give the old position.
		//  For added cells, the index is set to -1
		const labelListList cellZoneMap_;

		//- Reset patch flag
		//  True for the rebuilt patch (no mapping)
		boolList resetPatchFlag_;

		//- Pre-motion point positions.
		//  This specifies the correct way of blowing up zero-volume objects
		const pointField preMotionPoints_;

		//- List of the old patch sizes
		labelList oldPatchSizes_;

		//- List of the old patch start labels
		const labelList oldPatchStarts_;

		//- List of numbers of mesh points per old patch
		const labelList oldPatchNMeshPoints_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		mapPolyMesh(const mapPolyMesh&);

		//- Disallow default bitwise assignment
		void operator=(const mapPolyMesh&);


public:

	// Constructors

		//- Construct empty: no update
		mapPolyMesh(const polyMesh& mesh);

		//- Construct from components
		mapPolyMesh
		(
			const polyMesh& mesh,
			const label nOldPoints,
			const label nOldFaces,
			const label nOldCells,

			const labelList& pointMap,
			const List<objectMap>& pointsFromPoints,

			const labelList& faceMap,
			const List<objectMap>& facesFromPoints,
			const List<objectMap>& facesFromEdges,
			const List<objectMap>& facesFromFaces,

			const labelList& cellMap,
			const List<objectMap>& cellsFromPoints,
			const List<objectMap>& cellsFromEdges,
			const List<objectMap>& cellsFromFaces,
			const List<objectMap>& cellsFromCells,

			const labelList& reversePointMap,
			const labelList& reverseFaceMap,
			const labelList& reverseCellMap,

			const labelHashSet& flipFaceFlux,

			const labelListList& patchPointMap,

			const labelListList& pointZoneMap,
			const labelListList& faceZonePointMap,
			const labelListList& faceZoneFaceMap,
			const labelListList& cellZoneMap,

			const boolList& resetPatchFlag,

			const pointField& preMotionPoints,
			const labelList& oldPatchStarts,
			const labelList& oldPatchNMeshPoints
		);

		//- Construct from components and optionally reuse storage
		mapPolyMesh
		(
			const polyMesh& mesh,
			const label nOldPoints,
			const label nOldFaces,
			const label nOldCells,

			labelList& pointMap,
			List<objectMap>& pointsFromPoints,

			labelList& faceMap,
			List<objectMap>& facesFromPoints,
			List<objectMap>& facesFromEdges,
			List<objectMap>& facesFromFaces,

			labelList& cellMap,
			List<objectMap>& cellsFromPoints,
			List<objectMap>& cellsFromEdges,
			List<objectMap>& cellsFromFaces,
			List<objectMap>& cellsFromCells,

			labelList& reversePointMap,
			labelList& reverseFaceMap,
			labelList& reverseCellMap,

			labelHashSet& flipFaceFlux,

			labelListList& patchPointMap,

			labelListList& pointZoneMap,
			labelListList& faceZonePointMap,
			labelListList& faceZoneFaceMap,
			labelListList& cellZoneMap,

			boolList& resetPatchFlag,

			pointField& preMotionPoints,
			labelList& oldPatchStarts,
			labelList& oldPatchNMeshPoints,
			const bool reUse
		);


	// Member Functions

		// Access

			//- Return polyMesh
			const polyMesh& mesh() const
			{
				return mesh_;
			}

			//- Is the mesh morphing?
			bool morphing() const
			{
				return morphing_;
			}

			//- Number of old points
			label nOldPoints() const
			{
				return nOldPoints_;
			}

			//- Number of old internal faces
			label nOldInternalFaces() const
			{
				return oldPatchStarts_[0];
			}

			//- Number of old faces
			label nOldFaces() const
			{
				return nOldFaces_;
			}

			//- Number of old cells
			label nOldCells() const
			{
				return nOldCells_;
			}

			//- Old point map.
			//  Contains the old point label for all new points.
			//  For preserved points this is the old point label.
			//  For added points this is the master point ID
			const labelList& pointMap() const
			{
				return pointMap_;
			}

			//- Points originating from points
			const List<objectMap>& pointsFromPointsMap() const
			{
				return pointsFromPointsMap_;
			}

			//- Old face map.
			//  Contains a list of old face labels for every new face.
			//  Warning: this map contains invalid entries for new faces
			const labelList& faceMap() const
			{
				return faceMap_;
			}

			//- Faces inflated from points
			const List<objectMap>& facesFromPointsMap() const
			{
				return facesFromPointsMap_;
			}

			//- Faces inflated from edges
			const List<objectMap>& facesFromEdgesMap() const
			{
				return facesFromEdgesMap_;
			}

			//- Faces originating from faces
			const List<objectMap>& facesFromFacesMap() const
			{
				return facesFromFacesMap_;
			}

			//- Old cell map.
			//  Contains old cell label for all preserved cells.
			const labelList& cellMap() const
			{
				return cellMap_;
			}

			//- Cells inflated from points
			const List<objectMap>& cellsFromPointsMap() const
			{
				return cellsFromPointsMap_;
			}

			//- Cells inflated from edges
			const List<objectMap>& cellsFromEdgesMap() const
			{
				return cellsFromEdgesMap_;
			}

			//- Cells inflated from faces
			const List<objectMap>& cellsFromFacesMap() const
			{
				return cellsFromFacesMap_;
			}

			//- Cells originating from cells
			const List<objectMap>& cellsFromCellsMap() const
			{
				return cellsFromCellsMap_;
			}


			// Reverse maps

				//- Reverse point map
				//  Contains new point label for all old and added points
				const labelList& reversePointMap() const
				{
					return reversePointMap_;
				}

				//- If point is removed return point (on new mesh) it merged
				//  into
				label mergedPoint(const label oldPointI) const
				{
					label i = reversePointMap_[oldPointI];

					if (i == -1)
					{
					    return i;
					}
					else if (i < -1)
					{
					    return -i - 2;
					}
					else
					{
					    FatalErrorIn("mergedPoint(const label) const")
					        << "old point label " << oldPointI
					        << " has reverseMap " << i << endl
					        << "Only call mergedPoint for removed points."
					        << abort(FatalError);

					    return -1;
					}
				}

				//- Reverse face map
				//  Contains new face label for all old and added faces
				const labelList& reverseFaceMap() const
				{
					return reverseFaceMap_;
				}

				//- If face is removed return face (on new mesh) it merged into
				label mergedFace(const label oldFaceI) const
				{
					label i = reverseFaceMap_[oldFaceI];

					if (i == -1)
					{
					    return i;
					}
					else if (i < -1)
					{
					    return -i - 2;
					}
					else
					{
					    FatalErrorIn("mergedFace(const label) const")
					        << "old face label " << oldFaceI
					        << " has reverseMap " << i << endl
					        << "Only call mergedFace for removed faces."
					        << abort(FatalError);

					    return -1;
					}
				}

				//- Reverse cell map
				//  Contains new cell label for all old and added cells
				const labelList& reverseCellMap() const
				{
					return reverseCellMap_;
				}

				//- If cell is removed return cell (on new mesh) it merged into
				label mergedCell(const label oldCellI) const
				{
					label i = reverseCellMap_[oldCellI];

					if (i == -1)
					{
					    return i;
					}
					else if (i < -1)
					{
					    return -i - 2;
					}
					else
					{
					    FatalErrorIn("mergedCell(const label) const")
					        << "old cell label " << oldCellI
					        << " has reverseMap " << i << endl
					        << "Only call mergedCell for removed cells."
					        << abort(FatalError);

					    return -1;
					}
				}

				//- Map of flipped face flux faces
				const labelHashSet& flipFaceFlux() const
				{
					return flipFaceFlux_;
				}

				//- Patch point renumbering
				//  For every preserved point on a patch give the old position.
				//  For added points, the index is set to -1
				const labelListList& patchPointMap() const
				{
					return patchPointMap_;
				}


			// Zone mapping

				//- Point zone renumbering
				//  For every preserved point in zone give the old position.
				//  For added points, the index is set to -1
				const labelListList& pointZoneMap() const
				{
					return pointZoneMap_;
				}

				//- Face zone point renumbering
				//  For every preserved point in zone give the old position.
				//  For added points, the index is set to -1
				const labelListList& faceZonePointMap() const
				{
					return faceZonePointMap_;
				}

				//- Face zone face renumbering
				//  For every preserved face in zone give the old position.
				//  For added faces, the index is set to -1
				const labelListList& faceZoneFaceMap() const
				{
					return faceZoneFaceMap_;
				}

				//- Cell zone renumbering
				//  For every preserved cell in zone give the old position.
				//  For added cells, the index is set to -1
				const labelListList& cellZoneMap() const
				{
					return cellZoneMap_;
				}

				//- Reset patch flag
				//  True for the rebuilt patch (no mapping)
				const boolList& resetPatchFlag() const
				{
					return resetPatchFlag_;
				}


			//- Pre-motion point positions.
			//  This specifies the correct way of blowing up
			//  zero-volume objects
			const pointField& preMotionPoints() const
			{
				return preMotionPoints_;
			}

			//- Has valid preMotionPoints?
			bool hasMotionPoints() const
			{
				return preMotionPoints_.size() > 0;
			}


			//- Return list of the old patch sizes
			const labelList& oldPatchSizes() const
			{
				return oldPatchSizes_;
			}

			//- Return list of the old patch start labels
			const labelList& oldPatchStarts() const
			{
				return oldPatchStarts_;
			}

			//- Return numbers of mesh points per old patch
			const labelList& oldPatchNMeshPoints() const
			{
				return oldPatchNMeshPoints_;
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
