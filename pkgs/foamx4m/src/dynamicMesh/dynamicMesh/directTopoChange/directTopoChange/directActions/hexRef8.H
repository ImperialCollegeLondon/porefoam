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
	Foam::hexRef8

Description
	Refinement of (split) hexes using directTopoChange.

SourceFiles
	hexRef8.C

\*---------------------------------------------------------------------------*/

#ifndef hexRef8_H
#define hexRef8_H

#include "labelIOList.H"
#include "face.H"
#include "HashSet.H"
#include "DynamicList.H"
#include "primitivePatch.H"
#include "removeFaces.H"
#include "refinementHistory.H"
#include "PackedList.H"
#include "labelIOField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;
class directTopoChange;
class mapPolyMesh;
class mapDistributePolyMesh;


class hexRef8
{
	// Private data

		//- Reference to underlying mesh.
		const polyMesh& mesh_;

		//- Per cell the refinement level
		labelIOField cellLevel_;

		//- Per point the refinement level
		labelIOField pointLevel_;

		//- Typical edge length between unrefined points
		const scalar level0Edge_;

		//- Refinement history
		refinementHistory history_;

		//- Face remover engine
		removeFaces faceRemover_;

		//- Level of saved points
		Map<label> savedPointLevel_;

		//- Level of saved cells
		Map<label> savedCellLevel_;


	// Private Member Functions

		//- Reorder according to map.
		static void reorder
		(
			const labelList& map,
			const label len,
			const label null,
			labelList& elems
		);

		//- Get patch and zone info
		void getFaceInfo
		(
			const label faceI,
			label& patchID,
			label& zoneID,
			label& zoneFlip
		) const;

		//- Adds a face on top of existing faceI. Reverses if nessecary.
		label addFace
		(
			directTopoChange& meshMod,
			const label faceI,
			const face& newFace,
			const label own,
			const label nei
		) const;

		//- Adds internal face from point. No checks on reversal.
		label addInternalFace
		(
			directTopoChange& meshMod,
			const label meshFaceI,
			const label meshPointI,
			const face& newFace,
			const label own,
			const label nei
		) const;

		//- Modifies existing faceI for either new owner/neighbour or new face
		//  points. Reverses if nessecary.
		void modFace
		(
			directTopoChange& meshMod,
			const label faceI,
			const face& newFace,
			const label own,
			const label nei
		) const;

		scalar getLevel0EdgeLength() const;

		//- Get cell added to point of cellI (if any)
		label getAnchorCell
		(
			const labelListList& cellAnchorPoints,
			const labelListList& cellAddedCells,
			const label cellI,
			const label faceI,
			const label pointI
		) const;

		//- Get new owner and neighbour (in unspecified order) of pointI
		//  on faceI.
		void getFaceNeighbours
		(
			const labelListList& cellAnchorPoints,
			const labelListList& cellAddedCells,
			const label faceI,
			const label pointI,

			label& own,
			label& nei
		) const;


		//- Get index of minimum pointlevel.
		label findMinLevel(const labelList& f) const;

		//- Get maximum pointlevel.
		label findMaxLevel(const labelList& f) const;

		//- Count number of vertices <= anchorLevel
		label countAnchors(const labelList&, const label) const;

		//- Find index of point with wantedLevel, starting from fp.
		label findLevel
		(
			const face& f,
			const label startFp,
			const bool searchForward,
			const label wantedLevel
		) const;

		////- Print levels of list of points.
		//void printLevels(Ostream&, const labelList&) const;

		//- debug:check orientation of added internal face
		static void checkInternalOrientation
		(
			directTopoChange& meshMod,
			const label cellI,
			const label faceI,
			const point& ownPt,
			const point& neiPt,
			const face& newFace
		);

		//- debug:check orientation of new boundary face
		static void checkBoundaryOrientation
		(
			directTopoChange& meshMod,
			const label cellI,
			const label faceI,
			const point& ownPt,
			const point& boundaryPt,
			const face& newFace
		);

		//- If p0 and p1 are existing vertices check if edge is split and insert
		// splitPoint.
		void insertEdgeSplit
		(
			const labelList& edgeMidPoint,
			const label p0,
			const label p1,
			dynamicLabelList& verts
		) const;

		//- Store in maps correspondence from midpoint to anchors and faces.
		label storeMidPointInfo
		(
			const labelListList& cellAnchorPoints,
			const labelListList& cellAddedCells,
			const labelList& cellMidPoint,
			const labelList& edgeMidPoint,
			const label cellI,
			const label faceI,
			const bool faceOrder,
			const label midPointI,
			const label anchorPointI,
			const label faceMidPointI,

			Map<edge>& midPointToAnchors,
			Map<edge>& midPointToFaceMids,
			directTopoChange& meshMod
		) const;

		//- Create all internal faces from an unsplit face.
		void createInternalFromSplitFace
		(
			const labelListList& cellAnchorPoints,
			const labelListList& cellAddedCells,
			const labelList& cellMidPoint,
			const labelList& faceMidPoint,
			const labelList& edgeMidPoint,
			const label cellI,
			const label faceI,

			Map<edge>& midPointToAnchors,
			Map<edge>& midPointToFaceMids,
			directTopoChange& meshMod,
			label& nFacesAdded
		) const;

		//- Create all internal faces to split cellI into 8.
		void createInternalFaces
		(
			const labelListList& cellAnchorPoints,
			const labelListList& cellAddedCells,
			const labelList& cellMidPoint,
			const labelList& faceMidPoint,
			const labelList& faceAnchorLevel,
			const labelList& edgeMidPoint,
			const label cellI,
			directTopoChange& meshMod
		) const;

		//- Store vertices from startFp upto face split point.
		//  Used when splitting face into 4.
		void walkFaceToMid
		(
			const labelList& edgeMidPoint,
			const label cLevel,
			const label faceI,
			const label startFp,
			dynamicLabelList& faceVerts
		) const;

		//- Same as walkFaceToMid but now walk back.
		void walkFaceFromMid
		(
			const labelList& edgeMidPoint,
			const label cLevel,
			const label faceI,
			const label startFp,
			dynamicLabelList& faceVerts
		) const;

		//- Updates refineCell so consistent 2:1 refinement. Returns local
		//  number of cells changed.
		label faceConsistentRefinement
		(
			const bool maxSet,
			PackedList<1>& refineCell
		) const;

		//- Check wanted refinement for 2:1 consistency
		void checkWantedRefinementLevels(const labelList&) const;



		//- Disallow default bitwise copy construct
		hexRef8(const hexRef8&);

		//- Disallow default bitwise assignment
		void operator=(const hexRef8&);


public:

	//- Runtime type information
	ClassName("hexRef8");


	// Constructors

		//- Construct from mesh, read_if_present refinement data
		//  (from write below)
		hexRef8(const polyMesh& mesh);

		//- Construct from mesh and un/refinement data.
		hexRef8
		(
			const polyMesh& mesh,
			const labelField& cellLevel,
			const labelField& pointLevel,
			const refinementHistory& history
		);

		//- Construct from mesh and refinement data.
		hexRef8
		(
			const polyMesh& mesh,
			const labelField& cellLevel,
			const labelField& pointLevel
		);


	// Member Functions

		// Access

			const labelIOField& cellLevel() const
			{
				return cellLevel_;
			}

			const labelIOField& pointLevel() const
			{
				return pointLevel_;
			}

			const refinementHistory& history() const
			{
				return history_;
			}

			//- Typical edge length between unrefined points
			scalar level0EdgeLength() const
			{
				return level0Edge_;
			}

		// Refinement
			//- Gets level such that the face has four points <= level.
			label getAnchorLevel(const label faceI) const;

			//- Given valid mesh and current cell level and proposed
			//  cells to refine calculate any clashes (due to 2:1) and return
			//  ok list of cells to refine.
			//  Either adds cells to refine to set (maxSet = true) or
			//  removes cells to refine (maxSet = false)

			//- Helper:get points of a cell without using cellPoints addressing
			labelList cellPoints(const label cellI) const;

			//- Given valid mesh and current cell level and proposed
			//  cells to refine calculate any clashes (due to 2:1) and return
			//  ok list of cells to refine.
			//  Either adds cells to refine to set (maxSet = true) or
			//  removes cells to refine (maxSet = false)
			labelList consistentRefinement
			(
				const labelList& cellsToRefine,
				const bool maxSet
			) const;

			//- Like consistentRefinement but slower:
			//  - specify number of cells between consecutive refinement levels
			//    (consistentRefinement equivalent to 1)
			//  - specify max level difference between point-connected cells.
			//    (-1 to disable) Note that with normal 2:1 limitation
			//    (maxFaceDiff=1) there can be 8:1 size difference across point
			//    connected cells so maxPointDiff allows you to make that less.
			//  cellsToRefine : cells we're thinking about refining. It will
			//                  extend this set. All refinement levels will be
			//                  at least maxFaceDiff layers thick.
			//  facesToCheck  : additional faces where to implement the
			//                  maxFaceDiff thickness (usually only boundary
			//                  faces)
			labelList consistentSlowRefinement
			(
				const label maxFaceDiff,
				const labelList& cellsToRefine,
				const labelList& facesToCheck,
				const label maxPointDiff,
				const labelList& pointsToCheck
			) const;

			//- Like consistentSlowRefinement but uses different meshWave
			//  (proper distance instead of toplogical count). No point checks
			//  yet.
			labelList consistentSlowRefinement2
			(
				const label maxFaceDiff,
				const labelList& cellsToRefine,
				const labelList& facesToCheck
			) const;

			//- Insert refinement. All selected cells will be split into 8.
			//  Returns per element in cells the 8 cells they were split into.
			//  Guarantees that the 0th element is the original cell label.
			//  Mapping:
			//  -split cells: 7 new ones get added from original
			//  -split faces: original gets modified; 3 new ones get added
			//               from original
			//  -added internal faces: added from original cell face(if
			//   that was internal) or created out-of-nothing (so will not
			//   get mapped!). Note: could make this inflate from point but
			//   that will allocate interpolation.
			//  -points added to split edge: added from edge start()
			//  -midpoints added: added from cellPoints[0].
			labelListList setRefinement
			(
				const labelList& cells,
				directTopoChange&
			);

			//- Update local numbering for changed mesh.
			void updateMesh(const mapPolyMesh&);


			// Restoring : is where other processes delete and reinsert data.
			//             These callbacks allow this to restore the cellLevel
			//             and pointLevel for reintroduced points.
			//             Is not related to undoing my refinement

				//- Signal points/face/cells for which to store data
				void storeData
				(
					const labelList& pointsToStore,
					const labelList& facesToStore,
					const labelList& cellsToStore
				);

				//- Update local numbering + undo
				//  Data to restore given as new pointlabel + stored pointlabel
				//  (i.e. what was in pointsToStore)
				void updateMesh
				(
					const mapPolyMesh&,
					const Map<label>& pointsToRestore,
					const Map<label>& facesToRestore,
					const Map<label>& cellsToRestore
				);


			//- Update local numbering for subsetted mesh.
			//  Gets new-to-old maps. Not compatible with unrefinement.
			void subset
			(
				const labelList& pointMap,
				const labelList& faceMap,
				const labelList& cellMap
			);

			//- Update local numbering for mesh redistribution
			void distribute(const mapDistributePolyMesh&);

			//- Debug: Check coupled mesh for correctness
			void checkMesh() const;

			//- Debug: Check 2:1 consistency across faces.
			//  maxPointDiff==-1 : only check 2:1 across faces
			//  maxPointDiff!=-1 : check point-connected cells.
			void checkRefinementLevels
			(
				const label maxPointDiff,
				const labelList& pointsToCheck
			) const;

		// Unrefinement (undoing refinement, not arbitrary coarsening)

			//- Return the points at the centre of top-level split cells
			//  that can be unsplit.
			labelList getSplitPoints() const;

			//- Given proposed
			//  splitPoints to unrefine according to calculate any clashes
			//  (due to 2:1) and return ok list of points to unrefine.
			//  Either adds points to refine to set (maxSet = true) or
			//  removes points to refine (maxSet = false)
			labelList consistentUnrefinement
			(
				const labelList& pointsToUnrefine,
				const bool maxSet
			) const;

			//- Remove some refinement. Needs to be supplied output of
			//  consistentUnrefinement. Only call if undoable set.
			//  All 8 pointCells of a split point will be combined into
			//  the lowest numbered cell of those 8.
			void setUnrefinement
			(
				const labelList& splitPointLabels,
				directTopoChange&
			);

		// Write

			// Set instance for mesh files
			void setInstance(const fileName& inst);

			//- Force writing refinement+history to polyMesh directory.
			bool write() const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
