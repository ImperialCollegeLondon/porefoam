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
	Foam::refinement

Description
	Abstract base class for adaptive mesh refinement using the mesh modifier
	engine. The class provides common interface and functionalities for 3D
	polyhedral refinement and 2D prismatic refinement.

	The common interface includes (pure virtuals) following member functions:
	- setCellsToRefine
	- setSplitPointsToUnrefine

	With a lot of ordinary protected member functions which are used by both
	derived classes.

	Note: I've written it this way in order to avoid unnecesasry code
	duplication, but I'm 99% sure that if someone else wants to write additional
	refinement strategy (e.g. directional refinement) as derived class, the
	interface will need to change.

SourceFiles
	refinement.C

Author
	Vuko Vukcevic, Wikki Ltd.  All rights reserved.

Notes
	Generalisation of hexRef8 for polyhedral cells and refactorisation using
	polyMesh modifier engine.

\*---------------------------------------------------------------------------*/

#ifndef refinement_H
#define refinement_H

#include "polyMeshModifier.H"
#include "labelIOField.H"
#include "removeFaces.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class refinement
:
	public polyMeshModifier
{
protected:

	// Protected enumeration for refinement status
	enum refinementStatus
	{
		UNREFINED = -1,
		UNCHANGED = 0,
		REFINED = 1
	};


	// Protected data

		//- Reference to polyMesh for easy access in helper functions
		const polyMesh& mesh_;


		// Refinement control and handling

			//- List of cells to refine in this time step
			mutable labelList cellsToRefine_;

			//- List of split point labels to unrefine in this time step
			mutable labelList splitPointsToUnrefine_;

			//- Cell refinement level
			mutable labelIOField cellLevel_;

			//- Point refinement level
			mutable labelIOField pointLevel_;

			//- Helper list for original (old) cells that will be refined or
			//  unrefined. The list is updated in setPolyhedralRefinement and
			//  setPolyhedralUnfereinement and is used in updateMesh to update
			//  the cellLevel on the new mesh.
			//  Values stored in the list are:
			//  a) UNREFINED = -1 = cell is unrefined
			//  b) UNCHANGED =  0 = cell is untouched
			//  c) REFINED   = +1 = cell is refined
			mutable labelList refinementLevelIndicator_;

			//- Face remover engine
			mutable removeFaces faceRemover_;

			//- Maximum number of cells in the mesh. Note: not strictly enforced
			label maxCells_;

			//- Maximum number of refinement levels for a given cell
			label maxRefinementLevel_;

			//- Switch whether to use edge based consistency on refinement
			Switch edgeBasedConsistency_;

			//- Number of buffer layers for refinement
			label nRefinementBufferLayers_;

			//- Number of buffer layers for unrefinement, controlling how far
			//  the unrefinement region needs to be from current refinement
			//  region
			label nUnrefinementBufferLayers_;


	// Protected Pure Virtual Member Functions

		// Global topology modification functions (operate on whole polyMesh)

			//- Set refinement instruction
			virtual void setRefinementInstruction
			(
				polyTopoChange& ref
			) const = 0;

			//- Set unrefinement instruction
			virtual void setUnrefinementInstruction
			(
				polyTopoChange& ref
			) const = 0;


	// Protected Member Functions

		// Useful helper functions used by derived classes

			//- Set file instance for cellLevel_ and pointLevel_
			void setInstance(const fileName& inst) const;


		// Local topology modification functions (operate on cells/faces)

			//- Adds a face on top of existing faceI. Reverses if nessecary
			label addFace
			(
				polyTopoChange& ref,
				const label faceI,
				const face& newFace,
				const label own,
				const label nei
			) const;

			//- Adds internal face from point. No checks on reversal
			label addInternalFace
			(
				polyTopoChange& ref,
				const label meshFaceI,
				const label meshPointI,
				const face& newFace,
				const label own,
				const label nei
			) const;

			//- Modifies existing faceI for either new owner/neighbour or new
			//  face points. Reverses if nessecary
			void modifyFace
			(
				polyTopoChange& ref,
				const label faceI,
				const face& newFace,
				const label own,
				const label nei
			) const;


		// Topological change helper functions

			//- Store vertices from startFp up to face split point.
			//  Used when splitting face into n faces where n is the number of
			//  points in a face (or number of edges)
			void walkFaceToMid
			(
				const labelList& edgeMidPoint,
				const label cLevel,
				const label faceI,
				const label startFp,
				dynamicLabelList& faceVerts
			) const;

			//- Same as walkFaceToMid but now walk back
			void walkFaceFromMid
			(
				const labelList& edgeMidPoint,
				const label cLevel,
				const label faceI,
				const label startFp,
				dynamicLabelList& faceVerts
			) const;

			//- Get index of point with minimum point level
			label findMinLevel(const labelList& f) const;

			//- Get index of point with maximum point level
			label findMaxLevel(const labelList& f) const;

			//- Count number of vertices <= anchorLevel for a given face
			label countAnchors
			(
				const labelList& f,
				const label anchorLevel
			) const;

			//- Adjust cell refinement level after topo change
			void adjustRefLevel
			(
				label& curNewCellLevel,
				const label oldCellI
			);

	
		// Debug functions

			//- Check orientation of added internal face
			void checkInternalOrientation
			(
				polyTopoChange& ref,
				const label cellI,
				const label faceI,
				const point& ownPt,
				const point& neiPt,
				const face& newFace
			) const;

			//- Check orientation of a new boundary face
			void checkBoundaryOrientation
			(
				polyTopoChange& ref,
				const label cellI,
				const label faceI,
				const point& ownPt,
				const point& boundaryPt,
				const face& newFace
			) const;


		// Refinement/unrefinement consistency checks

			//- Updates cellsToRefine such that a face consistent 2:1 refinement
			//  is obtained. Returns local number of cells changed
			label faceConsistentRefinement(boolList& cellsToRefine) const;

			//- Updates cellsToRefine such that an edge consistent 4:1 refinement
			//  is obtained. Returns local number of cells changed
			label edgeConsistentRefinement(boolList& cellsToRefine) const;

			//- Updates cellsToUnrefine such that a face consistent 2:1
			//  unrefinement is obtained. Returns local number of cells changed
			label faceConsistentUnrefinement(boolList& cellsToUnrefine) const;

			//- Updates cellsToUnrefine such that an edge consistent 4:1
			//  unrefinement is obtained. Returns local number of cells changed
			label edgeConsistentUnrefinement(boolList& cellsToUnrefine) const;


		// Copy control

			//- Disallow default bitwise copy construct
			refinement(const refinement&);

			//- Disallow default bitwise assignment
			void operator=(const refinement&);


public:

	//- Runtime type information
	TypeName("refinement");


	// Constructors

		//- Construct from dictionary
		refinement
		(
			const word& name,
			const dictionary& dict,
			const label index,
			const polyTopoChanger& mme
		);


	//- Destructor
	virtual ~refinement();


	// Member Functions

		// Access

			//- Return refinement cell level
			const labelIOField& cellLevel() const
			{
				return cellLevel_;
			}

			//- Return refinement point level
			const labelIOField& pointLevel() const
			{
				return pointLevel_;
			}


		// Edit

			//- Set cells to refine given a list of refinement
			//  candidates. Refinement candidates are extended within the
			//  function due to possible 4:1 conflicts and specified number of
			//  buffer layers.
			//  Note: must be called BEFORE setSplitPointsToUnrefine
			virtual void setCellsToRefine
			(
				const labelList& refinementCellCandidates
			) = 0;

			//- Set split points to unrefine given a list of all mesh points
			//  that are candidates for unrefinement. Split points are
			//  determined as a subset of unrefinement candidates, avoiding
			//  splitting points of cells that are going to be refined at the
			//  same time and ensuring consistent unrefinement.
			//  Note: must be called AFTER setCellsToRefine
			virtual void setSplitPointsToUnrefine
			(
				const labelList& unrefinementPointCandidates
			) = 0;


		// Inherited interface from polyMeshModifier

			//- Check for topology change
			virtual bool changeTopology() const;

			//- Insert the polyhedral refinement/unrefinement into the
			//  topological change
			virtual void setRefinement(polyTopoChange&) const;

			//- Modify motion points to comply with the topological change
			virtual void modifyMotionPoints(pointField& motionPoints) const;

			//- Force recalculation of locally stored data on topological change
			virtual void updateMesh(const mapPolyMesh&);

			//- Write
			virtual void write(Ostream&) const;

			//- Write dictionary
			virtual void writeDict(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
