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
	Foam::FaceCellWave

Description
	Wave propagation of information through grid. Every iteration
	information goes through one layer of cells. Templated on information
	that is transferred.

	Handles parallel and cyclics and non-parallel cyclics.

	Note: whether to propagate depends on the return value of Type::update
	which returns true (i.e. propagate) if the value changes by more than a
	certain tolerance.
	This tolerance can be very strict for normal face-cell and parallel
	cyclics (we use a value of 0.01 just to limit propagation of small changes)
	but for non-parallel cyclics this tolerance can be critical and if chosen
	too small can lead to non-convergence.

SourceFiles
	FaceCellWave.C

\*---------------------------------------------------------------------------*/

#ifndef FaceCellWave_H
#define FaceCellWave_H

#include "boolList.H"
#include "labelList.H"
#include "primitiveFieldsFwd.H"
#include "tolerancesSwitch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class polyPatch;


TemplateName(FaceCellWave);



template <class Type>
class FaceCellWave
:
	public FaceCellWaveName
{
	// Private data

		//- Reference to mesh
		const polyMesh& mesh_;

		//- Information for all faces
		UList<Type>& allFaceInfo_;

		//- Information for all cells
		UList<Type>& allCellInfo_;

		//- Has face changed
		boolList changedFace_;

		//- List of changed faces
		labelList changedFaces_;

		//- Number of changed faces
		label nChangedFaces_;

		// Cells that have changed
		boolList changedCell_;
		labelList changedCells_;
		label nChangedCells_;

		//- Contains cyclics
		bool hasCyclicPatches_;

		//- Number of evaluations
		label nEvals_;

		//- Number of unvisited cells/faces
		label nUnvisitedCells_;
		label nUnvisitedFaces_;

		//- Iteration counter
		label iter_;


	// Static Functions

		//- Write faces info
		static Ostream& writeFaces
		(
			const label nFaces,
			const labelList& faceLabels,
			const List<Type>& faceInfo,
			Ostream& os
		);

		//- Read faces info
		static Istream& readFaces
		(
			label& nFaces,
			labelList& faceLabels,
			List<Type>& faceInfo,
			Istream& is
		);


	// Private Member Functions

		//- Disallow default bitwise copy construct
		FaceCellWave(const FaceCellWave&);

		//- Disallow default bitwise assignment
		void operator=(const FaceCellWave&);


		//- Updates cellInfo with information from neighbour. Updates all
		//  statistics.
		bool updateCell
		(
			const label cellI,
			const label neighbourFaceI,
			const Type& neighbourInfo,
			const scalar tol,
			Type& cellInfo
		);

		//- Updates faceInfo with information from neighbour. Updates all
		//  statistics.
		bool updateFace
		(
			const label faceI,
			const label neighbourCellI,
			const Type& neighbourInfo,
			const scalar tol,
			Type& faceInfo
		);

		//- Updates faceInfo with information from same face. Updates all
		//  statistics.
		bool updateFace
		(
			const label faceI,
			const Type& neighbourInfo,
			const scalar tol,
			Type& faceInfo
		);


		// Parallel, cyclic

			//- Debugging: check info on both sides of cyclic
			void checkCyclic(const polyPatch& pPatch) const;

			//- Has patches of certain type?
			bool hasPatchType(const word& nameOfType);

			//- Merge received patch data into global data
			void mergeFaceInfo
			(
				const polyPatch& patch,
				const label nFaces,
				const labelList&,
				const List<Type>&,
				const bool isParallel
			);

			//- Extract info for single patch only
			label getChangedPatchFaces
			(
				const polyPatch& patch,
				const label startFaceI,
				const label nFaces,
				labelList& changedPatchFaces,
				List<Type>& changedPatchFacesInfo
			) const;

			//- Handle leaving domain. Implementation referred to Type
			void leaveDomain
			(
				const polyPatch& patch,
				const label nFaces,
				const labelList& faceLabels,
				List<Type>& faceInfo
			) const;

			//- Handle leaving domain. Implementation referred to Type
			void enterDomain
			(
				const polyPatch& patch,
				const label nFaces,
				const labelList& faceLabels,
				List<Type>& faceInfo
			) const;

			//- Send info to neighbour
			void sendPatchInfo
			(
				const label neighbour,
				const label nFaces,
				const labelList&,
				const List<Type>&
			) const;

			//- Receive info from neighbour. Returns number of faces received.
			label receivePatchInfo
			(
				const label neighbour,
				labelList&,
				List<Type>&
			) const;

			//- Offset face labels by constant value
			static void offset
			(
				const polyPatch& patch,
				const label off,
				const label nFaces,
				labelList& faces
			);

			//- Apply transformation to Type
			void transform
			(
				const tensorField& rotTensor,
				const label nFaces,
				List<Type>& faceInfo
			);

			//- Merge data from across processor boundaries
			void handleProcPatches();

			//- Merge data from across cyclics
			void handleCyclicPatches();


	  // Private static data

			static debug::tolerancesSwitch geomTol_;
			static debug::tolerancesSwitch propagationTol_;

public:

	// Static Functions

		//- Access to tolerance
		static scalar propagationTol()
		{
			return propagationTol_();
		}

		//- Change tolerance
		static void setPropagationTol(const scalar tol)
		{
			propagationTol_ = tol;
		}


	// Constructors

		// Construct from mesh. Use setFaceInfo and iterate() to do actual
		// calculation.
		FaceCellWave
		(
			const polyMesh&,
			UList<Type>& allFaceInfo,
			UList<Type>& allCellInfo
		);

		//- Construct from mesh and list of changed faces with the Type
		//  for these faces. Iterates until nothing changes or maxIter reached.
		//  (maxIter can be 0)
		FaceCellWave
		(
			const polyMesh&,
			const labelList& initialChangedFaces,
			const List<Type>& changedFacesInfo,
			UList<Type>& allFaceInfo,
			UList<Type>& allCellInfo,
			const label maxIter
		);


	// Member Functions

		// Access

			//- Access allFaceInfo
			UList<Type>& allFaceInfo()
			{
				return allFaceInfo_;
			}

			//- Access allCellInfo
			UList<Type>& allCellInfo()
			{
				return allCellInfo_;
			}

			//- Access mesh
			const polyMesh& mesh() const
			{
				return mesh_;
			}

			//- Get number of unvisited cells, i.e. cells that were not (yet)
			//  reached from walking across mesh. This can happen from
			//  - not enough iterations done
			//  - a disconnected mesh
			//  - a mesh without walls in it
			label getUnsetCells() const;

			//- Get number of unvisited faces
			label getUnsetFaces() const;


		// Edit

			//- Set initial changed faces
			void setFaceInfo
			(
				const labelList& changedFaces,
				const List<Type>& changedFacesInfo
			);

			//- Propagate from face to cell. Returns total number of cells
			//  (over all processors) changed.
			label faceToCell();

			//- Propagate from cell to face. Returns total number of faces
			//  (over all processors) changed. (Faces on processorpatches are
			//  counted double)
			label cellToFace();

			//- Iterate until no changes or maxIter reached. Returns number of
			//  unset cells (see getUnsetCells)
			label iterate(const label maxIter);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "FaceCellWave.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
