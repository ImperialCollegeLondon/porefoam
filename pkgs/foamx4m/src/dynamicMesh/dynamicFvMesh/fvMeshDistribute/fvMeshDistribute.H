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
	Foam::fvMeshDistribute

Description
	Sends/receives parts of mesh+fvfields to neighbouring processors.
	Used in load balancing.

	Input is per local cell the processor it should move to. Moves meshes
	and volFields/surfaceFields and returns map which can be used to
	distribute other.

	Notes:
	- does not handle cyclics. Will probably handle separated proc patches.
	- if all cells move off processor also all its processor patches will
	  get deleted so comms might be screwed up (since e.g. globalMeshData
	  expects procPatches on all)
	- initial mesh has to have procPatches last and all normal patches common
	  to all processors and in the same order. This is checked.
	- faces are matched topologically but points on the faces are not. So
	  expect problems -on separated patches (cyclics?) -on zero sized processor
	  edges.

SourceFiles
	fvMeshDistribute.C
	fvMeshDistributeTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef fvMeshDistribute_H
#define fvMeshDistribute_H

#include "Field.H"
//#include "uLabel.H"
#include "fvMeshSubset.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class mapAddedPolyMesh;
class mapDistributePolyMesh;


class fvMeshDistribute
{
	// Private data

		//- Underlying fvMesh
		fvMesh& mesh_;

		//- Absolute merging tolerance (constructing meshes gets done using
		//  geometric matching)
		const scalar mergeTol_;


	// Private Member Functions

		//- Find indices with value
		static labelList select
		(
			const bool selectEqual,
			const labelList& values,
			const label value
		);

		//- Check all procs have same names and in exactly same order.
		static void checkEqualWordList(const string&, const wordList&);

		//- Merge wordlists over all processors
		static wordList mergeWordList(const wordList&);


		// Patch handling

			//- Find patch to put exposed faces into.
			label findNonEmptyPatch() const;

			//- Appends processorPolyPatch. Returns patchID.
			label addProcPatch(const word& patchName, const label nbrProc);

			//- Add patch field
			template<class GeoField>
			void addPatchFields(const word& patchFieldType);

			//- Deletes last patch.
			void deleteTrailingPatch();

			// Delete trailing patch fields
			template<class GeoField>
			void deleteTrailingPatchFields();

			//- Save boundary fields
			template <class T, class Mesh>
			void saveBoundaryFields
			(
				PtrList<FieldField<fvsPatchField, T> >& bflds
			) const;

			//- Map boundary fields
			template <class T, class Mesh>
			void mapBoundaryFields
			(
				const mapPolyMesh& map,
				const PtrList<FieldField<fvsPatchField, T> >& oldBflds
			);

			//- Init patch fields of certain type
			template<class GeoField, class PatchFieldType>
			void initPatchFields
			(
				const typename GeoField::value_type& initVal
			);

			//- Delete all processor patches. Move any processor faces into
			//  patchI.
			autoPtr<mapPolyMesh> deleteProcPatches(const label patchI);

			//- Repatch the mesh. This is only nessecary for the proc
			//  boundary faces. newPatchID is over all boundary faces: -1 or
			//  new patchID. constructFaceMap is being adapted for the
			//  possible new face position (since proc faces get automatically
			//  matched)
			autoPtr<mapPolyMesh> repatch
			(
				const labelList& newPatchID,
				labelListList& constructFaceMap
			);

			//- Merge any shared points that are geometrically shared. Needs
			//  parallel valid mesh - uses globalMeshData.
			//  constructPointMap is adapted for the new point labels.
			autoPtr<mapPolyMesh> mergeSharedPoints
			(
				labelListList& constructPointMap
			);

		// Coupling information

			//- Construct the local environment of all boundary faces.
			void getNeighbourData
			(
				const labelList& distribution,
				labelList& sourceFace,
				labelList& sourceProc,
				labelList& sourceNewProc
			) const;

			// Subset the neighbourCell/neighbourProc fields
			static void subsetBoundaryData
			(
				const fvMesh& mesh,
				const labelList& faceMap,
				const labelList& cellMap,

				const labelList& oldDistribution,
				const labelList& oldFaceOwner,
				const labelList& oldFaceNeighbour,
				const label oldInternalFaces,

				const labelList& sourceFace,
				const labelList& sourceProc,
				const labelList& sourceNewProc,

				labelList& subFace,
				labelList& subProc,
				labelList& subNewProc
			);

			//- Find cells on mesh whose faceID/procID match the neighbour
			//  cell/proc of domainMesh. Store the matching face.
			static void findCouples
			(
				const primitiveMesh&,
				const labelList& sourceFace,
				const labelList& sourceProc,

				const label domain,
				const primitiveMesh& domainMesh,
				const labelList& domainFace,
				const labelList& domainProc,

				labelList& masterCoupledFaces,
				labelList& slaveCoupledFaces
			);

			//- Map data on boundary faces to new mesh (resulting from adding
			//  two meshes)
			static labelList mapBoundaryData
			(
				const primitiveMesh& mesh,      // mesh after adding
				const mapAddedPolyMesh& map,
				const labelList& boundaryData0, // mesh before adding
				const label nInternalFaces1,
				const labelList& boundaryData1  // added mesh
			);


		// Other

			//- Remove cells. Add all exposed faces to patch oldInternalPatchI
			autoPtr<mapPolyMesh> doRemoveCells
			(
				const labelList& cellsToRemove,
				const label oldInternalPatchI
			);

			//- Add processor patches. Changes mesh and returns per neighbour
			//  proc the processor patchID.
			void addProcPatches
			(
				const labelList&, // processor that neighbour is on
				labelList& procPatchID
			);

			//- Get boundary faces to be repatched. Is -1 or new patchID
			static labelList getProcBoundaryPatch
			(
				const labelList& neighbourNewProc,// new processor per b. face
				const labelList& procPatchID      // patchID
			);

			//- Send mesh and coupling data.
			static void sendMesh
			(
				const label domain,
				const fvMesh& mesh,
				const wordList& pointZoneNames,
				const wordList& facesZoneNames,
				const wordList& cellZoneNames,
				const labelList& sourceFace,
				const labelList& sourceProc,
				const labelList& sourceNewProc,
				OSstream& toDomain
			);
			//- Send subset of fields
			template<class GeoField>
			static void sendFields
			(
				const label domain,
				const wordList& fieldNames,
				const fvMeshSubset&,
				OSstream& toNbr
			);

			//- Receive mesh. Opposite of sendMesh
			static autoPtr<fvMesh> receiveMesh
			(
				const label domain,
				const wordList& pointZoneNames,
				const wordList& facesZoneNames,
				const wordList& cellZoneNames,
				const Time& runTime,
				labelList& domainSourceFace,
				labelList& domainSourceProc,
				labelList& domainSourceNewProc,
				ISstream& fromNbr
			);

			//- Receive fields. Opposite of sendFields
			template<class GeoField>
			static void receiveFields
			(
				const label domain,
				const wordList& fieldNames,
				fvMesh&,
				PtrList<GeoField>&,
				const dictionary& fieldDicts
			);

			//- Do parallel exchange
			template <class Container, class T>
			static void exchange
			(
				const List<Container >& sendBufs,
				List<Container >& recvBufs,
				labelListList& sizes
			);

			//- Disallow default bitwise copy construct
			fvMeshDistribute(const fvMeshDistribute&);

			//- Disallow default bitwise assignment
			void operator=(const fvMeshDistribute&);

public:

	ClassName("fvMeshDistribute");


	// Constructors

		//- Construct from mesh and absolute merge tolerance
		fvMeshDistribute(fvMesh& mesh, const scalar mergeTol);


	// Member Functions

		//- Helper function: count cells per processor in wanted distribution
		static labelList countCells(const labelList&);

		//- Send cells to neighbours according to distribution
		//  (for every cell the new proc)
		autoPtr<mapDistributePolyMesh> distribute(const labelList& dist);

		// Debugging

			//- Print some info on coupling data
			static void printCoupleInfo
			(
				const primitiveMesh&,
				const labelList&,
				const labelList&,
				const labelList&
			);

			//- Print some field info
			template<class GeoField>
			static void printFieldInfo(const fvMesh&);

			//- Print some info on mesh.
			static void printMeshInfo(const fvMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fvMeshDistributeTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
