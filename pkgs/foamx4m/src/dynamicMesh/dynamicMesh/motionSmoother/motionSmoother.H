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
	Foam::motionSmoother

Description
	Given a displacement moves the mesh by scaling the displacement back
	until there are no more mesh errors.

	Holds displacement field (read upon construction since need boundary
	conditions) and scaling factor and optional patch number on which to
	scale back displacement.

	E.g.
	@verbatim
		// Construct iterative mesh mover.
		motionSmoother meshMover(mesh, labelList(1, patchI));

		// Set desired displacement:
		meshMover.displacement() = ..

		for (label iter = 0; iter < maxIter; iter++)
		{
			if (meshMover.scaleMesh(true))
			{
				Info<< "Successfully moved mesh" << endl;
				return true;
			}
		}
	@endverbatim

Note
	- Shared points (parallel): a processor can have points which are part of
	pp on another processor but have no pp itself (i.e. it has points
	and/or edges but no faces of pp). Hence we have to be careful when e.g.
	synchronising displacements that the value from the processor which has
	faces of pp get priority. This is currently handled in setDisplacement
	by resetting the internal displacement to zero before doing anything
	else. The combine operator used will give preference to non-zero
	values.

	- Various routines take baffles. These are sets of boundary faces that
	are treated as a single internal face. This is a hack used to apply
	movement to internal faces.

	- Mesh constraints are looked up from the supplied dictionary. (uses
	recursive lookup)

SourceFiles
	motionSmoother.C
	motionSmootherTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef motionSmoother_H
#define motionSmoother_H

#include "objectRegistry.H"
#include "foamTime.H"
#include "pointFields.H"
#include "HashSet.H"
#include "PackedBoolList.H"
#include "indirectPrimitivePatch.H"
#include "className.H"
#include "twoDPointCorrector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class polyMeshGeometry;
class faceSet;


class motionSmoother
{
	// Private class

		//- To synchronise displacements. We want max displacement since
		//  this is what is specified on pp and internal mesh will have
		//  zero displacement.
		class maxMagEqOp
		{

		public:

			void operator()(vector& x, const vector& y) const
			{
				for (direction i = 0; i < vector::nComponents; i++)
				{
					scalar magX = mag(x[i]);
					scalar magY = mag(y[i]);

					if (magX < magY)
					{
					    x[i] = y[i];
					}
					else if (magX == magY)
					{
					    if (y[i] > x[i])
					    {
					        x[i] = y[i];
					    }
					}
				}
			}
		};


	// Private data

		//- Reference to polyMesh. Non-const since we move mesh.
		polyMesh& mesh_;

		//- Reference to pointMesh
		pointMesh& pMesh_;

		//- Reference to face subset of all adaptPatchIDs
		indirectPrimitivePatch& pp_;

		//- Indices of fixedValue patches that we're allowed to modify the
		// displacement on.
		const labelList adaptPatchIDs_;


		// Smoothing and checking parameters
		dictionary paramDict_;

		// Internal data

		//- Displacement field
		pointVectorField displacement_;

		//- Scale factor for displacement
		pointScalarField scale_;

		//- Starting mesh position
		pointField oldPoints_;

		//- Is mesh point on boundary or not
		PackedBoolList isInternalPoint_;

		//- Is edge master (always except if on coupled boundary and on
		//  lower processor)
		PackedBoolList isMasterEdge_;

		//- 2-D motion corrector
		twoDPointCorrector twoDCorrector_;

		// Muli-patch constraints (from pointPatchInterpolation)

			labelList patchPatchPointConstraintPoints_;
			tensorField patchPatchPointConstraintTensors_;


	// Private Member Functions

		//- Average of connected points.
		template <class Type>
		tmp<GeometricField<Type, pointPatchField, pointMesh> > avg
		(
			const GeometricField<Type, pointPatchField, pointMesh>& fld,
			const scalarField& edgeWeight,
			const bool separation
		) const;

		//- Check constraints
		template<class Type>
		static void checkConstraints
		(
			GeometricField<Type, pointPatchField, pointMesh>&
		);

		//- Multi-patch constraints
		template<class Type>
		void applyCornerConstraints
		(
			GeometricField<Type, pointPatchField, pointMesh>&
		) const;

		//- Test synchronisation of pointField
		template<class Type, class CombineOp>
		void testSyncField
		(
			const Field<Type>&,
			const CombineOp& cop,
			const Type& zero,
			const bool separation,
			const scalar maxMag
		) const;

		//- Assemble tensors for multi-patch constraints
		void makePatchPatchAddressing();

		static void checkFld(const pointScalarField&);

		//- Get points used by given faces
		labelHashSet getPoints(const labelHashSet&) const;

		//- explicit smoothing and min on all affected internal points
		void minSmooth
		(
			const PackedBoolList& isAffectedPoint,
			const pointScalarField& fld,
			pointScalarField& newFld
		) const;

		//- same but only on selected points (usually patch points)
		void minSmooth
		(
			const PackedBoolList& isAffectedPoint,
			const labelList& meshPoints,
			const pointScalarField& fld,
			pointScalarField& newFld
		) const;

		//- Scale certain (internal) points of a field
		void scaleField
		(
			const labelHashSet& pointLabels,
			const scalar scale,
			pointScalarField&
		) const;

		//- As above but points have to be in meshPoints as well
		//  (usually to scale patch points)
		void scaleField
		(
			const labelList& meshPoints,
			const labelHashSet& pointLabels,
			const scalar scale,
			pointScalarField&
		) const;

		//- Helper function. Is point internal?
		bool isInternalPoint(const label pointI) const;

		//- Given a set of faces that cause smoothing and a number of
		//  iterations determine the maximum set of points who are affected
		//  and the accordingly affected faces.
		void getAffectedFacesAndPoints
		(
			const label nPointIter,
			const faceSet& wrongFaces,

			labelList& affectedFaces,
			PackedBoolList& isAffectedPoint
		) const;

		//- Disallow default bitwise copy construct
		motionSmoother(const motionSmoother&);

		//- Disallow default bitwise assignment
		void operator=(const motionSmoother&);


public:

	ClassName("motionSmoother");

	// Constructors

		//- Construct from mesh, patches to work on and smoothing parameters.
		//  Reads displacement field (only boundary conditions used)
		motionSmoother
		(
			polyMesh&,
			pointMesh&,
			indirectPrimitivePatch& pp,         // 'outside' points
			const labelList& adaptPatchIDs,     // patches forming 'outside'
			const dictionary& paramDict
		);

		//- Construct from mesh, patches to work on and smoothing parameters
		//   and displacement field (only boundary conditions used)
		motionSmoother
		(
			polyMesh&,
			indirectPrimitivePatch& pp,         // 'outside' points
			const labelList& adaptPatchIDs,     // patches forming 'outside'
			const pointVectorField&,
			const dictionary& paramDict
		);


	// Destructor

		~motionSmoother();


	// Member Functions

		// Access

			//- Reference to mesh
			const polyMesh& mesh() const;

			//- Reference to pointMesh
			const pointMesh& pMesh() const;

			//- Reference to patch
			const indirectPrimitivePatch& patch() const;

			//- Patch labels that are being adapted
			const labelList& adaptPatchIDs() const;

			const dictionary& paramDict() const;

			//- Reference to displacement field
			pointVectorField& displacement();

			//- Reference to displacement field
			const pointVectorField& displacement() const;

			//- Reference to scale field
			const pointScalarField& scale() const;

			//- Starting mesh position
			const pointField& oldPoints() const;

			//- Return reference to 2D point motion correction
			twoDPointCorrector& twoDCorrector()
			{
				return twoDCorrector_;
			}



		// Edit

			//- Take over existing mesh position.
			void correct();

			//- Set displacement field from displacement on patch points.
			//  Modify provided displacement to be consistent with actual
			//  boundary conditions on displacement. Note: resets the
			//  displacement to be 0 on coupled patches beforehand
			//  to make sure shared points
			//  partially on pp (on some processors) and partially not
			//  (on other processors) get the value from pp.
			void setDisplacement(pointField& patchDisp);

			//- Special correctBoundaryConditions which evaluates fixedValue
			//  patches first so they get overwritten with any constraint
			//  bc's.
			void correctBoundaryConditions(pointVectorField&) const;

			//- Move mesh. Does 2D correction (modifies passed pointField) and
			//  polyMesh::movePoints. Returns swept volumes.
			tmp<scalarField> movePoints(pointField&);

			//- Set the errorReduction (by how much to scale the displacement
			//  at error locations) parameter. Returns the old value.
			//  Set to 0 (so revert to old mesh) grows out one cell layer
			//  from error faces.
			scalar setErrorReduction(const scalar);

			//- Move mesh with given scale. Return true if mesh ok or has
			//  less than nAllow errors, false
			//  otherwise and locally update scale. Smoothmesh=false means only
			//  patch points get moved.
			//  Parallel ok (as long as displacement field is consistent
			//  across patches)
			bool scaleMesh
			(
				labelList& checkFaces,
				const bool smoothMesh = true,
				const label nAllow = 0
			);

			//- Move mesh (with baffles) with given scale.
			bool scaleMesh
			(
				labelList& checkFaces,
				const List<labelPair>& baffles,
				const bool smoothMesh = true,
				const label nAllow = 0
			);

			//- Move mesh with externally provided mesh constraints
			bool scaleMesh
			(
				labelList& checkFaces,
				const List<labelPair>& baffles,
				const dictionary& paramDict,
				const dictionary& meshQualityDict,
				const bool smoothMesh = true,
				const label nAllow = 0
			);

			//- Update topology
			void updateMesh();

			//- Check mesh with mesh settings in dict. Collects incorrect faces
			//  in set. Returns true if one or more faces in error.
			//  Parallel ok.
			static bool checkMesh
			(
				const bool report,
				const polyMesh& mesh,
				const dictionary& dict,
				labelHashSet& wrongFaces
			);

			//- Check (subset of mesh) with mesh settings in dict.
			//  Collects incorrect faces in set. Returns true if one
			//  or more faces in error. Parallel ok.
			static bool checkMesh
			(
				const bool report,
				const polyMesh& mesh,
				const dictionary& dict,
				const labelList& checkFaces,
				labelHashSet& wrongFaces
			);

			//- Check (subset of mesh including baffles) with mesh settings
			//  in dict. Collects incorrect faces in set. Returns true if one
			//  or more faces in error. Parallel ok.
			static bool checkMesh
			(
				const bool report,
				const polyMesh& mesh,
				const dictionary& dict,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet& wrongFaces
			);

			//- Check part of mesh with mesh settings in dict.
			//  Collects incorrect faces in set. Returns true if one or
			//  more faces in error. Parallel ok.
			static bool checkMesh
			(
				const bool report,
				const dictionary& dict,
				const polyMeshGeometry&,
				const labelList& checkFaces,
				labelHashSet& wrongFaces
			);

			//- Check part of mesh including baffles with mesh settings in dict.
			//  Collects incorrect faces in set. Returns true if one or
			//  more faces in error. Parallel ok.
			static bool checkMesh
			(
				const bool report,
				const dictionary& dict,
				const polyMeshGeometry&,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet& wrongFaces
			);

			// Helper functions to manipulate displacement vector.

				//- Fully explicit smoothing of internal points with varying
				//  diffusivity.
				template <class Type>
				void smooth
				(
					const GeometricField<Type, pointPatchField, pointMesh>& fld,
					const scalarField& edgeWeight,
					const bool separation,
					GeometricField<Type, pointPatchField, pointMesh>& newFld
				) const;
};


template<>
void motionSmoother::applyCornerConstraints<scalar>
(
	GeometricField<scalar, pointPatchField, pointMesh>& pf
) const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "motionSmootherTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
