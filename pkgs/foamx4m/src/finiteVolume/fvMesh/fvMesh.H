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
	Foam::fvMesh

Description
	Mesh data needed to do the Finite Volume discretisation.

	NOTE ON USAGE:
	fvMesh contains all the topological and geometric information
	related to the mesh.  It is also responsible for keeping the data
	up-to-date.  This is done by deleting the cell volume, face area,
	cell/face centre, addressing and other derived information as
	required and recalculating it as necessary.  The fvMesh therefore
	reserves the right to delete the derived information upon every
	topological (mesh refinement/morphing) or geometric change (mesh
	motion).  It is therefore unsafe to keep local references to the
	derived data outside of the time loop.

SourceFiles
	fvMesh.C
	fvMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef fvMesh_H
#define fvMesh_H

#include "polyMesh.H"
#include "lduMesh.H"
#include "primitiveMesh.H"
#include "fvBoundaryMesh.H"
#include "surfaceInterpolation.H"
#include "DimensionedField.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "pointFieldsFwd.H"
#include "slicedVolFieldsFwd.H"
#include "slicedSurfaceFieldsFwd.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMeshLduAddressing;
class volMesh;



class fvMesh
:
	public polyMesh,
	public lduMesh,
	public surfaceInterpolation
{
	// Private data

		//- Boundary mesh
		fvBoundaryMesh boundary_;


	// Demand-driven data

		mutable fvMeshLduAddressing* lduPtr_;

		//- Current time index for cell volumes
		//  Note.  The whole mechanism will be replaced once the
		//  dimensionedField is created and the dimensionedField
		//  will take care of the old-time levels.
		mutable label curTimeIndex_;

		//- Cell volumes
		mutable DimensionedField<scalar, volMesh>* VPtr_;

		//- Cell volumes old time level
		mutable DimensionedField<scalar, volMesh>* V0Ptr_;

		//- Cell volumes old-old time level
		mutable DimensionedField<scalar, volMesh>* V00Ptr_;

		//- Face area vectors
		mutable slicedSurfaceVectorField* SfPtr_;

		//- Mag face area vectors
		mutable surfaceScalarField* magSfPtr_;

		//- Cell centres
		mutable slicedVolVectorField* CPtr_;

		//- Face centres
		mutable slicedSurfaceVectorField* CfPtr_;

		//- Face motion fluxes
		mutable surfaceScalarField* phiPtr_;


	// Private Member Functions

	   // Make geometric data

			//- Make face areas
			void makeSf() const;

			//- Make face area magnitudes
			void makeMagSf() const;

			//- Make mesh motion fluxes
			void makePhi() const;

			//- Update mesh motion fluxes given swept volumes
			void updatePhi(const scalarField& sweptVols) const;

			//- Make cell centres
			void makeC() const;

			//- Make face centres
			void makeCf() const;


		//- Disallow construct as copy
		fvMesh(const fvMesh&);

		//- Disallow assignment
		void operator=(const fvMesh&);


		// Storage management

			//- Clear geometry but not the old-time cell volumes
			void clearGeomNotOldVol();

			//- Clear geometry
			void clearGeom();

			//- Clear addressing
			void clearAddressing();


public:

	// Public typedefs

		typedef fvMesh Mesh;
		typedef fvBoundaryMesh BoundaryMesh;


	// Declare name of the class and its debug switch
	TypeName("fvMesh");


	// Constructors

		//- Construct from IOobject
		explicit fvMesh(const IOobject& io);

		//- Construct from Istream
		fvMesh
		(
			const IOobject& io,
			Istream& is,
			const bool syncPar = true
		);

		//- Construct from components without boundary.
		//  Boundary is added using addFvPatches() member function
		fvMesh
		(
			const IOobject& io,
			const Xfer<pointField>& points,
			const Xfer<faceList>& faces,
			const Xfer<labelList>& allOwner,
			const Xfer<labelList>& allNeighbour,
			const bool syncPar = true
		);

		//- Construct without boundary from cells rather than owner/neighbour.
		//  Boundary is added using addPatches() member function
		fvMesh
		(
			const IOobject& io,
			const Xfer<pointField>& points,
			const Xfer<faceList>& faces,
			const Xfer<cellList>& cells,
			const bool syncPar = true
		);

		//- Construct from cell shapes
		fvMesh
		(
			const IOobject& io,
			const Xfer<pointField>& points,
			const cellShapeList& shapes,
			const faceListList& boundaryFaces,
			const wordList& boundaryPatchNames,
			const wordList& boundaryPatchTypes,
			const word& defaultBoundaryPatchName,
			const word& defaultBoundaryPatchType,
			const wordList& boundaryPatchPhysicalTypes,
			const bool syncPar = true
		);


	//- Destructor
	virtual ~fvMesh();


	// Member Functions

		// Helpers

			//- Add boundary patches. Constructor helper
			//  Note: adding polyPatches to polyMesh patches and
			//  fvPatches to fvMesh
			void addFvPatches
			(
				const List<polyPatch*>&,
				const bool validBoundary = true
			);

			//- Update the mesh based on the mesh files saved in time
			//  directories
			virtual readUpdateState readUpdate();


		// Access

			//- Return the top-level database
			const Time& time() const
			{
				return polyMesh::time();
			}

			//- Return the object registry - resolve conflict polyMesh/lduMesh
			virtual const objectRegistry& thisDb() const
			{
				return polyMesh::thisDb();
			}

			//- Return reference to name
			//  Note: name() is currently ambiguous due to derivation from
			//  surfaceInterpolation
			const word& name() const
			{
				return polyMesh::name();
			}

			//- Return reference to boundary mesh
			const fvBoundaryMesh& boundary() const;

			//- Return ldu addressing
			virtual const lduAddressing& lduAddr() const;

			//- Return a list of pointers for each patch
			//  with only those pointing to interfaces being set
			virtual lduInterfacePtrsList interfaces() const
			{
				return boundary().interfaces();
			}

			//- Internal face owner
			const unallocLabelList& owner() const
			{
				return lduAddr().lowerAddr();
			}

			//- Internal face neighbour
			const unallocLabelList& neighbour() const
			{
				return lduAddr().upperAddr();
			}


			// Communication support

				//- Return communicator used for parallel communication
				virtual int comm() const
				{
					return polyMesh::comm();
				}


			//- Return cell volumes
			const DimensionedField<scalar, volMesh>& V() const;

			//- Return old-time cell volumes
			const DimensionedField<scalar, volMesh>& V0() const;

			//- Return old-old-time cell volumes
			const DimensionedField<scalar, volMesh>& V00() const;

			//- Return sub-cycle cell volumes
			tmp<DimensionedField<scalar, volMesh> > Vsc() const;

			//- Return sub-cycl old-time cell volumes
			tmp<DimensionedField<scalar, volMesh> > Vsc0() const;

			//- Return cell face area vectors
			const surfaceVectorField& Sf() const;

			//- Return cell face area magnitudes
			const surfaceScalarField& magSf() const;

			//- Return cell face motion fluxes
			const surfaceScalarField& phi() const;

			//- Return cell centres as volVectorField
			const volVectorField& C() const;

			//- Return face centres as surfaceVectorField
			const surfaceVectorField& Cf() const;


		// Edit

			//- Clear all geometry and addressing
			void clearOut();

			//- Update mesh corresponding to the given map
			virtual void updateMesh(const mapPolyMesh& mpm);

			//- Sync mesh update for topo change on other processors
			//  Locally, there is no topological change
			virtual void syncUpdateMesh();

			//- Move points, returns volumes swept by faces in motion
			virtual tmp<scalarField> movePoints(const pointField&);

			//- Map all fields in time using given map.
			virtual void mapFields(const mapPolyMesh& mpm) const;

			//- Map cell volumes in time using given map.
			virtual void mapOldVolumes(const mapPolyMesh& mpm);

			//- Remove boundary patches. Warning: fvPatchFields hold ref to
			//  these fvPatches.
			void removeFvBoundary();

			//- Reset mesh primitive data. Assumes all patch info correct
			//  (so does e.g. parallel communication). If not, use
			//  validBoundary = false
			//  (still assumes patchStarts[0] = nInternalFaces and last
			//  patch ends at nActiveFaces) and change patches with addPatches.
			//  resetFvPatchFlag indicates which fvPatches will be rebuilt
			//  if the flag is set to false, original patch is preserved
			void resetFvPrimitives
			(
				const Xfer<pointField>& points,
				const Xfer<faceList>& faces,
				const Xfer<labelList>& owner,
				const Xfer<labelList>& neighbour,
				const labelList& patchSizes,
				const labelList& patchStarts,
				const boolList& resetFvPatchFlag,
				const bool validBoundary = true
			);

			//- Return cell face motion fluxes
			surfaceScalarField& setPhi();

			//- Set old-time cell volumes
			DimensionedField<scalar, volMesh>& setV0();


		// Write

			//- Write the underlying polyMesh and other data
			virtual bool writeObjects
			(
				IOstream::streamFormat fmt,
				IOstream::versionNumber ver,
				IOstream::compressionType cmp
			) const;

			//- Write mesh using IO settings from time
			virtual bool write() const;

			//- Return solver settings
			const dictionary& solverDict(const word& name) const
			{
				return solutionDict().solverDict(name);
			}

	// Member Operators

		bool operator!=(const fvMesh&) const;
		bool operator==(const fvMesh&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "fvPatchFvMeshTemplates.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
