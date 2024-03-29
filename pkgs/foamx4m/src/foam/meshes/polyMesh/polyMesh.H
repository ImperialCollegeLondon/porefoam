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
	Foam::polyMesh

Description
	Mesh consisting of general polyhedral cells.

SourceFiles
	polyMesh.C
	polyMeshInitMesh.C
	polyMeshClear.C
	polyMeshFromShapeMesh.C
	polyMeshIO.C
	polyMeshUpdate.C

\*---------------------------------------------------------------------------*/

#ifndef polyMesh_H
#define polyMesh_H

#include "objectRegistry.H"
#include "primitiveMesh.H"
#include "pointField.H"
#include "faceList.H"
#include "cellList.H"
#include "cellShapeList.H"
#include "pointIOField.H"
#include "faceIOList.H"
#include "labelIOList.H"
#include "polyBoundaryMesh.H"
#include "boundBox.H"
#include "pointZoneMesh.H"
#include "faceZoneMesh.H"
#include "cellZoneMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class globalMeshData;
class mapPolyMesh;

class polyMesh;
Ostream& operator<<(Ostream&, const polyMesh&);


class polyMesh
:
	public objectRegistry,
	public primitiveMesh
{
public:

	// Public data types

		//- Enumeration defining the state of the mesh after a read update.
		//  Used for post-processing applications, where the mesh
		//  needs to update based on the files written in time
		//  directores
		enum readUpdateState
		{
			UNCHANGED,
			POINTS_MOVED,
			TOPO_CHANGE,
			TOPO_PATCH_CHANGE
		};


private:

	// Permanent data

		// Primitive mesh data

			//- All points
			pointIOField allPoints_;

			//- Active points
			pointField::subField points_;

			//- All faces
			faceIOList allFaces_;

			//- Active faces
			faceList::subList faces_;

			//- Face owner
			labelIOList owner_;

			//- Face neighbour
			labelIOList neighbour_;

			//- Is the mesh in a parallel sync state?
			Switch syncPar_;

			//- Have the primitives been cleared
			bool clearedPrimitives_;


			//- Boundary mesh
			mutable polyBoundaryMesh boundary_;

			//- Mesh bounding-box.
			//  Created from points on construction, updated with mesh motion
			boundBox bounds_;

			//- Vector of non-constrained directions in mesh
			//  Defined according to the presence of empty and wedge patches
			mutable Vector<label> geometricD_;

			//- Vector of valid directions in mesh
			//  Defined according to the presence of empty patches
			mutable Vector<label> solutionD_;

			//- Communicator used for parallel communication
			label comm_;


		// Zoning information

			//- Point zones
			pointZoneMesh pointZones_;

			//- Face zones
			faceZoneMesh faceZones_;

			//- Cell zones
			cellZoneMesh cellZones_;


		//- Parallel info
		mutable globalMeshData* globalMeshDataPtr_;


		// Mesh motion related data

			//- Is the mesh moving
			bool moving_;

			//- Is the mesh changing (moving and/or topology changing)
			bool changing_;

			//- Current time index for mesh motion
			mutable label curMotionTimeIndex_;

			//- Old points (for the last mesh motion)
			mutable pointField* oldAllPointsPtr_;

			//- Old points
			mutable pointField::subField* oldPointsPtr_;


	// Private member functions

		//- Disallow construct as copy
		polyMesh(const polyMesh&);

		//- Disallow default bitwise assignment
		void operator=(const polyMesh&);

		//- Initialise the polyMesh from the primitive data
		void initMesh();

		//- Initialise the polyMesh from the given set of cells
		void initMesh(cellList& c);

		//- Calculate the valid directions in the mesh from the boundaries
		void calcDirections() const;

		//- Calculate the cell shapes from the primitive
		//  polyhedral information
		void calcCellShapes() const;


		// Helper functions for constructor from cell shapes

			labelListList cellShapePointCells(const cellShapeList&) const;

			labelList facePatchFaceCells
			(
				const faceList& patchFaces,
				const labelListList& pointCells,
				const faceListList& cellsFaceShapes,
				const label patchID
			) const;

			void setTopology
			(
				const cellShapeList& cellsAsShapes,
				const faceListList& boundaryFaces,
				const wordList& boundaryPatchNames,
				labelList& patchSizes,
				labelList& patchStarts,
				label& defaultPatchStart,
				label& nFaces,
				cellList& cells
			);


public:

	// Public typedefs

		typedef polyMesh Mesh;
		typedef polyBoundaryMesh BoundaryMesh;


	//- Runtime type information
	TypeName("polyMesh");


	//- Static mesh data

		//- Return the default region name
		static word defaultRegion;

		//- Return the mesh sub-directory name (usually "polyMesh")
		static word meshSubDir;


	//- Static data to control empty and wedge directions

		//- Empty direction tolerance
		static const debug::tolerancesSwitch emptyDirTol_;

		//- Wedge direction tolerance
		static const debug::tolerancesSwitch wedgeDirTol_;


	// Constructors

		//- Construct from IOobject
		explicit polyMesh(const IOobject& io);

		//- Construct from Istream
		polyMesh
		(
			const IOobject& io,
			Istream& is,
			const bool syncPar = true
		);

		//- Construct from components without boundary.
		//  Boundary is added using addPatches() member function
		polyMesh
		(
			const IOobject& io,
			const Xfer<pointField>& points,
			const Xfer<faceList>& faces,
			const Xfer<labelList>& owner,
			const Xfer<labelList>& neighbour,
			const bool syncPar = true
		);

		//- Construct without boundary with cells rather than owner/neighbour.
		//  Boundary is added using addPatches() member function
		polyMesh
		(
			const IOobject& io,
			const Xfer<pointField>& points,
			const Xfer<faceList>& faces,
			const Xfer<cellList>& cells,
			const bool syncPar = true
		);

		//- Construct from cell shapes
		polyMesh
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

		//- Construct from cell shapes with patch information in dictionary
		//  format.
		polyMesh
		(
			const IOobject& io,
			const Xfer<pointField>& points,
			const cellShapeList& shapes,
			const faceListList& boundaryFaces,
			const wordList& boundaryPatchNames,
			const PtrList<dictionary>& boundaryDicts,
			const word& defaultBoundaryPatchName,
			const word& defaultBoundaryPatchType,
			const bool syncPar = true
		);


	//- Destructor
	virtual ~polyMesh();


	// Member Functions

		// Database

			//- Override the objectRegistry dbDir for a single-region case
			virtual const fileName& dbDir() const;

			//- Return the local mesh directory (dbDir()/meshSubDir)
			fileName meshDir() const;

			//- Return the current instance directory for points
			//  Used in the consruction of gemometric mesh data dependent
			//  on points
			const fileName& pointsInstance() const;

			//- Return the current instance directory for faces
			const fileName& facesInstance() const;

			//- Set the instance for mesh files
			void setInstance(const fileName&);


			//- Set motion write option
			void setMotionWriteOpt(IOobject::writeOption);

			//- Set topological write option
			void setTopoWriteOpt(IOobject::writeOption);


		// Access

			//- Return all points, including inactive ones
			const pointField& allPoints() const;

			//- Return all faces, including inactive ones
			const faceList& allFaces() const;

			//- Return old mesh motion points, including inactive ones
			const pointField& oldAllPoints() const;


		// Active mesh data

			//- Return raw points
			virtual const pointField& points() const;

			//- Return raw faces
			virtual const faceList& faces() const;

			//- Return face owner
			virtual const labelList& faceOwner() const;

			//- Return face neighbour
			virtual const labelList& faceNeighbour() const;

			//- Return old points for mesh motion
			virtual const pointField& oldPoints() const;

			//- Return boundary mesh
			const polyBoundaryMesh& boundaryMesh() const
			{
				return boundary_;
			}


		// Mesh state and form

			//- Return mesh bounding box
			const boundBox& bounds() const
			{
				return bounds_;
			}

			//- Return parallel sync state
			const Switch& syncPar() const
			{
				return syncPar_;
			}

			//- Return access to parallel sync state
			Switch& syncPar()
			{
				return syncPar_;
			}

			//- Return the vector of geometric directions in mesh.
			//  Defined according to the presence of empty and wedge patches.
			//  1 indicates unconstrained direction and -1 a constrained
			//  direction.
			const Vector<label>& geometricD() const;

			//- Return the number of valid geometric dimensions in the mesh
			label nGeometricD() const;

			//- Return the vector of solved-for directions in mesh.
			//  Differs from geometricD in that it includes for wedge cases
			//  the circumferential direction in case of swirl.
			//  1 indicates valid direction and -1 an invalid direction.
			const Vector<label>& solutionD() const;

			//- Return the number of valid solved-for dimensions in the mesh
			label nSolutionD() const;


			// Communication support

				//- Return communicator used for parallel communication
				label comm() const;

				//- Return communicator used for parallel communication
				label& comm();


			// Point, face and cell zones

				//- Return point zone mesh
				const pointZoneMesh& pointZones() const
				{
					return pointZones_;
				}

				//- Return face zone mesh
				const faceZoneMesh& faceZones() const
				{
					return faceZones_;
				}

				//- Return cell zone mesh
				const cellZoneMesh& cellZones() const
				{
					return cellZones_;
				}


			//- Return parallel info
			const globalMeshData& globalData() const;

			//- Return the object registry
			const objectRegistry& thisDb() const
			{
				return *this;
			}


		// Mesh motion

			//- Is mesh moving
			bool moving() const
			{
				return moving_;
			}

			//- Return current motion time index
			label curMotionTimeIndex() const
			{
				return curMotionTimeIndex_;
			}

			//- Set the mesh to be moving
			bool moving(const bool m)
			{
				bool m0 = moving_;
				moving_ = m;
				changing_ = changing_ || moving_;
				return m0;
			}

			//- Is mesh changing (topology changing and/or moving)
			bool changing() const
			{
				return changing_;
			}

			//- Set the mesh to be changing
			bool changing(const bool c)
			{
				bool c0 = changing_;
				changing_ = c;
				return c0;
			}

			//- Move points, returns volumes swept by faces in motion
			virtual tmp<scalarField> movePoints(const pointField&);

			//- Reset motion
			void resetMotion() const;

			//- Set old points without executing motion
			void setOldPoints(const pointField& setPoints);


		// Topological change

			//- Return non-const access to the pointZones
			pointZoneMesh& pointZones()
			{
				return pointZones_;
			}

			//- Return non-const access to the faceZones
			faceZoneMesh& faceZones()
			{
				return faceZones_;
			}

			//- Return non-const access to the cellZones
			cellZoneMesh& cellZones()
			{
				return cellZones_;
			}

			//- Add boundary patches
			void addPatches
			(
				const List<polyPatch*>&,
				const bool validBoundary = true
			);

			//- Add mesh zones
			void addZones
			(
				const List<pointZone*>& pz,
				const List<faceZone*>& fz,
				const List<cellZone*>& cz
			);

			//- Update the mesh based on the mesh files saved in
			//  time directories
			virtual readUpdateState readUpdate();

			//- Update the mesh corresponding to given map
			virtual void updateMesh(const mapPolyMesh& mpm);

			//- Sync mesh update for topo change on other processors
			//  Locally, there is no topological change
			virtual void syncUpdateMesh();

			//- Remove boundary patches
			void removeBoundary();

			//- Remove point, face and cell zones
			void removeZones();

			//- Reset mesh primitive data. Assumes all patch info correct
			//  (so does e.g. parallel communication). If not, use
			//  validBoundary = false
			//  (still assumes patchStarts[0] = nInternalFaces and last
			//  patch ends at nActiveFaces) and change patches with addPatches.
			void resetPrimitives
			(
				const Xfer<pointField>& points,
				const Xfer<faceList>& faces,
				const Xfer<labelList>& owner,
				const Xfer<labelList>& neighbour,
				const labelList& patchSizes,
				const labelList& patchStarts,
				const bool validBoundary = true
			);


		//  Storage management

			//- Clear geometry
			void clearGeom();

			//- Clear addressing
			void clearAddressing();

			//- Clear all geometry and addressing unnecessary for CFD
			void clearOut();

			void clearPoints();

			//- Clear primitive data (points, faces and cells)
			void clearPrimitives();

			//- Remove all files from mesh instance
			void removeFiles(const fileName& instanceDir) const;

			//- Remove all files from mesh instance()
			void removeFiles() const;


		// Helper functions

			//- Find cell enclosing this location (-1 if not in mesh)
			label findCell(const point&) const;


		// IOstream Operators

			friend Ostream& operator<<(Ostream&, const polyMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
