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
	Foam::removePoints

Description
	Removes selected points from mesh and updates faces using these
	points.

SourceFiles
	removePoints.C

\*---------------------------------------------------------------------------*/

#ifndef removePoints_H
#define removePoints_H

#include "typeInfo.H"
#include "boolList.H"
#include "pointField.H"
#include "faceList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class polyMesh;
class directTopoChange;
class mapPolyMesh;
class face;


class removePoints
{

	// Private classes

		//- Combine-reduce operator to combine data on faces. Takes care
		//  of reverse orientation on coupled face.
		template <class T, template<class> class CombineOp>
		class faceEqOp
		{

		public:

			void operator()(List<T>& x, const List<T>& y) const
			{
				if (y.size() > 0)
				{
					if (x.size() == 0)
					{
					    x = y;
					}
					else
					{
					    label j = 0;
					    forAll(x, i)
					    {
					        CombineOp<T>()(x[i], y[j]);
					        j = y.rcIndex(j);
					    }
					}
				}
			}
		};


	// Private data

		//- Reference to mesh
		const polyMesh& mesh_;

		//- Whether undoable
		const bool undoable_;

		//- If undoable: deleted points
		pointField savedPoints_;

		//- If undoable: per stored face the original mesh face label
		labelList savedFaceLabels_;

		//- If undoable: per stored face the vertices. Negative indices
		//  refer to deletedPoints_
		faceList savedFaces_;


	// Private Member Functions

		//- Change the vertices of the face whilst keeping everything else
		//  (patch, zone) the same.
		void modifyFace
		(
			const label faceI,
			const face&,
			directTopoChange&
		) const;


		//- Disallow default bitwise copy construct
		removePoints(const removePoints&);

		//- Disallow default bitwise assignment
		void operator=(const removePoints&);

public:

	//- Runtime type information
	ClassName("removePoints");


	// Constructors

		//- Construct from mesh
		removePoints(const polyMesh& mesh, const bool undoable = false);


	// Member Functions

			//- If undoable: affected face labels. Already restored faces
			//  will be -1.
			const labelList& savedFaceLabels() const
			{
				return savedFaceLabels_;
			}


		// Helper functions

			//- Mark in pointCanBeDeleted the points that can be deleted
			//  (parallel synchronised) and returns the global number of these
			//  points. (this number is the global number before synchronisation
			//  so might be off!)
			//  A point can be deleted if
			//  - it is not used by any edge.
			//  or
			//  - is not used by an internal edge
			//  - is used by only two boundary edges. (note that these two
			//    edges will always be boundary ones!)
			//  - these two edges are sufficiently in line (cos > minCos)
			//  - all processors agree that point can be deleted.
			label countPointUsage
			(
				const scalar minCos,
				boolList& pointCanBeDeleted
			) const;

		// Topology changes

			//- Play commands into directTopoChange to remove points. Gets
			//  boolList (output of countPointUsage) as input.
			//  Does no check for whether resulting face is legal.
			//  Since pointCanBeDeleted is synced all coupled faces should
			//  decide the same.
			void setRefinement(const boolList&, directTopoChange&);

			//- Force recalculation of locally stored data on topological change
			void updateMesh(const mapPolyMesh&);

			//- Given set of faces to restore calculates a consistent set of
			//  saved faces (indices into savedFaces_) and saved vertices
			//  (indices into savedPoints_) to restore. The undoFaces have to
			//  be synced.
			void getUnrefimentSet
			(
				const labelList& undoFaces,
				labelList& localFaces,
				labelList& localPoints
			) const;

			//- Restore selected faces and vertices.
			void setUnrefinement
			(
				const labelList& localFaces,
				const labelList& localPoints,
				directTopoChange&
			);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
