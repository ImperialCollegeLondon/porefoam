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
	Foam::polyMeshGeometry

Description
	Updateable mesh geometry and checking routines.

	- non-ortho done across coupled faces.
	- faceWeight (delta factors) done across coupled faces.

SourceFiles
	polyMeshGeometry.C

\*---------------------------------------------------------------------------*/

#ifndef polyMeshGeometry_H
#define polyMeshGeometry_H

#include "objectRegistry.H"
#include "foamTime.H"
#include "pointFields.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class polyMeshGeometry
{
	// Private Data

		//- Reference to polyMesh.
		const polyMesh& mesh_;

		//- Uptodate copy of face areas
		vectorField faceAreas_;

		//- Uptodate copy of face centres
		vectorField faceCentres_;

		//- Uptodate copy of cell centres
		vectorField cellCentres_;

		//- Uptodate copy of cell volumes
		scalarField cellVolumes_;


	// Private Member Functions

		//- Update face areas and centres on selected faces.
		void updateFaceCentresAndAreas
		(
			const pointField& p,
			const labelList& changedFaces
		);

		//- Update cell volumes and centres on selected cells. Requires
		//  cells and faces to be consistent set.
		void updateCellCentresAndVols
		(
			const labelList& changedCells,
			const labelList& changedFaces
		);

		//- Detect&report non-ortho error for single face.
		static scalar checkNonOrtho
		(
			const polyMesh& mesh,
			const bool report,
			const scalar severeNonorthogonalityThreshold,
			const label faceI,
			const vector& s,    // face area vector
			const vector& d,    // cc-cc vector

			label& severeNonOrth,
			label& errorNonOrth,
			labelHashSet* setPtr
		);

		//- Calculate skewness given two cell centres and one face centre.
		static scalar calcSkewness
		(
			const point& ownCc,
			const point& neiCc,
			const point& fc
		);


public:

	ClassName("polyMeshGeometry");


	// Constructors

		//- Construct from mesh
		polyMeshGeometry(const polyMesh&);


	// Member Functions

		// Access

			const polyMesh& mesh() const
			{
				return mesh_;
			}

			const vectorField& faceAreas() const
			{
				return faceAreas_;
			}
			const vectorField& faceCentres() const
			{
				return faceCentres_;
			}
			const vectorField& cellCentres() const
			{
				return cellCentres_;
			}
			const scalarField& cellVolumes() const
			{
				return cellVolumes_;
			}

		// Edit

			//- Take over properties from mesh
			void correct();

			//- Recalculate on selected faces. Recalculates cell properties
			//  on owner and neighbour of these cells.
			void correct
			(
				const pointField& p,
				const labelList& changedFaces
			);

			//- Helper function: get affected cells from faces
			static labelList affectedCells
			(
				const polyMesh&,
				const labelList& changedFaces
			);


		// Checking of selected faces with supplied geometry (mesh only used for
		// topology). Coupled aware (coupled faces treated as internal ones)

			//- See primitiveMesh
			static bool checkFaceDotProduct
			(
				const bool report,
				const scalar orthWarn,
				const polyMesh&,
				const vectorField& cellCentres,
				const vectorField& faceAreas,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			);

			//- See primitiveMesh
			static bool checkFacePyramids
			(
				const bool report,
				const scalar minPyrVol,
				const polyMesh&,
				const vectorField& cellCentres,
				const pointField& p,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet*
			);

			//- See primitiveMesh
			static bool checkFaceSkewness
			(
				const bool report,
				const scalar internalSkew,
				const scalar boundarySkew,
				const polyMesh& mesh,
				const vectorField& cellCentres,
				const vectorField& faceCentres,
				const vectorField& faceAreas,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			);

			//- Interpolation weights (0.5 for regular mesh)
			static bool checkFaceWeights
			(
				const bool report,
				const scalar warnWeight,
				const polyMesh& mesh,
				const vectorField& cellCentres,
				const vectorField& faceCentres,
				const vectorField& faceAreas,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			);

			//- Cell volume ratio of neighbouring cells (1 for regular mesh)
			static bool checkVolRatio
			(
				const bool report,
				const scalar warnRatio,
				const polyMesh& mesh,
				const scalarField& cellVolumes,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			);

			//- See primitiveMesh
			static bool checkFaceAngles
			(
				const bool report,
				const scalar maxDeg,
				const polyMesh& mesh,
				const vectorField& faceAreas,
				const pointField& p,
				const labelList& checkFaces,
				labelHashSet* setPtr
			);

			//- Triangle (from face-centre decomposition) normal v.s.
			//  average face normal
			static bool checkFaceTwist
			(
				const bool report,
				const scalar minTwist,
				const polyMesh&,
				const vectorField& cellCentres,
				const vectorField& faceAreas,
				const vectorField& faceCentres,
				const pointField& p,
				const labelList& checkFaces,
				labelHashSet* setPtr
			);

			//- Consecutive triangle (from face-centre decomposition) normals
			static bool checkTriangleTwist
			(
				const bool report,
				const scalar minTwist,
				const polyMesh&,
				const vectorField& faceAreas,
				const vectorField& faceCentres,
				const pointField& p,
				const labelList& checkFaces,
				labelHashSet* setPtr
			);

			//- Small faces
			static bool checkFaceArea
			(
				const bool report,
				const scalar minArea,
				const polyMesh&,
				const vectorField& faceAreas,
				const labelList& checkFaces,
				labelHashSet* setPtr
			);

			static bool checkCellDeterminant
			(
				const bool report,
				const scalar minDet,
				const polyMesh&,
				const vectorField& faceAreas,
				const labelList& checkFaces,
				const labelList& affectedCells,
				labelHashSet* setPtr
			);


		// Checking of selected faces with local geometry. Uses above static
		// functions. Parallel aware.

			bool checkFaceDotProduct
			(
				const bool report,
				const scalar orthWarn,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			) const;

			bool checkFacePyramids
			(
				const bool report,
				const scalar minPyrVol,
				const pointField& p,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			) const;

			bool checkFaceSkewness
			(
				const bool report,
				const scalar internalSkew,
				const scalar boundarySkew,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			) const;

			bool checkFaceWeights
			(
				const bool report,
				const scalar warnWeight,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			) const;

			bool checkVolRatio
			(
				const bool report,
				const scalar warnRatio,
				const labelList& checkFaces,
				const List<labelPair>& baffles,
				labelHashSet* setPtr
			) const;

			bool checkFaceAngles
			(
				const bool report,
				const scalar maxDeg,
				const pointField& p,
				const labelList& checkFaces,
				labelHashSet* setPtr
			) const;

			bool checkFaceTwist
			(
				const bool report,
				const scalar minTwist,
				const pointField& p,
				const labelList& checkFaces,
				labelHashSet* setPtr
			) const;

			bool checkTriangleTwist
			(
				const bool report,
				const scalar minTwist,
				const pointField& p,
				const labelList& checkFaces,
				labelHashSet* setPtr
			) const;

			bool checkFaceArea
			(
				const bool report,
				const scalar minArea,
				const labelList& checkFaces,
				labelHashSet* setPtr
			) const;

			bool checkCellDeterminant
			(
				const bool report,
				const scalar warnDet,
				const labelList& checkFaces,
				const labelList& affectedCells,
				labelHashSet* setPtr
			) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
