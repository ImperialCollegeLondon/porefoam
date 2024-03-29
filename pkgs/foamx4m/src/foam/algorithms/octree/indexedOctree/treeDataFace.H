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
	Foam::treeDataFace

Description
	Encapsulation of data needed to search for faces.

SourceFiles
	treeDataFace.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataFace_H
#define treeDataFace_H

#include "face.H"
#include "indexedOctree.H"
#include "treeBoundBoxList.H"
#include "PackedBoolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
//template<class Type> class indexedOctree;
class polyPatch;


class treeDataFace
{
	// Static data

		//- tolerance on linear dimensions
		static scalar tolSqr;



	// Private data

		const primitiveMesh& mesh_;

		//- Subset of faces to work on
		const labelList faceLabels_;

		//- Inverse of faceLabels. For every mesh whether face is in faceLabels.
		PackedBoolList isTreeFace_;

		//- Whether to precalculate and store face bounding box
		const bool cacheBb_;

		//- face bounding boxes (valid only if cacheBb_)
		treeBoundBoxList bbs_;


	// Private Member Functions

		//- Calculate face bounding box
		treeBoundBox calcBb(const label cellI) const;

		//- Initialise all member data
		void update();

public:

	// Declare name of the class and its debug switch
	ClassName("treeDataFace");


	// Constructors

		//- Construct from mesh and subset of faces.
		treeDataFace
		(
			const bool cacheBb,
			const primitiveMesh&,
			const labelList&
		);

		//- Construct from mesh. Uses all faces in mesh.
		treeDataFace(const bool cacheBb, const primitiveMesh&);

		//- Construct from mesh. Uses all faces in patch.
		treeDataFace(const bool cacheBb, const polyPatch&);


	// Member Functions

		// Access

			const labelList& faceLabels() const
			{
				return faceLabels_;
			}

			const primitiveMesh& mesh() const
			{
				return mesh_;
			}

			label size() const
			{
				return faceLabels_.size();
			}

			//- Get representative point cloud for all shapes inside
			//  (one point per shape)
			pointField points() const;


		// Search

			//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
			//  Only makes sense for closed surfaces.
			label getVolumeType
			(
				const indexedOctree<treeDataFace>&,
				const point&
			) const;

			//- Does (bb of) shape at index overlap bb
			bool overlaps
			(
				const label index,
				const treeBoundBox& sampleBb
			) const;

			//- Calculates nearest (to sample) point in shape.
			//  Returns actual point and distance (squared)
			void findNearest
			(
				const labelList& indices,
				const point& sample,

				scalar& nearestDistSqr,
				label& nearestIndex,
				point& nearestPoint
			) const;

			//- Calculates nearest (to line) point in shape.
			//  Returns point and distance (squared)
			void findNearest
			(
				const labelList& indices,
				const linePointRef& ln,

				treeBoundBox& tightest,
				label& minIndex,
				point& linePoint,
				point& nearestPoint
			) const
			{
				notImplemented
				(
					"treeDataFace::findNearest"
					"(const labelList&, const linePointRef&, ..)"
				);
			}

			//- Calculate intersection of shape with ray. Sets result
			//  accordingly
			bool intersects
			(
				const label index,
				const point& start,
				const point& end,
				point& result
			) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
