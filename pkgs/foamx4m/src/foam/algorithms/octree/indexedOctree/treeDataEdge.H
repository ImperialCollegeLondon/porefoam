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
	Foam::treeDataEdge

Description
	Holds data for octree to work on an edges subset.

SourceFiles
	treeDataEdge.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataEdge_H
#define treeDataEdge_H

#include "treeBoundBoxList.H"
#include "linePointRef.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class indexedOctree;


class treeDataEdge
{
	// Static data

		//- tolerance on linear dimensions
		static scalar tol;


	// Private data

		//- Reference to edgeList
		const edgeList& edges_;

		//- Reference to points
		const pointField& points_;

		//- labels of edges
		const labelList edgeLabels_;

		//- Whether to precalculate and store face bounding box
		const bool cacheBb_;

		//- bbs for all above edges (valid only if cacheBb_)
		treeBoundBoxList bbs_;


	// Private Member Functions

		//- Calculate edge bounding box
		treeBoundBox calcBb(const label edgeI) const;

public:

	// Declare name of the class and its debug switch
	ClassName("treeDataEdge");


	// Constructors

		//- Construct from selected edges. !Holds references to edges and points
		treeDataEdge
		(
			const bool cacheBb,
			const edgeList& edges,
			const pointField& points,
			const labelList& edgeLabels
		);


	// Member Functions

		// Access

			label size() const
			{
				return edgeLabels_.size();
			}

			//- Get representative point cloud for all shapes inside
			//  (one point per shape)
			pointField points() const;


		// Search

			//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface.
			//  Only makes sense for closed surfaces.
			label getVolumeType
			(
				const indexedOctree<treeDataEdge>&,
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
			) const;

			//- Calculate intersection of shape with ray. Sets result
			//  accordingly
			bool intersects
			(
				const label index,
				const point& start,
				const point& end,
				point& result
			) const
			{
				notImplemented
				(
					"treeDataEdge::intersects(const label, const point&,"
					"const point&, point&)"
				);
				return false;
			}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif

// ************************************************************************* //
