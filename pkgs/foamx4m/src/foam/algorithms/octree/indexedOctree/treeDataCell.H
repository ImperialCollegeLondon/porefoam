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
	Foam::treeDataCell

Description
	Encapsulation of data needed to search in/for cells. Used to find the
	cell containing a point (e.g. cell-cell mapping).

SourceFiles
	treeDataCell.C

\*---------------------------------------------------------------------------*/

#ifndef treeDataCell_H
#define treeDataCell_H

#include "treeBoundBoxList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class primitiveMesh;
template<class Type> class indexedOctree;


class treeDataCell
{
	// Private data

		const primitiveMesh& mesh_;

		//- Subset of cells to work on
		const labelList cellLabels_;

		//- Whether to precalculate and store cell bounding box
		const bool cacheBb_;

		//- cell bounding boxes (valid only if cacheBb_)
		treeBoundBoxList bbs_;


	// Private Member Functions

		//- Calculate cell bounding box
		treeBoundBox calcCellBb(const label cellI) const;

public:

	// Declare name of the class and its debug switch
	ClassName("treeDataCell");


	// Constructors

		//- Construct from mesh and subset of cells.
		treeDataCell
		(
			const bool cacheBb,
			const primitiveMesh&,
			const labelList&
		);

		//- Construct from mesh. Uses all cells in mesh.
		treeDataCell(const bool cacheBb, const primitiveMesh&);


	// Member Functions

		// Access

			const labelList& cellLabels() const
			{
				return cellLabels_;
			}

			const primitiveMesh& mesh() const
			{
				return mesh_;
			}


			label size() const
			{
				return cellLabels_.size();
			}

			//- Get representative point cloud for all shapes inside
			//  (one point per shape)
			pointField points() const;


		// Search

			//- Get type (inside,outside,mixed,unknown) of point w.r.t. surface
			//  Only makes sense for closed surfaces
			label getVolumeType
			(
				const indexedOctree<treeDataCell>&,
				const point&
			) const
			{
				notImplemented
				(
					"treeDataCell::getVolumeType"
					"(const indexedOctree<treeDataCell>&, const point&)"
				);
				return -1;
			}

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
					"treeDataCell::findNearest"
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
