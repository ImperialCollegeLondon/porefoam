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
	Foam::indexedOctree

Description
	Non-pointer based hierarchical recursive searching

SourceFiles
	indexedOctree.C

\*---------------------------------------------------------------------------*/

#ifndef indexedOctree_H
#define indexedOctree_H

#include "treeBoundBox.H"
#include "pointIndexHit.H"
#include "FixedList.H"
#include "Ostream.H"
#include "HashSet.H"
#include "labelBits.H"
#include "PackedList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
template<class Type> class indexedOctree;
template<class Type> Ostream& operator<<(Ostream&, const indexedOctree<Type>&);
class Istream;



TemplateName(indexedOctree);



template <class Type>
class indexedOctree
:
	public indexedOctreeName
{
public:

	// Data types

		//- volume types
		enum volumeType
		{
			UNKNOWN = 0,
			MIXED = 1,
			INSIDE = 2,
			OUTSIDE = 3
		};


	// Tree node. Has up pointer and down pointers.

		class node
		{
		public:

			//- Bounding box of this node
			treeBoundBox bb_;

			//- Parent node (index into nodes_ of tree)
			label parent_;

			//- IDs of the 8 nodes on all sides of the mid point
			FixedList<labelBits, 8> subNodes_;

			friend Ostream& operator<< (Ostream& os, const node& n)
			{
				return os << n.bb_ << token::SPACE
					<< n.parent_ << token::SPACE << n.subNodes_;
			}

			friend Istream& operator>> (Istream& is, node& n)
			{
				return is >> n.bb_ >> n.parent_ >> n.subNodes_;
			}

			friend bool operator==(const node& a, const node& b)
			{
				return
					a.bb_ == b.bb_
				 && a.parent_ == b.parent_
				 && a.subNodes_ == b.subNodes_;
			}

			friend bool operator!=(const node& a, const node& b)
			{
				return !(a == b);
			}
		};


private:

	// Static data

		//- Relative peturbation tolerance. Determines when point is
		//  considered to be close to face/edge of bb of node.
		//  The tolerance is relative to the bounding box of the smallest
		//  node.
		static scalar perturbTol_;


	// Private data

		//- Underlying shapes for geometric queries.
		const Type shapes_;

		//- List of all nodes
		List<node> nodes_;

		//- List of all contents (referenced by those nodes that are contents)
		labelListList contents_;

		//- Per node per octant whether is fully inside/outside/mixed.
		mutable PackedList<2> nodeTypes_;

	// Private Member Functions

		//- Helper: does bb intersect a sphere around sample? Or is any
		//  corner point of bb closer than nearestDistSqr to sample.
		//  (bb is implicitly provided as parent bb + octant)
		static bool overlaps
		(
			const treeBoundBox& parentBb,
			const direction octant,
			const scalar nearestDistSqr,
			const point& sample
		);

		// Construction

			//- Split list of indices into 8 bins according to where they are
			//  in relation to mid.
			void divide
			(
				const labelList& indices,
				const treeBoundBox& bb,
				labelListList& result
			) const;

			//- Subdivide the contents node at position contentI. Appends to
			//  contents.
			node divide
			(
				const treeBoundBox& bb,
				DynamicList<labelList>& contents,
				const label contentI
			) const;

			//- Split any contents node with more than minSize elements.
			void splitNodes
			(
				const label minSize,
				DynamicList<node>& nodes,
				DynamicList<labelList>& contents
			) const;

			static label compactContents
			(
				DynamicList<node>& nodes,
				DynamicList<labelList>& contents,
				const label compactLevel,
				const label nodeI,
				const label level,
				List<labelList>& compactedContents,
				label& compactI
			);

			//- Determine inside/outside per node (mixed if cannot be
			//  determined). Only valid for closed shapes.
			volumeType calcVolumeType(const label nodeI) const;

			//- Search cached volume type.
			volumeType getVolumeType(const label nodeI, const point&) const;


		// Query

			//- Find nearest point to line.
			void findNearest
			(
				const label nodeI,
				const linePointRef& ln,

				treeBoundBox& tightest,
				label& nearestShapeI,
				point& linePoint,
				point& nearestPoint
			) const;

			//- Return bbox of octant
			treeBoundBox subBbox
			(
				const label parentNodeI,
				const direction octant
			) const;

			//- Helper: take a point on/close to face of bb and push it
			//  inside or outside of bb.
			static point pushPoint
			(
				const treeBoundBox&,
				const point&,
				const bool pushInside
			);

			//- Helper: take a point on face of bb and push it
			//  inside or outside of bb.
			static point pushPoint
			(
				const treeBoundBox&,
				const direction,
				const point&,
				const bool pushInside
			);

			//- Helper: take point on face(s) of bb and push it away from
			//  edges of face.
			static point pushPointIntoFace
			(
				const treeBoundBox& bb,
				const vector& dir,          // end-start
				const point& pt
			);

			////- Push point on multiple faces away from any corner/edge.
			//static void checkMultipleFaces
			//(
			//    const treeBoundBox& bb,
			//    const vector& dir,          // end-start
			//    pointIndexHit& faceHitInfo,
			//    direction& faceID
			//);

			//- Walk to parent of node+octant.
			bool walkToParent
			(
				const label nodeI,
				const direction octant,

				label& parentNodeI,
				label& parentOctant
			) const;

			//- Walk tree to neighbouring node. Return false if edge of tree
			//  hit.
			bool walkToNeighbour
			(
				const point& facePoint,
				const direction faceID,         // direction to walk in
				label& nodeI,
				direction& octant
			) const;

			//- Debug: return verbose the bounding box faces
			static word faceString(const direction faceID);

			//- Traverse a node. If intersects a triangle return first
			// intersection point.
			// findAny=true : return any intersection
			// findAny=false: return nearest (to start) intersection
			void traverseNode
			(
				const bool findAny,
				const point& treeStart,
				const vector& treeVec,

				const point& start,
				const point& end,
				const label nodeI,
				const direction octantI,

				pointIndexHit& hitInfo,
				direction& faceID
			) const;

			//- Find any or nearest intersection
			pointIndexHit findLine
			(
				const bool findAny,
				const point& treeStart,
				const point& treeEnd,
				const label startNodeI,
				const direction startOctantI,
				const bool verbose = false
			) const;

			//- Find any or nearest intersection of line between start and end.
			pointIndexHit findLine
			(
				const bool findAny,
				const point& start,
				const point& end
			) const;

			//- Find all elements intersecting box.
			void findBox
			(
				const label nodeI,
				const treeBoundBox& searchBox,
				labelHashSet& elements
			) const;

		// Other

			//- Count number of elements on this and sublevels
			label countElements(const labelBits index) const;

			//- Dump node+octant to an obj file
			void writeOBJ(const label nodeI, const direction octant) const;

			//- From index into contents_ to subNodes_ entry
			static labelBits contentPlusOctant
			(
				const label i,
				const direction octant
			)
			{
				return labelBits(-i - 1, octant);
			}

			//- From index into nodes_ to subNodes_ entry
			static labelBits nodePlusOctant
			(
				const label i,
				const direction octant
			)
			{
				return labelBits(i + 1, octant);
			}

			//- From empty to subNodes_ entry
			static labelBits emptyPlusOctant
			(
				const direction octant
			)
			{
				return labelBits(0, octant);
			}

public:

	// Static functions

		//- Get the perturbation tolerance
		static scalar& perturbTol();


	// Constructors

		//- Construct null
		indexedOctree(const Type& shapes);

		//- Construct from components
		indexedOctree
		(
			const Type& shapes,
			const List<node>& nodes,
			const labelListList& contents
		);

		//- Construct from shapes
		indexedOctree
		(
			const Type& shapes,
			const treeBoundBox& bb,
			const label maxLevels,          // maximum number of levels
			const scalar maxLeafRatio,      // how many elements per leaf
			const scalar maxDuplicity       // in how many leaves is a shape on
					                        // average
		);

		//- Construct from Istream
		indexedOctree(const Type& shapes, Istream& is);

		//- Clone
		autoPtr<indexedOctree<Type> > clone() const
		{
			return autoPtr<indexedOctree<Type> >
			(
				new indexedOctree<Type>(*this)
			);
		}

	// Member Functions

		// Access

			//- Reference to shape
			const Type& shapes() const
			{
				return shapes_;
			}

			//- List of all nodes
			const List<node>& nodes() const
			{
				return nodes_;
			}

			//- List of all contents (referenced by those nodes that are
			//  contents)
			const labelListList& contents() const
			{
				return contents_;
			}

			//- Top bounding box
			const treeBoundBox& bb() const
			{
				if (nodes_.empty())
				{
					FatalErrorIn("indexedOctree<Type>::bb() const")
					    << "Tree is empty" << abort(FatalError);
				}
				return nodes_[0].bb_;
			}


		// Conversions for entries in subNodes_.

			static bool isContent(const labelBits i)
			{
				return i.val() < 0;
			}

			static bool isEmpty(const labelBits i)
			{
				return i.val() == 0;
			}

			static bool isNode(const labelBits i)
			{
				return i.val() > 0;
			}

			static label getContent(const labelBits i)
			{
				if (!isContent(i))
				{
					FatalErrorIn("getContent(const label)")
					    << abort(FatalError);
				}
				return -i.val()-1;
			}

			static label getNode(const labelBits i)
			{
				if (!isNode(i))
				{
					FatalErrorIn("getNode(const label)")
					    << abort(FatalError);
				}
				return i.val() - 1;
			}

			static direction getOctant(const labelBits i)
			{
				return i.bits();
			}


		// Queries

			//- Calculate nearest point on nearest shape. Returns
			//  - bool : any point found nearer than nearestDistSqr
			//  - label: index in shapes
			//  - point: actual nearest point found
			pointIndexHit findNearest
			(
				const point& sample,
				const scalar nearestDistSqr
			) const;

			//- Low level: calculate nearest starting from subnode.
			void findNearest
			(
				const label nodeI,
				const point&,

				scalar& nearestDistSqr,
				label& nearestShapeI,
				point& nearestPoint
			) const;

			//- Find nearest to line. Returns
			//  - bool : any point found?
			//  - label: index in shapes
			//  - point: actual nearest point found
			//  sets:
			//  - linePoint : corresponding nearest point on line
			pointIndexHit findNearest
			(
				const linePointRef& ln,
				treeBoundBox& tightest,
				point& linePoint
			) const;

			//- Find nearest intersection of line between start and end.
			pointIndexHit findLine
			(
				const point& start,
				const point& end
			) const;

			//- Find any intersection of line between start and end.
			pointIndexHit findLineAny
			(
				const point& start,
				const point& end
			) const;

			//- Find (in no particular order) indices of all shapes inside or
			//  overlapping bounding box (i.e. all shapes not outside box)
			labelList findBox(const treeBoundBox& bb) const;

			//- Find deepest node (as parent+octant) containing point. Starts
			//  off from starting index in nodes_ (use 0 to start from top)
			//  Use getNode, getOctant to extract info.
			labelBits findNode(const label nodeI, const point&) const;

			//- Determine type (inside/outside/mixed) for point. unknown if
			//  cannot be determined (e.g. non-manifold surface)
			volumeType getVolumeType(const point&) const;

			//- Helper function to return the side. Returns outside if
			//  outsideNormal&vec >= 0, inside otherwise
			static volumeType getSide
			(
				const vector& outsideNormal,
				const vector& vec
			);

			//- Helper: does bb intersect a sphere around sample? Or is any
			//  corner point of bb closer than nearestDistSqr to sample.
			static bool overlaps
			(
				const point& bbMin,
				const point& bbMax,
				const scalar nearestDistSqr,
				const point& sample
			);


		// Write

			//- Print tree. Either print all indices (printContent = true) or
			//  just size of contents nodes.
			void print
			(
				prefixOSstream&,
				const bool printContents,
				const label
			) const;

			bool write(Ostream& os) const;

	// IOstream Operators

		friend Ostream& operator<< <Type>(Ostream&, const indexedOctree<Type>&);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "indexedOctree.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
