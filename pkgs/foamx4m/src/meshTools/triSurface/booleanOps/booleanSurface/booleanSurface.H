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
	Foam::booleanSurface

Description
	Surface-surface intersection.
	Given two surfaces construct a combined surface.

	Called 'boolean' since the volume of resulting surface will encompass
	the volumes of the original surface according to some boolean operation:
	- all which is in surface1 AND in surface2 (intersection)
	- all which is in surface1 AND NOT in surface2 ('chop' out surface2)
	- all which is in surface1 OR in surface2 (union)

	Algorithm:
	-# find edge-surface intersection. Class 'surfaceIntersection'.
	-# combine intersection with both surfaces. Class 'intersectedSurface'.
	-# subset surfaces upto intersection. The 'side' of the surface to
	   include is based on the faces that can be reached from a
	   user-supplied face index.
	-# merge surfaces. Only the points on the intersection are shared.

SourceFiles
	booleanSurface.C

\*---------------------------------------------------------------------------*/

#ifndef booleanSurface_H
#define booleanSurface_H

#include "triSurface.H"
#include "surfaceIntersection.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurfaceSearch;
class intersectedSurface;


class booleanSurface
:
	public triSurface
{
	// Data types

		//- Enumeration listing the status of a face (visible/invisible from
		//  outside)
		enum sideStat
		{
			UNVISITED,
			OUTSIDE,
			INSIDE
		};

	// Private data

		//- From new to old face + surface:
		//     >0 : to face on surface1
		//     <0 : to face on surface2. Negate and offset by one to get
		//          face2 (e.g. face2I = -faceMap[]-1)
		labelList faceMap_;

	// Private Member Functions

		//- Check whether subset of faces (from markZones) reaches up to
		//  the intersection.
		static void checkIncluded
		(
			const intersectedSurface& surf,
			const labelList& faceZone,
			const label includedFace
		);

		//- Get label in elems of elem.
		static label index(const labelList& elems, const label elem);

		//- Find index of edge e in subset edgeLabels.
		static label findEdge
		(
			const edgeList& edges,
			const labelList& edgeLabels,
			const edge& e
		);

		//- Get index of face in zoneI whose faceCentre is nearest farAwayPoint
		static label findNearest
		(
			const triSurface& surf,
			const labelList& faceZone,
			const label zoneI
		);

		//- Generate combined patchList (returned). Sets patchMap to map from
		// surf region numbers into combined patchList
		static geometricSurfacePatchList mergePatches
		(
			const triSurface& surf1,
			const triSurface& surf2,
			labelList& patchMap2
		);

		//- On edgeI, coming from face prevFace, determines visibility/side of
		// all the other faces using the edge.
		static void propagateEdgeSide
		(
			const triSurface& surf,
			const label prevVert0,
			const label prevFaceI,
			const label prevState,
			const label edgeI,
			labelList& side
		);

		//- Given in/outside status of face determines status for all
		//  neighbouring faces.
		static void propagateSide
		(
			const triSurface& surf,
			const label prevState,
			const label faceI,
			labelList& side
		);

public:

	ClassName("booleanSurface");


	// Data types

		//- Enumeration listing the possible volume operator types
		enum booleanOpType
		{
			OR,     // Union of volumes
			AND,    // Intersection of volumes
			XOR,    // Difference of volumes
			ALL     // Special: does not subset combined surface. (Produces
					// multiply connected surface)
		};


	// Constructors

		//- Construct null
		booleanSurface();

		//- Construct from surfaces and face labels to keep.
		//  Walks from provided seed faces without crossing intersection line
		//  to determine faces to keep.
		booleanSurface
		(
			const triSurface& surf1,
			const triSurface& surf2,
			const surfaceIntersection& inter,
			const label includeFace1,
			const label includeFace2
		);

		//- Construct from surfaces and operation. Surfaces need to be closed
		//  for this to make any sense since uses inside/outside to determine
		//  which part of combined surface to include.
		booleanSurface
		(
			const triSurface& surf1,
			const triSurface& surf2,
			const surfaceIntersection& inter,
			const label booleanOp
		);


	// Member Functions

		//- new to old face map. >0: surface 1 face label. <0: surface 2. Negate
		//  and subtract 1 to get face label on surface 2.
		const labelList& faceMap() const
		{
			return faceMap_;
		}

		bool from1(const label faceI) const
		{
			return faceMap_[faceI] >= 0;
		}

		bool surf1Face(const label faceI) const
		{
			if (!from1(faceI))
			{
				FatalErrorIn("booleanSurface::surf1Face(const label)")
					<< "face " << faceI << " not from surface 1"
					<< abort(FatalError);
			}
			return faceMap_[faceI];
		}

		bool surf2Face(const label faceI) const
		{
			if (from1(faceI))
			{
				FatalErrorIn("booleanSurface::surf2Face(const label)")
					<< "face " << faceI << " not from surface 2"
					<< abort(FatalError);
			}
			return -faceMap_[faceI]-1;
		}


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
