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
	Foam::PatchTools

Description
	A collection of tools for searching, sorting PrimitivePatch information.

	The class could also be extended to include more that just static methods.

SourceFiles
	PatchTools.C
	PatchToolsCheck.C
	PatchToolsEdgeOwner.C
	PatchToolsSearch.C
	PatchToolsSortEdges.C

\*---------------------------------------------------------------------------*/

#ifndef PatchTools_H
#define PatchTools_H

#include "PrimitivePatchTemplate.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class PatchTools
{
public:

	//- Check for orientation issues.
	//  Returns true if problems were found.
	//  If a normal flips across an edge, places it in the HashSet
	template
	<
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static bool checkOrientation
	(
		const PrimitivePatch<Face, FaceList, PointField, PointType>&,
		const bool report = false,
		labelHashSet* marked = 0
	);


	//- Fill faceZone with currentZone for every face reachable
	//  from faceI without crossing edge marked in borderEdge.
	//  Note: faceZone has to be sized nFaces before calling.
	template
	<
		class BoolListType,
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static void markZone
	(
		const PrimitivePatch<Face, FaceList, PointField, PointType>&,
		const BoolListType& borderEdge,
		const label faceI,
		const label currentZone,
		labelList& faceZone
	);

	//- Size and fills faceZone with zone of face.
	//  Zone is area reachable by edge crossing without crossing borderEdge.
	//  Returns number of zones.
	template
	<
		class BoolListType,
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static label markZones
	(
		const PrimitivePatch<Face, FaceList, PointField, PointType>&,
		const BoolListType& borderEdge,
		labelList& faceZone
	);

	//- Determine the mapping for a sub-patch.
	//  Only include faces for which bool-list entry is true.
	//  @param[in]  includeFaces faces to include
	//  @param[out] pointMap mapping new to old localPoints
	//  @param[out] faceMap  mapping new to old faces
	template
	<
		class BoolListType,
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static void subsetMap
	(
		const PrimitivePatch<Face, FaceList, PointField, PointType>&,
		const BoolListType& includeFaces,
		labelList& pointMap,
		labelList& faceMap
	);

	//- Return edge-face addressing sorted by angle around the edge.
	//  Orientation is anticlockwise looking from edge.vec(localPoints())
	template
	<
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static labelListList sortedEdgeFaces
	(
		const PrimitivePatch<Face, FaceList, PointField, PointType>&
	);


	//- If 2 face neighbours: label of face where ordering of edge
	//  is consistent with righthand walk.
	//  If 1 neighbour: label of only face.
	//  If >2 neighbours: undetermined.
	template
	<
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static labelList edgeOwner
	(
		const PrimitivePatch<Face, FaceList, PointField, PointType>&
	);


	//- Gather points and faces onto master and merge into single patch.
	//  Note: uses faces/points, not localFaces/localPoints.
	template
	<
		class Face,
		template<class> class FaceList,
		class PointField,
		class PointType
	>
	static void gatherAndMerge
	(
		const scalar mergeDist,
		const PrimitivePatch<Face, FaceList, PointField, PointType>& p,
		Field<PointType>& mergedPoints,
		List<Face>& mergedFaces,
		labelList& pointMergeMap
	);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "PatchTools.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
