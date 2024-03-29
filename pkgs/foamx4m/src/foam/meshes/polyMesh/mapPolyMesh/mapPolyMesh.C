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

Description

\*---------------------------------------------------------------------------*/

#include "mapPolyMesh.H"
#include "polyMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from mesh.  No morphing data: the mesh has not changed
// HJ, 27/Nov/2009
Foam::mapPolyMesh::mapPolyMesh(const polyMesh& mesh)
:
	mesh_(mesh),
	morphing_(false),
	nOldPoints_(mesh.nPoints()),
	nOldFaces_(mesh.nFaces()),
	nOldCells_(mesh.nCells())
{}


// Construct from components
Foam::mapPolyMesh::mapPolyMesh
(
	const polyMesh& mesh,
	const label nOldPoints,
	const label nOldFaces,
	const label nOldCells,
	const labelList& pointMap,
	const List<objectMap>& pointsFromPoints,
	const labelList& faceMap,
	const List<objectMap>& facesFromPoints,
	const List<objectMap>& facesFromEdges,
	const List<objectMap>& facesFromFaces,
	const labelList& cellMap,
	const List<objectMap>& cellsFromPoints,
	const List<objectMap>& cellsFromEdges,
	const List<objectMap>& cellsFromFaces,
	const List<objectMap>& cellsFromCells,
	const labelList& reversePointMap,
	const labelList& reverseFaceMap,
	const labelList& reverseCellMap,
	const labelHashSet& flipFaceFlux,
	const labelListList& patchPointMap,
	const labelListList& pointZoneMap,
	const labelListList& faceZonePointMap,
	const labelListList& faceZoneFaceMap,
	const labelListList& cellZoneMap,
	const boolList& resetPatchFlag,
	const pointField& preMotionPoints,
	const labelList& oldPatchStarts,
	const labelList& oldPatchNMeshPoints
)
:
	mesh_(mesh),
	morphing_(true),  // The mesh is morphing
	nOldPoints_(nOldPoints),
	nOldFaces_(nOldFaces),
	nOldCells_(nOldCells),
	pointMap_(pointMap),
	pointsFromPointsMap_(pointsFromPoints),
	faceMap_(faceMap),
	facesFromPointsMap_(facesFromPoints),
	facesFromEdgesMap_(facesFromEdges),
	facesFromFacesMap_(facesFromFaces),
	cellMap_(cellMap),
	cellsFromPointsMap_(cellsFromPoints),
	cellsFromEdgesMap_(cellsFromEdges),
	cellsFromFacesMap_(cellsFromFaces),
	cellsFromCellsMap_(cellsFromCells),
	reversePointMap_(reversePointMap),
	reverseFaceMap_(reverseFaceMap),
	reverseCellMap_(reverseCellMap),
	flipFaceFlux_(flipFaceFlux),
	patchPointMap_(patchPointMap),
	pointZoneMap_(pointZoneMap),
	faceZonePointMap_(faceZonePointMap),
	faceZoneFaceMap_(faceZoneFaceMap),
	cellZoneMap_(cellZoneMap),
	resetPatchFlag_(resetPatchFlag),
	preMotionPoints_(preMotionPoints),
	oldPatchSizes_(oldPatchStarts.size()),
	oldPatchStarts_(oldPatchStarts),
	oldPatchNMeshPoints_(oldPatchNMeshPoints)
{
	// Calculate old patch sizes
	for (label patchI = 0; patchI < oldPatchStarts_.size() - 1; patchI++)
	{
		oldPatchSizes_[patchI] =
			oldPatchStarts_[patchI + 1] - oldPatchStarts_[patchI];
	}

	// Set the last one by hand
	const label lastPatchID = oldPatchStarts_.size() - 1;

	oldPatchSizes_[lastPatchID] = nOldFaces_ - oldPatchStarts_[lastPatchID];

	if (polyMesh::debug)
	{
		if (min(oldPatchSizes_) < 0)
		{
			FatalErrorIn("mapPolyMesh::mapPolyMesh(...)")
				<< "Calculated negative old patch size.  Error in mapping data"
				<< abort(FatalError);
		}
	}
}


// Construct from components and optionally reuse storage
Foam::mapPolyMesh::mapPolyMesh
(
	const polyMesh& mesh,
	const label nOldPoints,
	const label nOldFaces,
	const label nOldCells,
	labelList& pointMap,
	List<objectMap>& pointsFromPoints,
	labelList& faceMap,
	List<objectMap>& facesFromPoints,
	List<objectMap>& facesFromEdges,
	List<objectMap>& facesFromFaces,
	labelList& cellMap,
	List<objectMap>& cellsFromPoints,
	List<objectMap>& cellsFromEdges,
	List<objectMap>& cellsFromFaces,
	List<objectMap>& cellsFromCells,
	labelList& reversePointMap,
	labelList& reverseFaceMap,
	labelList& reverseCellMap,
	labelHashSet& flipFaceFlux,
	labelListList& patchPointMap,
	labelListList& pointZoneMap,
	labelListList& faceZonePointMap,
	labelListList& faceZoneFaceMap,
	labelListList& cellZoneMap,
	boolList& resetPatchFlag,
	pointField& preMotionPoints,
	labelList& oldPatchStarts,
	labelList& oldPatchNMeshPoints,
	const bool reUse
)
:
	mesh_(mesh),
	morphing_(true),  // The mesh is morphing
	nOldPoints_(nOldPoints),
	nOldFaces_(nOldFaces),
	nOldCells_(nOldCells),
	pointMap_(pointMap, reUse),
	pointsFromPointsMap_(pointsFromPoints, reUse),
	faceMap_(faceMap, reUse),
	facesFromPointsMap_(facesFromPoints, reUse),
	facesFromEdgesMap_(facesFromEdges, reUse),
	facesFromFacesMap_(facesFromFaces, reUse),
	cellMap_(cellMap, reUse),
	cellsFromPointsMap_(cellsFromPoints, reUse),
	cellsFromEdgesMap_(cellsFromEdges, reUse),
	cellsFromFacesMap_(cellsFromFaces, reUse),
	cellsFromCellsMap_(cellsFromCells, reUse),
	reversePointMap_(reversePointMap, reUse),
	reverseFaceMap_(reverseFaceMap, reUse),
	reverseCellMap_(reverseCellMap, reUse),
	flipFaceFlux_(flipFaceFlux),
	patchPointMap_(patchPointMap, reUse),
	pointZoneMap_(pointZoneMap, reUse),
	faceZonePointMap_(faceZonePointMap, reUse),
	faceZoneFaceMap_(faceZoneFaceMap, reUse),
	cellZoneMap_(cellZoneMap, reUse),
	resetPatchFlag_(resetPatchFlag, reUse),
	preMotionPoints_(preMotionPoints, reUse),
	oldPatchSizes_(oldPatchStarts.size()),
	oldPatchStarts_(oldPatchStarts, reUse),
	oldPatchNMeshPoints_(oldPatchNMeshPoints, reUse)
{
	// Calculate old patch sizes
	for (label patchI = 0; patchI < oldPatchStarts_.size() - 1; patchI++)
	{
		oldPatchSizes_[patchI] =
			oldPatchStarts_[patchI + 1] - oldPatchStarts_[patchI];
	}

	// Set the last one by hand
	const label lastPatchID = oldPatchStarts_.size() - 1;

	oldPatchSizes_[lastPatchID] = nOldFaces_ - oldPatchStarts_[lastPatchID];

	if (polyMesh::debug)
	{
		if (min(oldPatchSizes_) < 0)
		{
			FatalErrorIn("mapPolyMesh::mapPolyMesh(...)")
				<< "Calculated negative old patch size.  Error in mapping data"
				<< abort(FatalError);
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// ************************************************************************* //
