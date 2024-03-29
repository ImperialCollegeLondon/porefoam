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

\*---------------------------------------------------------------------------*/

#include "inversePointDistanceDiffusivity.H"
#include "addToRunTimeSelectionTable.H"
#include "HashSet.H"
#include "pointEdgePoint.H"
#include "PointEdgeWave.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(inversePointDistanceDiffusivity, 0);

	addToRunTimeSelectionTable
	(
		motionDiffusivity,
		inversePointDistanceDiffusivity,
		Istream
	);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inversePointDistanceDiffusivity::inversePointDistanceDiffusivity
(
	const fvMotionSolver& mSolver,
	Istream& mdData
)
:
	uniformDiffusivity(mSolver, mdData),
	patchNames_(mdData)
{
	correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::inversePointDistanceDiffusivity::~inversePointDistanceDiffusivity()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inversePointDistanceDiffusivity::correct()
{
	const polyMesh& mesh = mSolver().mesh();
	const polyBoundaryMesh& bdry = mesh.boundaryMesh();

	labelHashSet patchSet(bdry.size());

	label nPatchEdges = 0;

	forAll (patchNames_, i)
	{
		label pID = bdry.findPatchID(patchNames_[i]);

		if (pID > -1)
		{
			patchSet.insert(pID);
			nPatchEdges += bdry[pID].nEdges();
		}
	}

	// Distance to wall on points and edges.
	List<pointEdgePoint> pointWallDist(mesh.nPoints());
	List<pointEdgePoint> edgeWallDist(mesh.nEdges());

	{
		// Seeds
		List<pointEdgePoint> seedInfo(nPatchEdges);
		labelList seedPoints(nPatchEdges);

		nPatchEdges = 0;

		forAllConstIter(labelHashSet, patchSet, iter)
		{
			const polyPatch& patch = bdry[iter.key()];

			const labelList& meshPoints = patch.meshPoints();

			forAll(meshPoints, i)
			{
				label pointI = meshPoints[i];

				if (!pointWallDist[pointI].valid())
				{
					// Not yet seeded
					seedInfo[nPatchEdges] = pointEdgePoint
					(
						mesh.points()[pointI],
						0.0
					);
					seedPoints[nPatchEdges] = pointI;
					pointWallDist[pointI] = seedInfo[nPatchEdges];

					nPatchEdges++;
				}
			}
		}
		seedInfo.setSize(nPatchEdges);
		seedPoints.setSize(nPatchEdges);

		// Do calculations
		PointEdgeWave<pointEdgePoint> waveInfo
		(
			mesh,
			seedPoints,
			seedInfo,

			pointWallDist,
			edgeWallDist,
			mesh.globalData().nTotalPoints() // max iterations
		);
	}


	for (label faceI=0; faceI<mesh.nInternalFaces(); faceI++)
	{
		const face& f = mesh.faces()[faceI];

		scalar dist = 0;

		forAll(f, fp)
		{
			dist += sqrt(pointWallDist[f[fp]].distSqr());
		}
		dist /= f.size();

		faceDiffusivity_[faceI] = 1.0/dist;
	}

	forAll(faceDiffusivity_.boundaryField(), patchI)
	{
		fvsPatchScalarField& bfld = faceDiffusivity_.boundaryField()[patchI];

		if (patchSet.found(patchI))
		{
			const unallocLabelList& faceCells = bfld.patch().faceCells();

			forAll(bfld, i)
			{
				const cell& ownFaces = mesh.cells()[faceCells[i]];

				labelHashSet cPoints(4*ownFaces.size());

				scalar dist = 0;

				forAll(ownFaces, ownFaceI)
				{
					const face& f = mesh.faces()[ownFaces[ownFaceI]];

					forAll(f, fp)
					{
						if (cPoints.insert(f[fp]))
						{
							dist += sqrt(pointWallDist[f[fp]].distSqr());
						}
					}
				}
				dist /= cPoints.size();

				bfld[i] = 1.0/dist;
			}
		}
		else
		{
			label start = bfld.patch().patch().start();

			forAll(bfld, i)
			{
				const face& f = mesh.faces()[start+i];

				scalar dist = 0;

				forAll(f, fp)
				{
					dist += sqrt(pointWallDist[f[fp]].distSqr());
				}
				dist /= f.size();

				bfld[i] = 1.0/dist;
			}
		}
	}
}


// ************************************************************************* //
