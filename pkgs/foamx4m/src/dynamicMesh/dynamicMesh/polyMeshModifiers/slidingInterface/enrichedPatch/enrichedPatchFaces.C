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
	Calculate enriched faces

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.  Copyright Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#include "enrichedPatch.H"
#include "DynamicList.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::enrichedPatch::enrichedFaceRatio_ = 3;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcEnrichedFaces
(
	const labelListList& pointsIntoMasterEdges,
	const labelListList& pointsIntoSlaveEdges,
	const pointField& projectedSlavePoints
)
{
	if (enrichedFacesPtr_)
	{
		FatalErrorIn
		(
			"void enrichedPatch::calcEnrichedFaces\n"
			"(\n"
			"    const labelListList& pointsIntoMasterEdges,\n"
			"    const labelListList& pointsIntoSlaveEdges,\n"
			"    const pointField& projectedSlavePoints\n"
			")"
		)   << "Enriched faces already calculated."
			<< abort(FatalError);
	}

	// Create a list of enriched faces
	// Algorithm:
	// 1) Grab the original face and start from point zero.
	// 2) If the point has been merged away, grab the merge label;
	//	otherwise, keep the original label.
	// 3) Go to the next edge. Collect all the points to be added along
	//	the edge; order them in the necessary direction and insert onto the
	//	face.
	// 4) Grab the next point and return on step 2.
	enrichedFacesPtr_ = new faceList(masterPatch_.size() + slavePatch_.size());
	faceList& enrichedFaces = *enrichedFacesPtr_;

	label nEnrichedFaces = 0;

	const pointField& masterLocalPoints = masterPatch_.localPoints();
	const faceList& masterLocalFaces = masterPatch_.localFaces();
	const labelListList& masterFaceEdges = masterPatch_.faceEdges();

	const faceList& slaveLocalFaces = slavePatch_.localFaces();
	const labelListList& slaveFaceEdges = slavePatch_.faceEdges();

	// For correct functioning of the enrichedPatch class, the slave
	// faces need to be inserted first.  See comments in
	// enrichedPatch.H

	// Get reference to the point merge map
	const Map<label>& pmm = pointMergeMap();

	// Add slave faces into the enriched faces list

	forAll (slavePatch_, faceI)
	{
		const face oldFace = slavePatch_[faceI];
		const face oldLocalFace = slaveLocalFaces[faceI];
//		 Pout<< "old slave face " << faceI << ": " << oldFace
//			 << oldLocalFace.points(slavePatch_.localPoints()) << endl;
		const labelList& curEdges = slaveFaceEdges[faceI];

		dynamicLabelList newFace(oldFace.size()*enrichedFaceRatio_);

		// Note: The number of points and edges in a face is always identical
		// so both can be done is the same loop
		forAll (oldFace, i)
		{
			// Add the point
			Map<label>::const_iterator mpIter =
				pmm.find(oldFace[i]);

			if (mpIter == pmm.end())
			{
				// Point not mapped
				newFace.append(oldFace[i]);

				// Add the projected point into the patch support
				if (!pointMap().found(oldFace[i]))
				{
//					 Pout << "Inserting point from face " << i << " " << oldFace[i] << projectedSlavePoints[oldLocalFace[i]] << endl;
					pointMap().insert
					(
						oldFace[i],	// Global label of point
						projectedSlavePoints[oldLocalFace[i]] // Projected posn
					);
				}
			}
			else
			{
				// Point mapped
				newFace.append(mpIter());

				// Add the projected point into the patch support
				if (!pointMap().found(mpIter()))
				{
//					 Pout << "Inserting point " << mpIter() << " at " << projectedSlavePoints[oldLocalFace[i]] << endl;
					pointMap().insert
					(
						mpIter(),	// Merged global label of point
						projectedSlavePoints[oldLocalFace[i]] // Projected posn
					);
				}
			}

			// Grab the edge points

			const labelList& slavePointsOnEdge =
				pointsIntoSlaveEdges[curEdges[i]];
//			 Pout << "slavePointsOnEdge for " << curEdges[i] << ": " << slavePointsOnEdge << endl;
			// If there are no points on the edge, skip everything
			// If there is only one point, no need for sorting
			if (slavePointsOnEdge.size() > 0)
			{
				// Sort edge points in order
				scalarField edgePointWeights(slavePointsOnEdge.size());
				const point& startPoint = projectedSlavePoints[oldLocalFace[i]];

				vector e =
					projectedSlavePoints[oldLocalFace.nextLabel(i)]
				  - startPoint;

				scalar magSqrE = magSqr(e);

				if (magSqrE > SMALL)
				{
					e /= magSqrE;
				}
				else
				{
					FatalErrorIn
					(
						"void enrichedPatch::calcEnrichedFaces\n"
						"(\n"
						"    const labelListList& pointsIntoMasterEdges,\n"
						"    const labelListList& pointsIntoSlaveEdges,\n"
						"    const pointField& projectedSlavePoints\n"
						")"
					)   << "Zero length edge in slave patch for face " << i
						<< ".  This is not allowed."
						<< abort(FatalError);
				}

				pointField slavePosOnEdge(slavePointsOnEdge.size());

				forAll (slavePointsOnEdge, edgePointI)
				{
					slavePosOnEdge[edgePointI] =
						pointMap().find(slavePointsOnEdge[edgePointI])();

					edgePointWeights[edgePointI] =
						(e & (slavePosOnEdge[edgePointI] - startPoint));
				}

				if (debug)
				{
					// Check weights: all new points should be on the edge
					if (min(edgePointWeights) < 0 || max(edgePointWeights) > 1)
					{
						FatalErrorIn
						(
							"void enrichedPatch::calcEnrichedFaces\n"
							"(\n"
							"    const labelListList& pointsIntoMasterEdges,\n"
							"    const labelListList& pointsIntoSlaveEdges,\n"
							"    const pointField& projectedSlavePoints\n"
							")"
						)   << "Invalid point edge weights.  Some of points are"
							<< " not on the edge for edge " << curEdges[i]
							<< " of face " << faceI << " in slave patch." << nl
							<< "Min weight: " << min(edgePointWeights)
							<< " Max weight: " << max(edgePointWeights)
							<< abort(FatalError);
					}
				}

				// Go through the points and collect them based on
				// weights from lower to higher.  This gives the
				// correct order of points along the edge.
				for (label passI = 0; passI < edgePointWeights.size(); passI++)
				{
					// Max weight can only be one, so the sorting is
					// done by elimination.
					label nextPoint = -1;
					scalar dist = 2;

					forAll (edgePointWeights, wI)
					{
						if (edgePointWeights[wI] < dist)
						{
							dist = edgePointWeights[wI];
							nextPoint = wI;
						}
					}

					// Insert the next point and reset its weight to exclude it
					// from future picks
					newFace.append(slavePointsOnEdge[nextPoint]);
					edgePointWeights[nextPoint] = GREAT;

					// Add the point into patch support if not already there
					if (!pointMap().found(slavePointsOnEdge[nextPoint]))
					{
//						 Pout << "Inserting point, pass 1: " << slavePointsOnEdge[nextPoint] << " " << slavePosOnEdge[nextPoint] << endl;
						pointMap().insert
						(
							slavePointsOnEdge[nextPoint],
							slavePosOnEdge[nextPoint]
						);
					}
				}
			}
		}
//		 Pout<< "New slave face " << faceI << ": " << newFace << endl;

		// Add the new face to the list
		enrichedFaces[nEnrichedFaces].transfer(newFace.shrink());
		nEnrichedFaces++;
	}

	// Add master faces into the enriched faces list

	forAll (masterPatch_, faceI)
	{
		const face& oldFace = masterPatch_[faceI];
		const face& oldLocalFace = masterLocalFaces[faceI];
//		 Pout << "old master face: " << oldFace << endl;
		const labelList& curEdges = masterFaceEdges[faceI];

		dynamicLabelList newFace(oldFace.size()*enrichedFaceRatio_);

		// Note: The number of points and edges in a face is always identical
		// so both can be done is the same loop
		forAll (oldFace, i)
		{
			// Add the point
			Map<label>::const_iterator mpIter =
				pmm.find(oldFace[i]);

			if (mpIter == pmm.end())
			{
				// Point not mapped
				newFace.append(oldFace[i]);

				// Add the point into patch support
				if (!pointMap().found(oldFace[i]))
				{
//					 Pout << "Inserting master point: " << oldFace[i] << " " <<
//						 masterLocalPoints[oldLocalFace[i]] << endl;
					pointMap().insert
					(
						oldFace[i],
						masterLocalPoints[oldLocalFace[i]]
					);
				}
			}
			else
			{
				// Point mapped
				newFace.append(mpIter());

				// Add the point into support
				if (!pointMap().found(mpIter()))
				{
//					 Pout << "Inserting master 2 " << mpIter() << " " << masterLocalPoints[oldLocalFace[i]] << endl;
					pointMap().insert
					(
						mpIter(),
						masterLocalPoints[oldLocalFace[i]]
					);
				}
			}

			// Grab the edge points

			const labelList& masterPointsOnEdge =
				pointsIntoMasterEdges[curEdges[i]];

			// If there are no points on the edge, skip everything
			// If there is only one point, no need for sorting
			if (masterPointsOnEdge.size() > 0)
			{
				// Sort edge points in order
				scalarField edgePointWeights(masterPointsOnEdge.size());
				const point& startPoint = masterLocalPoints[oldLocalFace[i]];

				vector e =
					masterLocalPoints[oldLocalFace.nextLabel(i)]
				  - startPoint;

				scalar magSqrE = magSqr(e);

				if (magSqrE > SMALL)
				{
					e /= magSqrE;
				}
				else
				{
					FatalErrorIn
					(
						"void enrichedPatch::calcEnrichedFaces\n"
						"(\n"
						"    const labelListList& pointsIntoMasterEdges,\n"
						"    const labelListList& pointsIntoSlaveEdges,\n"
						"    const pointField& projectedSlavePoints\n"
						")"
					)   << "Zero length edge in master patch for face " << i
						<< ".  This is not allowed."
						<< abort(FatalError);
				}

				pointField masterPosOnEdge(masterPointsOnEdge.size());

				forAll (masterPointsOnEdge, edgePointI)
				{
					masterPosOnEdge[edgePointI] =
						pointMap().find(masterPointsOnEdge[edgePointI])();

					edgePointWeights[edgePointI] =
						(e & (masterPosOnEdge[edgePointI] - startPoint));
				}

				if (debug)
				{
					// Check weights: all new points should be on the edge
					if (min(edgePointWeights) < 0 || max(edgePointWeights) > 1)
					{
						FatalErrorIn
						(
							"void enrichedPatch::calcEnrichedFaces\n"
							"(\n"
							"    const labelListList& pointsIntoMasterEdges,\n"
							"    const labelListList& pointsIntoSlaveEdges,\n"
							"    const pointField& projectedSlavePoints\n"
							")"
						)   << "Invalid point edge weights.  Some of points are"
							<< " not on the edge for edge " << curEdges[i]
							<< " of face " << faceI << " in master patch." << nl
							<< "Min weight: " << min(edgePointWeights)
							<< " Max weight: " << max(edgePointWeights)
							<< abort(FatalError);
					}
				}

				// Go through the points and collect them based on
				// weights from lower to higher.  This gives the
				// correct order of points along the edge.
				for (label pass = 0; pass < edgePointWeights.size(); pass++)
				{
					// Max weight can only be one, so the sorting is
					// done by elimination.
					label nextPoint = -1;
					scalar dist = 2;

					forAll (edgePointWeights, wI)
					{
						if (edgePointWeights[wI] < dist)
						{
							dist = edgePointWeights[wI];
							nextPoint = wI;
						}
					}

					// Insert the next point and reset its weight to exclude it
					// from future picks
					newFace.append(masterPointsOnEdge[nextPoint]);
					edgePointWeights[nextPoint] = GREAT;

					// Add the point into patch support
					if (!pointMap().found(masterPointsOnEdge[nextPoint]))
					{
						pointMap().insert
						(
							masterPointsOnEdge[nextPoint],
							masterPosOnEdge[nextPoint]
						);
					}
				}
			}
		}
//		 Pout << "New master face: " << newFace << endl;

		// Add the new face to the list
		enrichedFaces[nEnrichedFaces].transfer(newFace.shrink());
		nEnrichedFaces++;
	}

	// Check the support for the enriched patch
	if (debug)
	{
		if (!checkSupport())
		{
			Info<< "Enriched patch support OK. Slave faces: "
				<< slavePatch_.size() << " Master faces: "
				<< masterPatch_.size() << endl;
		}
		else
		{
			FatalErrorIn
			(
				"void enrichedPatch::calcEnrichedFaces\n"
				"(\n"
				"    const labelListList& pointsIntoMasterEdges,\n"
				"    const labelListList& pointsIntoSlaveEdges,\n"
				"    const pointField& projectedSlavePoints\n"
				")"
			)   << "Error in enriched patch support"
				<< abort(FatalError);
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::faceList& Foam::enrichedPatch::enrichedFaces() const
{
	if (!enrichedFacesPtr_)
	{
		FatalErrorIn("const faceList& enrichedPatch::enrichedFaces() const")
			<< "Enriched faces not available yet.  Please use "
			<< "void enrichedPatch::calcEnrichedFaces\n"
			<< "(\n"
			<< "    const labelListList& pointsIntoMasterEdges,\n"
			<< "    const labelListList& pointsIntoSlaveEdges,\n"
			<< "    const pointField& projectedSlavePoints\n"
			<< ")"
			<< " before trying to access faces."
			<< abort(FatalError);
	}

	return *enrichedFacesPtr_;
}


// ************************************************************************* //
