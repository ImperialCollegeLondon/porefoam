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
	Create the list of loops of outside vertices. Goes wrong on multiply
	connected edges (loops will be unclosed).

\*---------------------------------------------------------------------------*/

#include "PrimitivePatchTemplate.H"


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template
<
	class Face,
	template<class> class FaceList,
	class PointField,
	class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcEdgeLoops() const
{
	if (debug)
	{
		Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
			<< "calcEdgeLoops() : "
			<< "calculating boundary edge loops"
			<< endl;
	}

	if (edgeLoopsPtr_)
	{
		// it is considered an error to attempt to recalculate
		// if already allocated
		FatalErrorIn
		(
			"PrimitivePatch<Face, FaceList, PointField, PointType>::"
			"calcEdgeLoops()"
		)   << "edge loops already calculated"
			<< abort(FatalError);
	}

	const edgeList& patchEdges = edges();
	label nIntEdges = nInternalEdges();
	label nBdryEdges = patchEdges.size() - nIntEdges;

	if (nBdryEdges == 0)
	{
		edgeLoopsPtr_ = new labelListList(0);
		return;
	}

	const labelListList& patchPointEdges = pointEdges();


	//
	// Walk point-edge-point and assign loop number
	//

	// Loop per (boundary) edge.
	labelList loopNumber(nBdryEdges, -1);

	// Size return list plenty big
	edgeLoopsPtr_ = new labelListList(nBdryEdges);
	labelListList& edgeLoops = *edgeLoopsPtr_;


	// Current loop number.
	label loopI = 0;

	while (true)
	{
		// Find edge not yet given a loop number.
		label currentEdgeI = -1;

		for (label edgeI = nIntEdges; edgeI < patchEdges.size(); edgeI++)
		{
			if (loopNumber[edgeI-nIntEdges] == -1)
			{
				currentEdgeI = edgeI;
				break;
			}
		}

		if (currentEdgeI == -1)
		{
			// Did not find edge not yet assigned a loop number so done all.
			break;
		}

		// Temporary storage for vertices of current loop
		dynamicLabelList loop(nBdryEdges);

		// Walk from first all the way round, assigning loops
		label currentVertI = patchEdges[currentEdgeI].start();

		do
		{
			loop.append(currentVertI);

			loopNumber[currentEdgeI - nIntEdges] = loopI;

			// Step to next vertex
			currentVertI = patchEdges[currentEdgeI].otherVertex(currentVertI);

			// Step to next (unmarked, boundary) edge.
			const labelList& curEdges = patchPointEdges[currentVertI];

			currentEdgeI = -1;

			forAll(curEdges, pI)
			{
				label edgeI = curEdges[pI];

				if (edgeI >= nIntEdges && (loopNumber[edgeI - nIntEdges] == -1))
				{
					// Unassigned boundary edge.
					currentEdgeI = edgeI;

					break;
				}
			}
		}
		while (currentEdgeI != -1);

		// Done all for current loop. Transfer to edgeLoops.
		edgeLoops[loopI].transfer(loop);

		loopI++;
	}

	edgeLoops.setSize(loopI);

	if (debug)
	{
		Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
			<< "calcEdgeLoops() : "
			<< "finished calculating boundary edge loops"
			<< endl;
	}
}


template
<
	class Face,
	template<class> class FaceList,
	class PointField,
	class PointType
>
const Foam::labelListList&
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
edgeLoops() const
{
	if (!edgeLoopsPtr_)
	{
		calcEdgeLoops();
	}

	return *edgeLoopsPtr_;
}


// ************************************************************************* //
