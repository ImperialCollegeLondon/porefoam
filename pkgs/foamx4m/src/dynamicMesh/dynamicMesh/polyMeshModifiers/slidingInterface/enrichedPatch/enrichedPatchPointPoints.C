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
	Enriched patch addressing

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.  Copyright Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#include "enrichedPatch.H"
#include "primitiveMesh.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcPointPoints() const
{
	// Calculate point-point addressing
	if (pointPointsPtr_)
	{
		FatalErrorIn("void enrichedPatch::calcPointPoints() const")
			<< "Point-point addressing already calculated."
			<< abort(FatalError);
	}

	// Algorithm:
	// Go through all faces and add the previous and next point as the
	// neighbour for each point. While inserting points, reject the
	// duplicates (as every internal edge will be visited twice).
	List<DynamicList<label, primitiveMesh::edgesPerPoint_> >
		pp(meshPoints().size());

	const faceList& lf = localFaces();

	 bool found = false;

	forAll (lf, faceI)
	{
		const face& curFace = lf[faceI];

		forAll (curFace, pointI)
		{
			DynamicList<label, primitiveMesh::edgesPerPoint_>&
				curPp = pp[curFace[pointI]];

			// Do next label
			label next = curFace.nextLabel(pointI);

			found = false;

			forAll (curPp, i)
			{
				if (curPp[i] == next)
				{
					found = true;
					break;
				}
			}

			if (!found)
			{
				curPp.append(next);
			}

			// Do previous label
			label prev = curFace.prevLabel(pointI);
			found = false;

			forAll (curPp, i)
			{
				if (curPp[i] == prev)
				{
					found = true;
					break;
				}
			}

			if (!found)
			{
				curPp.append(prev);
			}
		}
	}

	// Re-pack the list
	pointPointsPtr_ = new labelListList(pp.size());
	labelListList& ppAddr = *pointPointsPtr_;

	forAll (pp, pointI)
	{
		ppAddr[pointI].transfer(pp[pointI].shrink());
	}
}


// ************************************************************************* //
