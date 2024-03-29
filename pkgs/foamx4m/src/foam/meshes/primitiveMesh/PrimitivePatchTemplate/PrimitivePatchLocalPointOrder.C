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
	 Orders the local points on the patch for most efficient search

\*---------------------------------------------------------------------------*/

#include "SLList.H"
#include "boolList.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template
<
	class Face,
	template<class> class FaceList,
	class PointField,
	class PointType
>
void
Foam::PrimitivePatch<Face, FaceList, PointField, PointType>::
calcLocalPointOrder() const
{
	// Note: Cannot use bandCompressing as point-point addressing does
	// not exist and is not considered generally useful.
	//

	if (debug)
	{
		Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
			<< "calcLocalPointOrder() : "
			<< "calculating local point order"
			<< endl;
	}

	if (localPointOrderPtr_)
	{
		// it is considered an error to attempt to recalculate
		// if already allocated
		FatalErrorIn
		(
			"PrimitivePatch<Face, FaceList, PointField, PointType>::"
			"calcLocalPointOrder()"
		)   << "local point order already calculated"
			<< abort(FatalError);
	}

	const List<Face>& lf = localFaces();

	const labelListList& ff = faceFaces();

	boolList visitedFace(lf.size(), false);

	localPointOrderPtr_ = new labelList(meshPoints().size(), -1);

	labelList& pointOrder = *localPointOrderPtr_;

	boolList visitedPoint(pointOrder.size(), false);

	label nPoints = 0;

	forAll (lf, faceI)
	{
		if (!visitedFace[faceI])
		{
			SLList<label> faceOrder(faceI);

			do
			{
				const label curFace = faceOrder.first();

				faceOrder.removeHead();

				if (!visitedFace[curFace])
				{
					visitedFace[curFace] = true;

					const labelList& curPoints = lf[curFace];

					// mark points
					forAll (curPoints, pointI)
					{
					    if (!visitedPoint[curPoints[pointI]])
					    {
					        visitedPoint[curPoints[pointI]] = true;

					        pointOrder[nPoints] = curPoints[pointI];

					        nPoints++;
					    }
					}

					// add face neighbours to the list
					const labelList& nbrs = ff[curFace];

					forAll (nbrs, nbrI)
					{
					    if (!visitedFace[nbrs[nbrI]])
					    {
					        faceOrder.append(nbrs[nbrI]);
					    }
					}
				}
			} while (faceOrder.size());
		}
	}

	if (debug)
	{
		Pout<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
			<< "calcLocalPointOrder() "
			<< "finished calculating local point order"
			<< endl;
	}
}


// ************************************************************************* //
