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

#include "PrimitivePatchTemplate.H"
#include "HashSet.H"

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
calcBdryPoints() const
{
	if (debug)
	{
		Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
			<< "calcBdryPoints() : "
			<< "calculating boundary points"
			<< endl;
	}

	if (boundaryPointsPtr_)
	{
		// it is considered an error to attempt to recalculate
		// if already allocated
		FatalErrorIn
		(
			"PrimitivePatch<Face, FaceList, PointField, PointType>::"
			"calcBdryPoints()"
		)   << "edge types already calculated"
			<< abort(FatalError);
	}

	const edgeList& e = edges();

	labelHashSet bp(2*e.size());

	for (label edgeI = nInternalEdges_; edgeI < e.size(); edgeI++)
	{
		const edge& curEdge = e[edgeI];

		bp.insert(curEdge.start());
		bp.insert(curEdge.end());
	}

	boundaryPointsPtr_ = new labelList(bp.toc());
	sort(*boundaryPointsPtr_);

	if (debug)
	{
		Info<< "PrimitivePatch<Face, FaceList, PointField, PointType>::"
			<< "calcBdryPoints() : "
			<< "finished calculating boundary points"
			<< endl;
	}
}


// ************************************************************************* //
