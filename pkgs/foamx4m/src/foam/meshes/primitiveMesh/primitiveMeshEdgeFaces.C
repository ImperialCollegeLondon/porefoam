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

#include "primitiveMesh.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::labelListList& Foam::primitiveMesh::edgeFaces() const
{
	if (!efPtr_)
	{
		if (debug)
		{
			Pout<< "primitiveMesh::edgeFaces() : calculating edgeFaces"
				<< endl;

			if (debug == -1)
			{
				// For checking calls:abort so we can quickly hunt down
				// origin of call
				FatalErrorIn("primitiveMesh::edgeFaces()")
					<< abort(FatalError);
			}
		}

		// Invert faceEdges
		efPtr_ = new labelListList(nEdges());
		invertManyToMany(nEdges(), faceEdges(), *efPtr_);
	}

	return *efPtr_;
}


const Foam::labelList& Foam::primitiveMesh::edgeFaces
(
	const label edgeI,
	dynamicLabelList& storage
) const
{
	if (hasEdgeFaces())
	{
		return edgeFaces()[edgeI];
	}
	else
	{
		// Use the fact that pointEdges are sorted in incrementing edge order
		const edge& e = edges()[edgeI];
		const labelList& pFaces0 = pointFaces()[e[0]];
		const labelList& pFaces1 = pointFaces()[e[1]];

		label i0 = 0;
		label i1 = 0;

		storage.clear();

		while (i0 < pFaces0.size() && i1 < pFaces1.size())
		{
			if (pFaces0[i0] < pFaces1[i1])
			{
				++i0;
			}
			else if (pFaces0[i0] > pFaces1[i1])
			{
				++i1;
			}
			else
			{
				// Equal. Append.
				storage.append(pFaces0[i0]);
				++i0;
				++i1;
			}
		}

		return storage;
	}
}


const Foam::labelList& Foam::primitiveMesh::edgeFaces(const label edgeI) const
{
	return edgeFaces(edgeI, labels_);
}


// ************************************************************************* //
