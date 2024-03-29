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

#include "cellMatcher.H"

#include "primitiveMesh.H"
#include "Map.H"
#include "faceList.H"
#include "labelList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::labelList Foam::cellMatcher::makeIdentity(const label nElems)
{
	labelList result(nElems);

	forAll(result, elemI)
	{
		result[elemI] = elemI;
	}
	return result;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from components
Foam::cellMatcher::cellMatcher
(
	const label vertPerCell,
	const label facePerCell,
	const label maxVertPerFace,
	const word& cellModelName
)
:
	localPoint_(100),
	localFaces_(facePerCell),
	faceSize_(facePerCell, -1),
	pointMap_(vertPerCell),
	faceMap_(facePerCell),
	edgeFaces_(2*vertPerCell*vertPerCell),
	pointFaceIndex_(vertPerCell),
	vertLabels_(vertPerCell),
	faceLabels_(facePerCell),
	cellModelName_(cellModelName),
	cellModelPtr_(nullptr)
{
	forAll(localFaces_, faceI)
	{
		face& f = localFaces_[faceI];

		f.setSize(maxVertPerFace);
	}

	forAll(pointFaceIndex_, vertI)
	{
		pointFaceIndex_[vertI].setSize(facePerCell);
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Create localFaces_ , pointMap_ , faceMap_
Foam::label Foam::cellMatcher::calcLocalFaces
(
	const faceList& faces,
	const labelList& myFaces
)
{
	// Clear map from global to cell numbering
	localPoint_.clear();

	// Renumber face vertices and insert directly into localFaces_
	label newVertI = 0;
	forAll(myFaces, myFaceI)
	{
		label faceI = myFaces[myFaceI];

		const face& f = faces[faceI];
		face& localFace = localFaces_[myFaceI];

		// Size of localFace
		faceSize_[myFaceI] = f.size();

		forAll(f, localVertI)
		{
			label vertI = f[localVertI];

			Map<label>::iterator iter = localPoint_.find(vertI);
			if (iter == localPoint_.end())
			{
				// Not found. Assign local vertex number.

				if (newVertI >= pointMap_.size())
				{
					// Illegal face: more unique vertices than vertPerCell
					return -1;
				}

				localFace[localVertI] = newVertI;
				localPoint_.insert(vertI, newVertI);
				newVertI++;
			}
			else
			{
				// Reuse local vertex number.
				localFace[localVertI] = *iter;
			}
		}

		// Create face from localvertex labels
		faceMap_[myFaceI] = faceI;
	}

	// Create local to global vertex mapping
	for
	(
		Map<label>::iterator iter = localPoint_.begin();
		iter != localPoint_.end();
		++iter
	)
	{
		label fp = iter();
		pointMap_[fp] = iter.key();
	}

	////debug
	//write(Info);

	return newVertI;
}


// Create edgeFaces_ : map from edge to two localFaces for single cell.
void Foam::cellMatcher::calcEdgeAddressing(const label numVert)
{
	edgeFaces_ = -1;

	forAll(localFaces_, localFaceI)
	{
		const face& f = localFaces_[localFaceI];

		label prevVertI = faceSize_[localFaceI] - 1;
		//forAll(f, fp)
		for
		(
			label fp = 0;
			fp < faceSize_[localFaceI];
			fp++
		)
		{
			label start = f[prevVertI];
			label end = f[fp];

			label key1 = edgeKey(numVert, start, end);
			label key2 = edgeKey(numVert, end, start);

			if (edgeFaces_[key1] == -1)
			{
				// Entry key1 unoccupied. Store both permutations.
				edgeFaces_[key1] = localFaceI;
				edgeFaces_[key2] = localFaceI;
			}
			else if (edgeFaces_[key1+1] == -1)
			{
				// Entry key1+1 unoccupied
				edgeFaces_[key1+1] = localFaceI;
				edgeFaces_[key2+1] = localFaceI;
			}
			else
			{
				FatalErrorIn
				(
					"calcEdgeAddressing"
					"(const faceList&, const label)"
				)   << "edgeFaces_ full at entry:" << key1
					<< " for edge " << start << " " << end
					<< abort(FatalError);
			}

			prevVertI = fp;
		}
	}
}


// Create pointFaceIndex_ : map from vertI, faceI to index of vertI on faceI.
void Foam::cellMatcher::calcPointFaceIndex()
{
	// Fill pointFaceIndex_ with -1
	forAll(pointFaceIndex_, i)
	{
		labelList& faceIndices = pointFaceIndex_[i];

		faceIndices = -1;
	}

	forAll(localFaces_, localFaceI)
	{
		const face& f = localFaces_[localFaceI];

		for
		(
			label fp = 0;
			fp < faceSize_[localFaceI];
			fp++
		)
		{
			label vert = f[fp];
			pointFaceIndex_[vert][localFaceI] = fp;
		}
	}
}


// Given edge(v0,v1) and (local)faceI return the other face
Foam::label Foam::cellMatcher::otherFace
(
	const label numVert,
	const label v0,
	const label v1,
	const label localFaceI
) const
{
	label key = edgeKey(numVert, v0, v1);

	if (edgeFaces_[key] == localFaceI)
	{
		return edgeFaces_[key+1];
	}
	else if (edgeFaces_[key+1] == localFaceI)
	{
		return edgeFaces_[key];
	}
	else
	{
		FatalErrorIn
		(
			"otherFace"
			"(const label, const labelList&, const label, const label, "
			"const label)"
		)   << "edgeFaces_ does not contain:" << localFaceI
			<< " for edge " << v0 << " " << v1 << " at key " << key
			<< " edgeFaces_[key, key+1]:" <<  edgeFaces_[key]
			<< " , " << edgeFaces_[key+1]
			<< abort(FatalError);

		return -1;
	}
}


void Foam::cellMatcher::write(Foam::Ostream& os) const
{
	os  << "Faces:" << endl;

	forAll(localFaces_, faceI)
	{
		os  << "    ";

		for(label fp = 0; fp < faceSize_[faceI]; fp++)
		{
			os  << ' ' << localFaces_[faceI][fp];
		}
		os  << endl;
	}

	os  <<  "Face map  : " << faceMap_ << endl;
	os  <<  "Point map : " << pointMap_ << endl;
}


// ************************************************************************* //
