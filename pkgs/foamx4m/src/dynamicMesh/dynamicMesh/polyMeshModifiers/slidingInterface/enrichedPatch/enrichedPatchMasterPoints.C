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
	Enriched patch master points

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.  Copyright Hrvoje Jasak.

\*---------------------------------------------------------------------------*/

#include "enrichedPatch.H"
#include "primitiveMesh.H"
#include "demandDrivenData.H"
#include "DynamicList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::label Foam::enrichedPatch::nFaceHits_ = 4;

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::enrichedPatch::calcMasterPointFaces() const
{
	if (masterPointFacesPtr_)
	{
		FatalErrorIn("void enrichedPatch::calcMasterPointFaces() const")
			<< "Master point face addressing already calculated."
			<< abort(FatalError);
	}

	// Note:
	// Master point face addressing lists the master faces for all points
	// in the enriched patch support (if there are no master faces, which is
	// normal, the list will be empty).  The index represents the index of
	// the master face rather than the index from the enriched patch
	// Master face points lists the points of the enriched master face plus
	// points projected into the master face

	Map<dynamicLabelList > mpf(meshPoints().size());

	const faceList& ef = enrichedFaces();

	// Add the original face points
	forAll (masterPatch_, faceI)
	{
		const face& curFace = ef[faceI + slavePatch_.size()];
//		 Pout << "Cur face in pfAddr: " << curFace << endl;
		forAll (curFace, pointI)
		{
			Map<dynamicLabelList >::iterator mpfIter =
				mpf.find(curFace[pointI]);

			if (mpfIter == mpf.end())
			{
				// Not found, add new dynamic list
				mpf.insert
				(
					curFace[pointI],
					dynamicLabelList(primitiveMesh::facesPerPoint_)
				);

				// Iterator is invalidated - have to find again
				mpf.find(curFace[pointI])().append(faceI);
			}
			else
			{
				mpfIter().append(faceI);
			}
		}
	}

	// Add the projected points which hit the face
	const labelList& slaveMeshPoints = slavePatch_.meshPoints();

	forAll (slavePointFaceHits_, pointI)
	{
		if
		(
			slavePointPointHits_[pointI] < 0
		 && slavePointEdgeHits_[pointI] < 0
		 && slavePointFaceHits_[pointI].hit()
		)
		{
			// Get the index of projected point corresponding to this slave
			// point
			const label mergedSmp =
				pointMergeMap().find(slaveMeshPoints[pointI])();

			Map<dynamicLabelList >::iterator mpfIter =
				mpf.find(mergedSmp);

			if (mpfIter == mpf.end())
			{
				// Not found, add new dynamic list
				mpf.insert
				(
					mergedSmp,
					dynamicLabelList(primitiveMesh::facesPerPoint_)
				);

				// Iterator is invalidated - have to find again
				mpf.find(mergedSmp)().append
				(
					slavePointFaceHits_[pointI].hitObject()
				);
			}
			else
			{
				mpfIter().append(slavePointFaceHits_[pointI].hitObject());
			}
		}
	}

	// Re-pack dynamic lists into normal lists
	const labelList mpfToc = mpf.toc();

	masterPointFacesPtr_ = new Map<labelList>(2*mpfToc.size());
	Map<labelList>& masterPointFaceAddr = *masterPointFacesPtr_;

	forAll (mpfToc, mpfTocI)
	{
		labelList l;
		l.transfer(mpf.find(mpfToc[mpfTocI])().shrink());

		masterPointFaceAddr.insert(mpfToc[mpfTocI], l);
	}
//	 Pout << "masterPointFaceAddr: " << masterPointFaceAddr << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::Map<Foam::labelList>& Foam::enrichedPatch::masterPointFaces() const
{
	if (!masterPointFacesPtr_)
	{
		calcMasterPointFaces();
	}

	return *masterPointFacesPtr_;
}


// ************************************************************************* //
