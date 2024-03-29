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

#include "error.H"
#include "cellModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline const word& cellModel::name() const
{
	return name_;
}


inline label cellModel::index() const
{
	return index_;
}


inline label cellModel::nPoints() const
{
	return nPoints_;
}


inline label cellModel::nEdges() const
{
	return edges_.size();
}


inline label cellModel::nFaces() const
{
	return faces_.size();
}


//  Return the faces of a cellModel by untangling the geometry
//  supplied in terms of the face labels
inline edgeList cellModel::edges(const labelList& pointLabels) const
{
	edgeList e(edges_.size());

	// Translate model lebels into global labels
	forAll(edges_, edgeI)
	{
		 e[edgeI] =
			 edge
			 (
				 pointLabels[edges_[edgeI].start()],
				 pointLabels[edges_[edgeI].end()]
			 );
	}

	return e;
}


// Return a raw list of model faces
inline const faceList& cellModel::modelFaces() const
{
	return faces_;
}

//  Return the faces of a cellModel by untangling the geometry
//  supplied in terms of the face labels
inline faceList cellModel::faces(const labelList& pointLabels) const
{
	faceList f(faces_.size());

	// Translate model lebels into global labels
	forAll(faces_, faceI)
	{
		 const labelList& curModelLabels = faces_[faceI];

		 face& curFace = f[faceI];

		 curFace.setSize(curModelLabels.size());

		 forAll (curModelLabels, labelI)
		 {
			 curFace[labelI] = pointLabels[curModelLabels[labelI]];
		 }
	}

	return f;
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

// Equality operator: true => ptr to models are equal !
inline bool operator==(const cellModel& m1, const cellModel& m2)
{
	return (&m1 == &m2);
}

// Inequality operator: true => ptr to models are not equal !
inline bool operator!=(const cellModel& m1, const cellModel& m2)
{
	return (&m1 != &m2);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
