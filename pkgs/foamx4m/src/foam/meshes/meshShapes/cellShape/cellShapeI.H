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

#include "Istream.H"
#include "cell.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::cellShape::cellShape()
:
	m_(nullptr)
{}


inline Foam::cellShape::cellShape
(
	const cellModel& M,
	const labelList& l,
	const bool doCollapse
)
:
	labelList(l),
	m_(&M)
{
	if (doCollapse)
	{
		collapse();
	}
}


inline Foam::cellShape::cellShape(Istream& is)
{
	is >> *this;
}


inline Foam::autoPtr<Foam::cellShape> Foam::cellShape::clone() const
{
	return autoPtr<cellShape>(new cellShape(*this));
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::pointField Foam::cellShape::points
(
	const pointField& meshPoints
) const
{
	// There are as many points as there labels for them
	pointField p(size());

	// For each point in list, set it to the point in 'pnts' addressed
	// by 'labs'
	forAll(p, i)
	{
		p[i] = meshPoints[operator[](i)];
	}

	// Return list
	return p;
}


inline bool Foam::cellShape::validModel() const
{
	return m_ != nullptr;
}


inline const Foam::cellModel& Foam::cellShape::model() const
{
	return *m_;
}


inline Foam::labelList Foam::cellShape::meshFaces
(
	const faceList& allFaces,
	const cell& cFaces
) const
{
	// Faces in model order
	faceList localFaces(faces());

	// Do linear match (usually cell shape is low complexity)

	labelList modelToMesh(localFaces.size(), -1);

	forAll(localFaces, i)
	{
		const face& localF = localFaces[i];

		forAll(cFaces, j)
		{
			label meshFaceI = cFaces[j];

			if (allFaces[meshFaceI] == localF)
			{
				modelToMesh[i] = meshFaceI;

				break;
			}
		}
	}

	return modelToMesh;
}


inline Foam::labelList Foam::cellShape::meshEdges
(
	const edgeList& allEdges,
	const labelList& cEdges
) const
{
	// Edges in model order
	edgeList localEdges(edges());

	// Do linear match (usually cell shape is low complexity)

	labelList modelToMesh(localEdges.size(), -1);

	forAll(localEdges, i)
	{
		const edge& e = localEdges[i];

		forAll(cEdges, j)
		{
			label edgeI = cEdges[j];

			if (allEdges[edgeI] == e)
			{
				modelToMesh[i] = edgeI;

				break;
			}
		}
	}

	return modelToMesh;
}


inline Foam::faceList Foam::cellShape::faces() const
{
	return m_->faces(*this);
}


inline Foam::faceList Foam::cellShape::collapsedFaces() const
{
	faceList oldFaces(faces());

	faceList newFaces(oldFaces.size());
	label newFaceI = 0;

	forAll(oldFaces, oldFaceI)
	{
		const face& f = oldFaces[oldFaceI];

		face& newF = newFaces[newFaceI];

		newF.setSize(f.size());

		label newFp = 0;
		label prevVertI = -1;

		forAll(f, fp)
		{
			label vertI = f[fp];

			if (vertI != prevVertI)
			{
				newF[newFp++] = vertI;

				prevVertI = vertI;
			}
		}

		if ((newFp > 1) && (newF[newFp-1] == newF[0]))
		{
			--newFp;
		}

		if (newFp > 2)
		{
			// Size face and go to next one
			newF.setSize(newFp);

			newFaceI++;
		}
	}
	newFaces.setSize(newFaceI);

	return newFaces;
}


inline Foam::label Foam::cellShape::nFaces() const
{
	return m_->nFaces();
}


inline Foam::edgeList Foam::cellShape::edges() const
{
	return m_->edges(*this);
}


inline Foam::label Foam::cellShape::nEdges() const
{
	return m_->nEdges();
}


inline Foam::label Foam::cellShape::nPoints() const
{
	return size();
}


inline Foam::point Foam::cellShape::centre(const pointField& points) const
{
	return m_->centre(*this, points);
}


inline Foam::scalar Foam::cellShape::mag(const pointField& points) const
{
	return m_->mag(*this, points);
}


// ************************************************************************* //
