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

Class
	Foam::edgeMesh

Description
	points connected by edges.

SourceFiles
	edgeMeshI.H
	edgeMesh.C
	edgeMeshIO.C

\*---------------------------------------------------------------------------*/

#ifndef edgeMesh_H
#define edgeMesh_H

#include "pointField.H"
#include "edgeList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class edgeMesh;

Ostream& operator<<(Ostream&, const edgeMesh&);
Istream& operator>>(Istream&, edgeMesh&);



class edgeMesh
{
	// Private data

		//- Vertices of the edges
		pointField points_;

		//- The edges defining the boundary
		edgeList edges_;

		//- From point to edges
		mutable autoPtr<labelListList> pointEdgesPtr_;

	// Private Member Functions

		//- Calculate point-edge addressing (inverse of edges)
		void calcPointEdges() const;

public:

	// Constructors

		//- Construct from components
		edgeMesh(const pointField&, const edgeList&);

		//- Construct from file
		edgeMesh(const fileName&);

		//- Construct from Istream
		edgeMesh(Istream&);

		//- Construct as copy
		edgeMesh(const edgeMesh&);



	// Member Functions

		inline const pointField& points() const;

		inline const edgeList& edges() const;

		inline const labelListList& pointEdges() const;

		//- Find connected regions. Set region number per edge.
		//  Returns number of regions.
		label regions(labelList& edgeRegion) const;

		//- Merge common points (points within mergeDist)
		void mergePoints(const scalar mergeDist);

	// Member Operators

		inline void operator=(const edgeMesh&);

	// Ostream Operator

		friend Ostream& operator<<(Ostream&, const edgeMesh&);
		friend Istream& operator>>(Istream&, edgeMesh&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "edgeMeshI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
