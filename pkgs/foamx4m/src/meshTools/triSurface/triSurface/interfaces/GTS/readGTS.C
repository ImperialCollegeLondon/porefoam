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

#include "triSurface.H"
#include "IFstream.H"
#include "IStringStream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool triSurface::readGTS(const fileName& GTSfileName)
{
	IFstream GTSfile(GTSfileName);

	if (!GTSfile.good())
	{
		FatalErrorIn("triSurface::readGTS(const fileName&)")
			<< "Cannot read file " << GTSfileName
			<< exit(FatalError);
	}

	// Read header
	label nPoints, nEdges, nElems;

	string line = getLineNoComment(GTSfile);
	{
		IStringStream lineStream(line);
		lineStream >> nPoints >> nEdges >> nElems;
	}

	// Read points
	pointField& points_ = const_cast<pointField&>(points());
	points_.setSize(nPoints);

	forAll(points_, pointi)
	{
		scalar x, y, z;
		line = getLineNoComment(GTSfile);
		{
			IStringStream lineStream(line);
			lineStream >> x >> y >> z;
		}
		points_[pointi] = point(x, y, z);
	}

	// Read edges (Foam indexing)
	edgeList edges(nEdges);
	forAll(edges, edgei)
	{
		label start, end;
		line = getLineNoComment(GTSfile);
		{
			IStringStream lineStream(line);
			lineStream >> start >> end;
		}
		edges[edgei] = edge(start - 1, end - 1);
	}

	// Read triangles. Convert references to edges into pointlabels
	setSize(nElems);
	forAll(*this, trianglei)
	{
		label e0Label, e1Label, e2Label;
		label region = 0;

		line = getLineNoComment(GTSfile);
		{
			IStringStream lineStream(line);
			lineStream >> e0Label >> e1Label >> e2Label;

			// Optional region number: read first, then check state on stream
			if (lineStream)
			{
				label num;
				lineStream >> num;
				if (!lineStream.bad())
				{
					region = num;
				}
			}
		}

		// Determine ordering of edges e0, e1
		//  common:common vertex, shared by e0 and e1
		//  e0Far:vertex on e0 which is not common
		//  e1Far: ,,       e1  ,,
		const edge& e0 = edges[e0Label - 1];
		const edge& e1 = edges[e1Label - 1];
		const edge& e2 = edges[e2Label - 1];

		label common01 = e0.commonVertex(e1);
		if (common01 == -1)
		{
			FatalErrorIn("triSurface::readGTS(const fileName&)")
				<< "Edges 0 and 1 of triangle " << trianglei
				<< " do not share a point.\n"
				<< "    edge0:" << e0 << endl
				<< "    edge1:" << e1
				<< exit(FatalError);
		}

		label e0Far = e0.otherVertex(common01);
		label e1Far = e1.otherVertex(common01);

		label common12 = e1.commonVertex(e2);
		if (common12 == -1)
		{
			FatalErrorIn("triSurface::readGTS(const fileName&)")
				<< "Edges 1 and 2 of triangle " << trianglei
				<< " do not share a point.\n"
				<< "    edge1:" << e1 << endl
				<< "    edge2:" << e2
				<< exit(FatalError);
		}
		label e2Far = e2.otherVertex(common12);

		// Does edge2 sit between edge1 and 0?
		if ((common12 != e1Far) || (e2Far != e0Far))
		{
			FatalErrorIn("triSurface::readGTS(const fileName&)")
				<< "Edges of triangle " << trianglei
				<< " reference more than three points.\n"
				<< "    edge0:" << e0 << endl
				<< "    edge1:" << e1 << endl
				<< "    edge2:" << e2 << endl
				<< exit(FatalError);
		}

		operator[](trianglei) = labelledTri(e0Far, common01, e1Far, region);
	}

	// Construct patch names
	setDefaultPatches();

	return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
