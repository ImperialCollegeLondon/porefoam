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
	Foam::triSurfaceSearch

Description
	Helper class to search on triSurface.

SourceFiles
	triSurfaceSearch.C

\*---------------------------------------------------------------------------*/

#ifndef triSurfaceSearch_H
#define triSurfaceSearch_H

#include "pointField.H"
#include "boolList.H"
#include "pointIndexHit.H"
#include "indexedOctree.H"
#include "treeDataTriSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class triSurface;
class treeDataTriSurface;
template<class Type> class indexedOctree;


class triSurfaceSearch
{
	// Private data

		//- Reference to surface to work on
		const triSurface& surface_;

		//- Octree for searches
		autoPtr<indexedOctree<treeDataTriSurface> > treePtr_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		triSurfaceSearch(const triSurfaceSearch&);

		//- Disallow default bitwise assignment
		void operator=(const triSurfaceSearch&);


public:

	// Static data members

		//- Point far away; used for illegal finds
		static const point greatPoint;


	// Constructors

		//- Construct from surface. Holds reference to surface!
		explicit triSurfaceSearch(const triSurface&);


	// Member Functions

		const indexedOctree<treeDataTriSurface>& tree() const
		{
			return treePtr_();
		}

		const triSurface& surface() const
		{
			return surface_;
		}

		//- Calculate for each searchPoint inside/outside status.
		bool calcInside(const point& p) const;

		//- Calculate for each searchPoint inside/outside status.
		boolList calcInside(const pointField& searchPoints) const;

		//- Calculate index of nearest triangle (or -1) for each sample.
		//  Looks only in box of size 2*span around sample.
		labelList calcNearestTri
		(
			const pointField& samples,
			const vector& span
		) const;

		//- Calculate nearest points (to searchPoints) on surface.
		//  Looks only in box of size 2*span around sample. Returns greatPoint
		//  if not found.
		tmp<pointField> calcNearest
		(
			const pointField& samples,
			const vector& span
		) const;

		//- Calculate nearest point on surface for single searchPoint. Returns
		//  in pointIndexHit:
		//  - hit()      : whether nearest point found within bounding box
		//  - hitPoint() : coordinate of nearest point
		//  - index()    : surface triangle label
		pointIndexHit nearest(const point&, const vector& span) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
