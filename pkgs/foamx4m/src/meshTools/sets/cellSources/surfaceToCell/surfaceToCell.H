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
	Foam::surfaceToCell

Description
	A topoSetSource to select cells based on relation to surface.

	Selects:
	- all cells inside/outside/cut by surface
	- cells with centre nearer than XXX to surface
	- cells with centre nearer than XXX to surface @b and with normal
	  at nearest point to centre and cell-corners differing by
	  more than YYY (i.e., point of high curvature)

SourceFiles
	surfaceToCell.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceToCell_H
#define surfaceToCell_H

#include "topoSetSource.H"
#include "Map.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
class triSurfaceSearch;
class triSurface;


class surfaceToCell
:
	public topoSetSource
{

	// Private data

		//- Add usage string
		static addToUsageTable usage_;

		//- Name of surface file
		fileName surfName_;

		//- Points which are outside
		pointField outsidePoints_;

		//- Include cut cells
		bool includeCut_;

		//- Include inside cells
		bool includeInside_;

		//- Include outside cells
		bool includeOutside_;

		//- if > 0 : include cells with distance from cellCentre to surface
		//  less than nearDist.
		scalar nearDist_;

		//- if > -1 : include cells with normals at nearest surface points
		//  varying more than curvature_.
		scalar curvature_;

		//- triSurface to search on. On pointer since can be external.
		const triSurface* surfPtr_;

		//- search engine on surface.
		const triSurfaceSearch* querySurfPtr_;

		//- whether I allocated above surface ptrs or whether they are
		//  external.
		bool IOwnPtrs_;


	// Private Member Functions

		//- Find index of nearest triangle to point. Returns triangle or -1 if
		//  not found within search span.
		//  Cache result under pointI.
		static label getNearest
		(
			const triSurfaceSearch& querySurf,
			const label pointI,
			const point& pt,
			const vector& searchSpan,
			Map<label>& cache
		);

		//- Return true if surface normal of nearest points to vertices on
		//  cell differ from that on cell centre. Points cached in
		//  pointToNearest.
		bool differingPointNormals
		(
			const triSurfaceSearch& querySurf,
			const vector& span,
			const label cellI,
			const label cellTriI,
			Map<label>& pointToNearest
		) const;


		//- Depending on surface add to or delete from cellSet.
		void combine(topoSet& set, const bool add) const;

		//- Check values at construction time.
		void checkSettings() const;

		const triSurfaceSearch& querySurf() const
		{
			return *querySurfPtr_;
		}


public:

	//- Runtime type information
	TypeName("surfaceToCell");

	// Constructors

		//- Construct from components
		surfaceToCell
		(
			const polyMesh& mesh,
			const fileName& surfName,
			const pointField& outsidePoints,
			const bool includeCut,
			const bool includeInside,
			const bool includeOutside,
			const scalar nearDist,
			const scalar curvature
		);

		//- Construct from components (supplied surface, surfaceSearch)
		surfaceToCell
		(
			const polyMesh& mesh,
			const fileName& surfName,
			const triSurface& surf,
			const triSurfaceSearch& querySurf,
			const pointField& outsidePoints,
			const bool includeCut,
			const bool includeInside,
			const bool includeOutside,
			const scalar nearDist,
			const scalar curvature
		);

		//- Construct from dictionary
		surfaceToCell
		(
			const polyMesh& mesh,
			const dictionary& dict
		);

		//- Construct from Istream
		surfaceToCell
		(
			const polyMesh& mesh,
			Istream&
		);


	// Destructor

		virtual ~surfaceToCell();


	// Member Functions

		virtual void applyToSet
		(
			const topoSetSource::setAction action,
			topoSet&
		) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
