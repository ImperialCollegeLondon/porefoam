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
	Foam::searchableCylinder

Description
	Searching on cylinder

SourceFiles
	searchableCylinder.C

\*---------------------------------------------------------------------------*/

#ifndef searchableCylinder_H
#define searchableCylinder_H

#include "treeBoundBox.H"
#include "searchableSurface.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes


class searchableCylinder
:
	public searchableSurface
{
private:

	// Private Member Data

		//- 'left' point
		const point point1_;

		//- 'right' point
		const point point2_;

		//- length of vector point2-point1
		const scalar magDir_;

		//- normalised vector point2-point1
		const vector unitDir_;

		//- Radius squared
		const scalar radius_;

		//- Names of regions
		mutable wordList regions_;


	// Private Member Functions

		//- Find nearest point on cylinder.
		pointIndexHit findNearest
		(
			const point& sample,
			const scalar nearestDistSqr
		) const;

		scalar radius2(const point& pt) const;

		//- Find intersection with cylinder
		void findLineAll
		(
			const point& start,
			const point& end,
			pointIndexHit& near,
			pointIndexHit& far
		) const;


		//- Disallow default bitwise copy construct
		searchableCylinder(const searchableCylinder&);

		//- Disallow default bitwise assignment
		void operator=(const searchableCylinder&);


public:

	//- Runtime type information
	TypeName("searchableCylinder");


	// Constructors

		//- Construct from components
		searchableCylinder
		(
			const IOobject& io,
			const point&,
			const point&,
			const scalar radius
		);

		//- Construct from dictionary (used by searchableSurface)
		searchableCylinder
		(
			const IOobject& io,
			const dictionary& dict
		);

	// Destructor

		virtual ~searchableCylinder();


	// Member Functions

		virtual const wordList& regions() const;

		//- Whether supports volume type below
		virtual bool hasVolumeType() const
		{
			return true;
		}

		//- Range of local indices that can be returned.
		virtual label size() const
		{
			return 1;
		}

		//- Get representative set of element coordinates
		//  Usually the element centres (should be of length size()).
		virtual pointField coordinates() const;


		// Multiple point queries.

			virtual void findNearest
			(
				const pointField& sample,
				const scalarField& nearestDistSqr,
				List<pointIndexHit>&
			) const;

			virtual void findLine
			(
				const pointField& start,
				const pointField& end,
				List<pointIndexHit>&
			) const;

			virtual void findLineAny
			(
				const pointField& start,
				const pointField& end,
				List<pointIndexHit>&
			) const;

			//- Get all intersections in order from start to end.
			virtual void findLineAll
			(
				const pointField& start,
				const pointField& end,
				List<List<pointIndexHit> >&
			) const;

			//- From a set of points and indices get the region
			virtual void getRegion
			(
				const List<pointIndexHit>&,
				labelList& region
			) const;

			//- From a set of points and indices get the normal
			virtual void getNormal
			(
				const List<pointIndexHit>&,
				vectorField& normal
			) const;

			//- Determine type (inside/outside/mixed) for point. unknown if
			//  cannot be determined (e.g. non-manifold surface)
			virtual void getVolumeType
			(
				const pointField&,
				List<volumeType>&
			) const;


		// regIOobject implementation

			bool writeData(Ostream&) const
			{
				notImplemented("searchableCylinder::writeData(Ostream&) const");
				return false;
			}

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
