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
	Foam::pyramid

Description
	A geometric pyramid primitive with a base of 'n' sides:
	i.e. a parametric pyramid. A pyramid is constructed from
	a base polygon and an apex point.

SourceFiles
	pyramidI.H

\*---------------------------------------------------------------------------*/

#ifndef pyramid_H
#define pyramid_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Point, class PointRef, class polygonRef>
class pyramid;

template<class Point, class PointRef, class polygonRef>
inline Istream& operator>>
(
	Istream&,
	pyramid<Point, PointRef, polygonRef>&
);

template<class Point, class PointRef, class polygonRef>
inline Ostream& operator<<
(
	Ostream&,
	const pyramid<Point, PointRef, polygonRef>&
);



template<class Point, class PointRef, class polygonRef>
class pyramid
{
	// Private data

		polygonRef base_;
		PointRef apex_;


public:

	// Constructors

		//- Construct from base polygon and apex point
		inline pyramid(polygonRef base, const Point& apex);

		//- Construct from Istream
		inline pyramid(Istream&);


	// Member functions

		// Access

			//- Return apex point
			inline const Point& apex() const;

			//- Return base polygon
			inline polygonRef base() const;


		// Properties

			//- Return centre (centroid)
			inline Point centre(const pointField& points) const;

			//- Return height vector
			inline vector height(const pointField& points) const;

			//- Return scalar magnitude - returns volume of pyramid
			inline scalar mag(const pointField& points) const;


	// IOstream operators

		friend Istream& operator>> <Point, PointRef, polygonRef>
		(
			Istream&,
			pyramid&
		);

		friend Ostream& operator<< <Point, PointRef, polygonRef>
		(
			Ostream&,
			const pyramid&
		);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "pyramidI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
