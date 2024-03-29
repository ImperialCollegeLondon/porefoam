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

#include "IOstreams.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Point, class PointRef, class polygonRef>
inline pyramid<Point, PointRef, polygonRef>::pyramid
(
	polygonRef base,
	const Point& apex
)
:
	base_(base),
	apex_(apex)
{}


template<class Point, class PointRef, class polygonRef>
inline pyramid<Point, PointRef, polygonRef>::pyramid(Istream& is)
{
	is >> base_ >> apex_;
	is.check("pyramid::pyramid(Istream& is)");
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Point, class PointRef, class polygonRef>
inline const Point& pyramid<Point, PointRef, polygonRef>::apex() const
{
	return apex_;
}

template<class Point, class PointRef, class polygonRef>
inline polygonRef pyramid<Point, PointRef, polygonRef>::base() const
{
	return base_;
}


template<class Point, class PointRef, class polygonRef>
inline Point pyramid<Point, PointRef, polygonRef>::centre
(
	const pointField& points
) const
{
	return (3.0/4.0)*base_.centre(points) + (1.0/4.0)*apex_;
}


template<class Point, class PointRef, class polygonRef>
inline vector pyramid<Point, PointRef, polygonRef>::height
(
	const pointField& points
) const
{
	// Height = apex - baseCentroid
	return (apex_ - base_.centre(points));
}


template<class Point, class PointRef, class polygonRef>
inline scalar pyramid<Point, PointRef, polygonRef>::mag
(
	const pointField& points
) const
{
	return (1.0/3.0)*(base_.normal(points)&(height(points)));
}


// * * * * * * * * * * * * * * IOstream Operators  * * * * * * * * * * * * * //

template<class Point, class PointRef, class polygonRef>
inline Istream& operator>>
(
	Istream& is,
	pyramid<Point, PointRef, polygonRef>& p
)
{
	is >> p.base_ >> p.apex_;
	is.check("Istream& operator>>(Istream&, pyramid&)");
	return is;
}


template<class Point, class PointRef, class polygonRef>
inline Ostream& operator<<
(
	Ostream& os,
	const pyramid<Point, PointRef, polygonRef>& p
)
{
	os << p.base_ << tab << p.apex_ << nl;
	return os;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
