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

#include "IOstreams.H"
#include "face.H"
#include "triPointRef.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

inline int Foam::triFace::compare(const triFace& a, const triFace& b)
{
	if
	(
		(a[0] == b[0] && a[1] == b[1] && a[2] == b[2])
	 || (a[0] == b[1] && a[1] == b[2] && a[2] == b[0])
	 || (a[0] == b[2] && a[1] == b[0] && a[2] == b[1])
	)
	{
		// identical
		return 1;
	}
	else if
	(
		(a[0] == b[2] && a[1] == b[1] && a[2] == b[0])
	 || (a[0] == b[1] && a[1] == b[0] && a[2] == b[2])
	 || (a[0] == b[0] && a[1] == b[2] && a[2] == b[1])
	)
	{
		// same face, but reversed orientation
		return -1;
	}
	else
	{
		return 0;
	}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

inline Foam::triFace::triFace()
{}


inline Foam::triFace::triFace
(
	const label a,
	const label b,
	const label c
)
{
	operator[](0) = a;
	operator[](1) = b;
	operator[](2) = c;
}


inline Foam::triFace::triFace(const UList<label>& lst)
:
	FixedList<label, 3>(lst)
{}


inline Foam::triFace::triFace(Istream& is)
:
	FixedList<label, 3>(is)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline Foam::label Foam::triFace::collapse()
{
	// we cannot resize a FixedList, so mark duplicates with '-1'
	// (the lower vertex is retained)
	// catch any '-1' (eg, if called twice)

	label n = 3;
	if (operator[](0) == operator[](1) || operator[](1) == -1)
	{
		operator[](1) = -1;
		n--;
	}
	else if (operator[](1) == operator[](2) || operator[](2) == -1)
	{
		operator[](2) = -1;
		n--;
	}
	if (operator[](0) == operator[](2))
	{
		operator[](2) = -1;
		n--;
	}

	return n;
}


inline Foam::pointField Foam::triFace::points(const pointField& points) const
{
	pointField p(3);

	p[0] = points[operator[](0)];
	p[1] = points[operator[](1)];
	p[2] = points[operator[](2)];

	return p;
}


inline Foam::face Foam::triFace::triFaceFace() const
{
	face f(3);

	f[0] = operator[](0);
	f[1] = operator[](1);
	f[2] = operator[](2);

	return f;
}


inline Foam::label Foam::triFace::nEdges() const
{
	return 3;
}


inline Foam::edgeList Foam::triFace::edges() const
{
	edgeList e(3);

	e[0].start() = operator[](0);
	e[0].end() = operator[](1);

	e[1].start() = operator[](1);
	e[1].end() = operator[](2);

	e[2].start() = operator[](2);
	e[2].end() = operator[](0);

	return e;
}


// return
//  - +1: forward (counter-clockwise) on the face
//  - -1: reverse (clockwise) on the face
//  -  0: edge not found on the face
inline int Foam::triFace::edgeDirection(const edge& e) const
{
	if
	(
		(operator[](0) == e.start() && operator[](1) == e.end())
	 || (operator[](1) == e.start() && operator[](2) == e.end())
	 || (operator[](2) == e.start() && operator[](0) == e.end())
	)
	{
		return 1;
	}
	else if
	(
		(operator[](0) == e.end() && operator[](1) == e.start())
	 || (operator[](1) == e.end() && operator[](2) == e.start())
	 || (operator[](2) == e.end() && operator[](0) == e.start())
	)
	{
		return -1;
	}
	else
	{
		return 0;
	}
}


inline Foam::point Foam::triFace::centre(const pointField& points) const
{
	return (1.0/3.0)*
	(
		points[operator[](0)]
	  + points[operator[](1)]
	  + points[operator[](2)]
	);
}


inline Foam::scalar Foam::triFace::mag(const pointField& points) const
{
	return ::Foam::mag(normal(points));
}

// could also delegate to triPointRef(...).normal()
inline Foam::vector Foam::triFace::normal(const pointField& points) const
{
	return 0.5*
	(
		(points[operator[](1)] - points[operator[](0)])
	   ^(points[operator[](2)] - points[operator[](0)])
	);
}


inline Foam::label Foam::triFace::nTriangles() const
{
	return 1;
}


inline Foam::triFace Foam::triFace::reverseFace() const
{
	// The starting points of the original and reverse face are identical.
	return triFace(operator[](0), operator[](2), operator[](1));
}


inline Foam::scalar Foam::triFace::sweptVol
(
	const pointField& opts,
	const pointField& npts
) const
{
	return (1.0/6.0)*
	(
		(
			(npts[operator[](0)] - opts[operator[](0)])
		  & (
				(opts[operator[](1)] - opts[operator[](0)])
			  ^ (opts[operator[](2)] - opts[operator[](0)])
			)
		)
	  + (
			(npts[operator[](1)] - opts[operator[](1)])
		  & (
				(opts[operator[](2)] - opts[operator[](1)])
			  ^ (npts[operator[](0)] - opts[operator[](1)])
			)
		)
	  + (
			(opts[operator[](2)] - npts[operator[](2)])
		  & (
				(npts[operator[](1)] - npts[operator[](2)])
			  ^ (npts[operator[](0)] - npts[operator[](2)])
			)
		)
	);
}


inline Foam::pointHit Foam::triFace::ray
(
	const point& p,
	const vector& q,
	const pointField& points,
	const intersection::algorithm alg,
	const intersection::direction dir
) const
{
	return triPointRef
	(
		points[operator[](0)],
		points[operator[](1)],
		points[operator[](2)]
	).ray(p, q, alg, dir);
}


inline Foam::scalar Foam::triFace::contactSphereDiameter
(
	const point& p,
	const vector& n,
	const pointField& meshPoints
) const
{
	vector n1 = n/(Foam::mag(n) + SMALL);
	vector n2 = normal(meshPoints);

	n2 /= Foam::mag(n2) + SMALL;

	return 2*((centre(meshPoints) - p) & n2)/((n1 & n2) - 1.0);
}


inline Foam::triPointRef Foam::triFace::tri(const pointField& points) const
{
	return triPointRef
	(
		points[operator[](0)],
		points[operator[](1)],
		points[operator[](2)]
	);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

inline bool Foam::operator==(const triFace& a, const triFace& b)
{
	return triFace::compare(a,b) != 0;
}


inline bool Foam::operator!=(const triFace& a, const triFace& b)
{
	return triFace::compare(a,b) == 0;
}


// ************************************************************************* //
