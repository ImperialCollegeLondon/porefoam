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
	Foam::triFace

Description
	A triangular face using a FixedList of labels corresponding to mesh
	vertices.

SeeAlso
	Foam::face, Foam::triangle

SourceFiles
	triFaceI.H
	triFaceTemplates.C

\*---------------------------------------------------------------------------*/

#ifndef triFace_H
#define triFace_H

#include "FixedList.H"
#include "edgeList.H"
#include "pointHit.H"
#include "intersection.H"
#include "pointField.H"
#include "triPointRef.H"
#include "ListListOps.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class face;

// Forward declaration of friend functions and operators
class triFace;

inline bool operator==(const triFace&, const triFace&);
inline bool operator!=(const triFace&, const triFace&);



class triFace
:
	public FixedList<label, 3>
{

public:

	// Constructors

		//- Construct null
		inline triFace();

		//- Construct from three point labels
		inline triFace
		(
			const label a,
			const label b,
			const label c
		);

		//- Construct from a UList<label>
		explicit inline triFace(const UList<label>&);

		//- Construct from Istream
		inline triFace(Istream&);


	// Member Functions

		//- Collapse face by removing duplicate point labels
		//  return the collapsed size, set collapsed point labels to -1
		inline label collapse();

		//- Return the edge direction on the face
		//  - +1: forward (counter-clockwise) on the face
		//  - -1: reverse (clockwise) on the face
		//  -  0: edge not found on the face
		inline int edgeDirection(const edge&) const;


		// Properties

			//- Return the points corresponding to this face
			inline pointField points(const pointField& points) const;

			//- Return triangle as a face
			inline face triFaceFace() const;

			//- Return number of edges
			inline label nEdges() const;

			//- Return edges
			inline edgeList edges() const;

			//- Return centre (centroid)
			inline point centre(const pointField&) const;

			//- Return scalar magnitude
			inline scalar mag(const pointField&) const;

			//- Return vector normal
			inline vector normal(const pointField&) const;

			//- Number of triangles after splitting
			inline label nTriangles() const;

			//- Return face with reverse direction
			inline triFace reverseFace() const;

			//- Return swept-volume
			inline scalar sweptVol
			(
				const pointField& oldPoints,
				const pointField& newPoints
			) const;

			//- Return point intersection with a ray starting at p, with
			//  direction n.
			inline pointHit ray
			(
				const point& p,
				const vector& q,
				const pointField& points,
				const intersection::algorithm = intersection::FULL_RAY,
				const intersection::direction dir = intersection::VECTOR
			) const;

			//- Return contact sphere diameter
			inline scalar contactSphereDiameter
			(
				const point& p,
				const vector& n,
				const pointField& meshPoints
			) const;

			//- Return the triangle
			inline triPointRef tri(const pointField&) const;

		//- compare triFaces
		//  -  0: different
		//  - +1: identical
		//  - -1: same face, but different orientation
		static inline int compare(const triFace&, const triFace&);

	// Friend Operators

		friend bool operator==(const triFace&, const triFace&);
		friend bool operator!=(const triFace&, const triFace&);
};


//- Hash specialization for hashing triFace - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<triFace>::operator()
(
	const triFace& t,
	unsigned seed
) const
{
	// Fortunately we don't need this very often
	const uLabel t0(t[0]);
	const uLabel t1(t[1]);
	const uLabel t2(t[2]);

	const uLabel val = (t0*t1*t2 + t0+t1+t2);

	return Hash<uLabel>()(val, seed);
}


//- Hash specialization for hashing triFace - a commutative hash value.
//  Hash incrementally.
template<>
inline unsigned Hash<triFace>::operator()(const triFace& t) const
{
	return Hash<triFace>::operator()(t, 0);
}


template<>
inline bool contiguous<triFace>()  {return true;}


//- Hash specialization to offset faces in ListListOps::combineOffset
template<>
class offsetOp<triFace>
{

public:

	inline triFace operator()
	(
		const triFace& x,
		const label offset
	) const
	{
		triFace result(x);

		forAll(x, xI)
		{
			result[xI] = x[xI] + offset;
		}
		return result;
	}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "triFaceI.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
