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
	Foam::triangleFuncs

Description
	Various triangle functions.

SourceFiles
	triangleFuncs.C

\*---------------------------------------------------------------------------*/

#ifndef triangleFuncs_H
#define triangleFuncs_H

#include "point.H"
#include "label.H"
#include "scalar.H"
#include "pointField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class treeBoundBox;


class triangleFuncs
{

public:

	// Public data types

		//- Enumeration defining nearness classification
		enum proxType
		{
			NONE,
			POINT,
			EDGE
		};


private:

	// Private Member Functions

		//- Helper function for intersect. Sets pt to be anywhere on the edge
		//  between oppositeSidePt and thisSidePt depending on both signs.
		static void setIntersection
		(
			const point& oppositeSidePt,
			const scalar oppositeSign,
			const point& thisSidePt,
			const scalar thisSign,
			const scalar tol,
			point& pt
		);

		//- Helper function.
		static void selectPt
		(
			const bool select0,
			const point& p0,
			const point& p1,
			point& min
		);

public:

	//- Intersect triangle with parallel edges aligned with axis i0.
	//  Returns true (and intersection in pInter) if any of them intersects
	//  triangle. Used in intersectBb.
	static bool intersectAxesBundle
	(
		const point& V0,
		const point& V10,
		const point& V20,
		const label i0,
		const pointField& origin,
		const scalar maxLength,
		point& pInter
	);

	//- Does triangle intersect bounding box.
	static bool intersectBb
	(
		const point& p0,
		const point& p1,
		const point& p2,
		const treeBoundBox& cubeBb
	);

	//- Does triangle intersect plane. Return bool and set intersection segment
	static bool intersect
	(
		const point& va0,
		const point& va10,
		const point& va20,

		const point& basePoint,
		const vector& normal,

		point& pInter0,
		point& pInter1
	);

	//- Do triangles intersect. Return bool and set intersection segment
	static bool intersect
	(
		const point& va0,
		const point& va10,
		const point& va20,

		const point& vb0,
		const point& vb10,
		const point& vb20,

		point& pInter0,
		point& pInter1
	);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
