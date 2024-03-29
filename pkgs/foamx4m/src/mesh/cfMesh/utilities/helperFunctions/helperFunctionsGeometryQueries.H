/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
                     Author | F.Juretic (franjo.juretic@c-fields.com)
                  Copyright | Copyright (C) Creative Fields, Ltd.
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam::help

Description
    Geometry queries useful for mesh generation

SourceFiles

\*---------------------------------------------------------------------------*/

#ifndef helperFunctionsGeometryQueries_H
#define helperFunctionsGeometryQueries_H

#include "DynList.H"
#include "plane.H"
#include "face.H"

#include "triSurf.H"
#include "triangle.H"
#include "tetrahedron.H"
#include "boolList.H"
#include "HashSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class boundBox;

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Namespace help functions Declaration
\*---------------------------------------------------------------------------*/

namespace help
{
    //- check if a list has nan entries
    template<class ListType>
    bool isnan(const ListType&);

    //- check if a list has inf entries
    template<class ListType>
    bool isinf(const ListType&);

    //- calculate eigenvalues
    inline vector eigenValues(const tensor&);

    //- calculate eigenvector associated with a given eigenvalue
    inline vector eigenVector(const tensor&, const scalar eigenValue);

    //- check if the faces share a convex edge
    template<class Face1, class Face2>
    inline bool isSharedEdgeConvex
    (
        const pointField& points,
        const Face1& f1,
        const Face2& f2
    );

    //- angle between the two faces in radians
    template<class Face1, class Face2>
    inline scalar angleBetweenFaces
    (
        const pointField& points,
        const Face1& f1,
        const Face2& f2
    );

    // merge faces which are in the same patch
    // This is need to merge the triagles which are generated for cells
    // cut by more than on boundary regions
    inline faceList mergePatchFaces
    (
        const List< DynList<label> >& pfcs,
        const pointField& polyPoints
    );

    //- check if the point p belongs to the edge e
    inline bool vertexOnLine
    (
        const point& p,
        const edge& e,
        const pointField& ep
    );

    //- check if the point p belongs to the plane
    inline bool vertexInPlane(const point& p, const plane& pl);

    //- find the vertex on the line of the edge nearest to the point p
    inline point nearestPointOnTheEdge
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& p
    );

    //- find the vertex on the edge nearest to the point p
    inline point nearestPointOnTheEdgeExact
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& p
    );

    //- find and return the distance between the edge and the point p
    inline scalar distanceOfPointFromTheEdge
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& p
    );

    //- find the nearest points on the edge and the line
    inline bool nearestEdgePointToTheLine
    (
        const point& edgePoint0,
        const point& edgePoint1,
        const point& lp0,
        const point& lp1,
        point& nearestOnEdge,
        point& nearestOnLine
    );

    //- check if the edge intersects the plane
    inline bool planeIntersectsEdge
    (
        const point& start,
        const point& end,
        const plane& pl,
        point& intersection
    );

    //- check if a vertex lies inside the tetrahedron
    inline bool pointInTetrahedron
    (
        const point& p,
        const tetrahedron<point, point>& tet
    );

    //- check if a line intersects the triangle, and return the intersection
    inline bool triLineIntersection
    (
        const triangle<point, point>& tria,
        const point& lineStart,
        const point& lineEnd,
        point& intersection
    );

    //- check if a line intersects the triangle and return the intersection
    inline bool triLineIntersection
    (
        const triSurf&,
        const label,
        const point&,
        const point&,
        point&
    );

    //- check if the line intersects the bounding box
    inline bool boundBoxLineIntersection
    (
        const point&,
        const point&,
        const boundBox&
    );

    //- check if the line and the face intersect
    inline bool lineFaceIntersection
    (
        const point&,
        const point&,
        const face&,
        const pointField& fp,
        point& intersection
    );

    //- check if the surface triangle and the face intersect
    inline bool doFaceAndTriangleIntersect
    (
        const triSurf& surface,
        const label triI,
        const face& f,
        const pointField& facePoints
    );

    //- find the nearest point on the triangle to the given point
    inline point nearestPointOnTheTriangle
    (
        const triangle<point, point>& tri,
        const point&
    );

    //- find the nearest vertex on the surface triangle to the given point
    inline point nearestPointOnTheTriangle
    (
        const label,
        const triSurf&,
        const point&
    );

    //- find the minimiser point from a point and the given planes
    //- returns true if the minimizer exists
    inline bool findMinimizerPoint
    (
        const DynList<point>& origins,
        const DynList<vector>& normals,
        point& pMin
    );

    //- check the existence of overlap between the two edges
    inline bool doEdgesOverlap
    (
        const point& e0p0,
        const point& e0p1,
        const point& e1p0,
        const point& e1p1,
        FixedList<point, 2>& overlappingPart,
        const scalar distTol = -1.0,
        const scalar cosTol = Foam::cos(5.0*(M_PI/180.0)) // cosine tolerance
    );

    //- check the existence of overlap between the two triangles
    inline bool doTrianglesOverlap
    (
        const triangle<point, point>& tri0,
        const triangle<point, point>& tri1,
        DynList<point>& overlappingPolygon,
        const scalar distTol = -1.0,
        const scalar cosTol = Foam::cos(5.0*(M_PI/180.0)) // cosine tolerance
    );

    //- check the existence of intersection between the two triangles
    inline bool doTrianglesIntersect
    (
        const triangle<point, point> &tri0,
        const triangle<point, point> &tri1,
        const scalar distTol = -1.0
    );

    inline bool doTrianglesIntersect
    (
        const triangle<point, point>& tri0,
        const triangle<point, point>& tri1,
        DynList<point>& intersectionPoints,
        const scalar distTol = -1.0
    );

    //- check if the point is inside or outside the face
    inline bool pointInsideFace
    (
        const point& p,
        const face& f,
        const vector& n,
        const pointField& fp,
        const scalar distTol = SMALL
    );

    //- check if the point is inside or outside the face
    inline bool pointInsideFace
    (
        const point& p,
        const face& f,
        const pointField& fp,
        const scalar distTol = SMALL
    );

    //- check if the face is convex. Concave points are flagged false
    inline bool isFaceConvexAndOk
    (
        const face& f,
        const pointField& fp,
        DynList<bool>& OkPoints
    );

    //- calculate quality metric of a tetrahedron
    inline scalar tetQuality(const tetrahedron<point, point>& tet);

    //- check if the vertex is on the positive side of the face plane
    inline bool isVertexVisible(const point& p, const plane& pl);

    //- find number of face groups within a given range
    inline label numberOfFaceGroups
    (
        const labelHashSet& containedElements,
        const point& centre,
        const scalar range,
        const triSurf& surface
    );

    //- find the number of edge groups within the given range
    inline label numberOfEdgeGroups
    (
        const labelHashSet& containedEdges,
        const point& centre,
        const scalar range,
        const triSurf& surface
    );

} // End namespace help

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "helperFunctionsGeometryQueriesI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
