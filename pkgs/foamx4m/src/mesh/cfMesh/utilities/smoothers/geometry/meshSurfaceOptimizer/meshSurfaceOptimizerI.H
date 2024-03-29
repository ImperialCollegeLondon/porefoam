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

Description

\*---------------------------------------------------------------------------*/

#include "demandDrivenData.H"
#include "meshSurfaceOptimizer.H"
#include "plane.H"
#include "Map.H"
#include "surfaceOptimizer.H"
#include "helperFunctions.H"
#include "partTriMeshSimplex.H"

//#define DEBUGSmooth

# ifdef DEBUGSmooth
#include "triSurf.H"
#include "triSurfModifier.H"
#include "helperFunctions.H"
# endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

inline const partTriMesh& meshSurfaceOptimizer::triMesh() const
{
    if( !triMeshPtr_ )
        calculateTrianglesAndAddressing();

    return *triMeshPtr_;
}

inline void meshSurfaceOptimizer::updateTriMesh(const labelLongList& selPoints)
{
    if( !triMeshPtr_ )
        FatalErrorIn
        (
            "inline void meshSurfaceOptimizer::updateTriMesh"
            "(const labelLongList&)"
        ) << "triMeshPtr_ is not allocated " << abort(FatalError);

    triMeshPtr_->updateVertices(selPoints);
}

inline void meshSurfaceOptimizer::updateTriMesh()
{
    if( !triMeshPtr_ )
        FatalErrorIn
        (
            "inline void meshSurfaceOptimizer::updateTriMesh()"
        ) << "triMeshPtr_ is not allocated " << abort(FatalError);

    triMeshPtr_->updateVertices();
}

inline bool meshSurfaceOptimizer::transformIntoPlane
(
    const label bpI,
    const plane& pl,
    vector& vecX,
    vector& vecY,
    DynList<point>& pts,
    DynList<triFace>& trias
) const
{
    # ifdef DEBUGSmooth
    Pout << "Transforming boundary node " << bpI << endl;
    # endif

    const partTriMesh& triMesh = this->triMesh();
    const label triPointI = triMesh.meshSurfacePointLabelInTriMesh()[bpI];
    if( triPointI < 0 )
        return false;

    partTriMeshSimplex simplex(triMesh, triPointI);

    const DynList<point, 32>& sPts = simplex.pts();

    //- create vecX and vecY
    const point& p = simplex.centrePoint();
    pts.setSize(0);
    bool found(false);
    forAll(sPts, pI)
    {
        const point sp = pl.nearestPoint(sPts[pI]);
        const scalar d = mag(sp - p);
        if( d > VSMALL )
        {
            vecX = sp - p;
            vecX /= d;
            vecY = pl.normal() ^ vecX;
            vecY /= mag(vecY);
            found = true;
            break;
        }
    }

    if( !found )
        return false;

    trias = simplex.triangles();

    # ifdef DEBUGSmooth
    Pout << "VecX " << vecX << endl;
    Pout << "vecY " << vecY << endl;
    # endif

    //- transform the vertices
    pts.setSize(sPts.size());

    forAll(sPts, pI)
    {
        const point pt
        (
            ((sPts[pI] - p) & vecX),
            ((sPts[pI] - p) & vecY),
            0.0
        );

        pts[pI] = pt;
    }

    # ifdef DEBUGSmooth
    Info << "Original triangles " << endl;
    forAll(simplex.triangles(), triI)
        Info << "Tri " << triI << " is " << simplex.triangles()[triI] << endl;
    Info << "Transformed triangles are " << trias << endl;
    Info << "Transformed vertices " << pts << endl;

    triSurf surf;
    triSurfModifier sMod(surf);
    pointField& sPoints = sMod.pointsAccess();
    sPoints.setSize(pts.size());
    forAll(sPoints, i)
        sPoints[i] = pts[i];
    LongList<labelledTri>& sTrias = sMod.facetsAccess();
    sTrias.setSize(trias.size());
    forAll(sTrias, i)
    {
        labelledTri& tf = sTrias[i];
        tf[0] = trias[i][0];
        tf[1] = trias[i][1];
        tf[2] = trias[i][2];

        tf.region() = 0;
    }
    sMod.patchesAccess().setSize(1);
    sMod.patchesAccess()[0].name() = "bnd";
    sMod.patchesAccess()[0].geometricType() = "patch";

    fileName sName("bndSimplex_");
    sName += help::scalarToText(bpI);
    sName += ".stl";
    surf.writeSurface(sName);
    # endif

    return found;
}

inline point meshSurfaceOptimizer::newPositionLaplacian
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const VRWGraph& pPoints = surfaceEngine_.pointPoints();
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    vector newP(vector::zero);
    if( transformIntoPlane )
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if( magSqr(pNormal) < VSMALL )
            return points[bPoints[bpI]];

        plane pl(points[bPoints[bpI]], pNormal);

        DynList<point> projectedPoints;
        projectedPoints.setSize(pPoints.sizeOfRow(bpI));
        forAllRow(pPoints, bpI, pI)
            projectedPoints[pI] =
                pl.nearestPoint(points[bPoints[pPoints(bpI, pI)]]);

        forAll(projectedPoints, pI)
            newP += projectedPoints[pI];

        newP /= projectedPoints.size();
    }
    else
    {
        forAllRow(pPoints, bpI, pI)
            newP += points[bPoints[pPoints(bpI, pI)]];

        newP /= pPoints.sizeOfRow(bpI);
    }

    return newP;
}

inline point meshSurfaceOptimizer::newPositionLaplacianFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    vector newP(vector::zero);

    if( transformIntoPlane )
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if( magSqr(pNormal) < VSMALL )
            return points[bPoints[bpI]];

        plane pl(points[bPoints[bpI]], pNormal);

        DynList<point> projectedPoints;
        projectedPoints.setSize(pointFaces.sizeOfRow(bpI));
        forAllRow(pointFaces, bpI, pfI)
            projectedPoints[pfI] =
                pl.nearestPoint(faceCentres[pointFaces(bpI, pfI)]);

        forAll(projectedPoints, pI)
            newP += projectedPoints[pI];

        newP /= projectedPoints.size();
    }
    else
    {
        forAllRow(pointFaces, bpI, pfI)
            newP += faceCentres[pointFaces(bpI, pfI)];

        newP /= pointFaces.sizeOfRow(bpI);
    }

    return newP;
}

inline point meshSurfaceOptimizer::newPositionLaplacianWFC
(
    const label bpI,
    const bool transformIntoPlane
) const
{
    const VRWGraph& pointFaces = surfaceEngine_.pointFaces();
    const pointFieldPMG& points = surfaceEngine_.points();
    const vectorField& faceAreas = surfaceEngine_.faceNormals();
    const vectorField& faceCentres = surfaceEngine_.faceCentres();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    vector newP(vector::zero);

    if( transformIntoPlane )
    {
        const vector& pNormal = surfaceEngine_.pointNormals()[bpI];

        if( magSqr(pNormal) < VSMALL )
            return points[bPoints[bpI]];

        plane pl(points[bPoints[bpI]], pNormal);

        DynList<point> projectedPoints;
        projectedPoints.setSize(pointFaces.sizeOfRow(bpI));
        forAllRow(pointFaces, bpI, pfI)
            projectedPoints[pfI] =
                pl.nearestPoint(faceCentres[pointFaces(bpI, pfI)]);

        scalar sumW(0.0);
        forAll(projectedPoints, pI)
        {
            const label bfI = pointFaces(bpI, pI);
            const scalar w = (faceAreas[bfI] & faceAreas[bfI]);
            newP += w * projectedPoints[pI];
            sumW += w;
        }

        newP /= Foam::max(sumW, VSMALL);
    }
    else
    {
        scalar sumW(0.0);
        forAllRow(pointFaces, bpI, pfI)
        {
            const label bfI = pointFaces(bpI, pfI);
            const scalar w = (faceAreas[bfI] & faceAreas[bfI]);
            newP += w * faceCentres[pointFaces(bpI, pfI)];
            sumW += w;
        }

        newP /= Foam::max(sumW, VSMALL);
    }

    return newP;
}

inline point meshSurfaceOptimizer::newPositionSurfaceOptimizer
(
    const label bpI,
    const scalar tol
) const
{
    const pointFieldPMG& points = surfaceEngine_.points();
    const labelList& bPoints = surfaceEngine_.boundaryPoints();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    # ifdef DEBUGSmooth
    Pout << "Smoothing boundary node " << bpI << endl;
    Pout << "Node label in the mesh is " << bPoints[bpI] << endl;
    Pout << "Point coordinates " << points[bPoints[bpI]] << endl;
    # endif

    //- project vertices onto the plane
    const vector& pNormal = surfaceEngine_.pointNormals()[bpI];
    if( magSqr(pNormal) < VSMALL )
        return points[bPoints[bpI]];

    const plane pl(points[bPoints[bpI]], pNormal);

    DynList<point> pts;
    DynList<triFace> trias;
    vector vecX, vecY;
    bool success = this->transformIntoPlane(bpI, pl, vecX, vecY, pts, trias);
    if( !success )
    {
        //Warning << "Cannot transform to plane" << endl;
        return points[bPoints[bpI]];
    }

    surfaceOptimizer so(pts, trias);
    const point newPoint = so.optimizePoint(tol);

    const point newP
    (
        points[bPoints[bpI]] +
        vecX * newPoint.x() +
        vecY * newPoint.y()
    );

    if( help::isnan(newP) || help::isinf(newP) )
    {
        WarningIn
        (
            "inline point meshSurfaceOptimizer::newPositionSurfaceOptimizer"
            "(const label, const scalar) const"
        ) << "Cannot move point " << bpI << endl;

        return points[bPoints[bpI]];
    }

    return newP;
}

inline point meshSurfaceOptimizer::newEdgePositionLaplacian
(
    const label bpI
) const
{
    const labelList& bPoints = surfaceEngine_.boundaryPoints();
    const edgeList& edges = surfaceEngine_.edges();
    const VRWGraph& bpEdges = surfaceEngine_.boundaryPointEdges();
    const pointFieldPMG& points = surfaceEngine_.points();

    if( vertexType_[bpI] & LOCKED )
        return points[bPoints[bpI]];

    const labelHashSet& featureEdges = partitionerPtr_->featureEdges();

    DynList<label> edgePoints;

    forAllRow(bpEdges, bpI, i)
    {
        const label beI = bpEdges(bpI, i);

        if( featureEdges.found(beI) )
        {
            edgePoints.append(edges[beI].otherVertex(bPoints[bpI]));
        }
    }

    if( edgePoints.size() != 2 )
        return points[bPoints[bpI]];

    # ifdef DEBUGSearch
    Info << "Edge points " << edgePoints << endl;
    # endif

    vector pos(vector::zero);
    forAll(edgePoints, epI)
        pos += points[edgePoints[epI]];

    pos /= edgePoints.size();

    return pos;
}

template<class labelListType>
inline void meshSurfaceOptimizer::lockBoundaryFaces(const labelListType& l)
{
    lockedSurfaceFaces_ = l;

    const faceList::subList& bFaces = surfaceEngine_.boundaryFaces();
    const labelList& bp = surfaceEngine_.bp();

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 20)
    # endif
    forAll(lockedSurfaceFaces_, lfI)
    {
        const face& bf = bFaces[lockedSurfaceFaces_[lfI]];

        forAll(bf, pI)
            vertexType_[bp[bf[pI]]] |= LOCKED;
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const DynList<label>& bpNeiProcs = surfaceEngine_.bpNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData[bpNeiProcs[i]].clear();

        //- prepare data which will be sent to other processors
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( vertexType_[bpI] & LOCKED )
            {
                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(it.key());
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i]];

            vertexType_[bpI] |= LOCKED;
        }
    }
}

template<class labelListType>
inline void meshSurfaceOptimizer::lockBoundaryPoints(const labelListType& l)
{
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(l, i)
    {
        const label bpI = l[i];

        vertexType_[bpI] |= LOCKED;
    }

    if( Pstream::parRun() )
    {
        const Map<label>& globalToLocal =
            surfaceEngine_.globalToLocalBndPointAddressing();
        const VRWGraph& bpAtProcs = surfaceEngine_.bpAtProcs();
        const DynList<label>& bpNeiProcs = surfaceEngine_.bpNeiProcs();

        std::map<label, labelLongList> exchangeData;
        forAll(bpNeiProcs, i)
            exchangeData[bpNeiProcs[i]].clear();

        //- prepare data which will be sent to other processors
        forAllConstIter(Map<label>, globalToLocal, it)
        {
            const label bpI = it();

            if( vertexType_[bpI] & LOCKED )
            {
                forAllRow(bpAtProcs, bpI, i)
                {
                    const label neiProc = bpAtProcs(bpI, i);

                    if( neiProc == Pstream::myProcNo() )
                        continue;

                    exchangeData[neiProc].append(it.key());
                }
            }
        }

        labelLongList receivedData;
        help::exchangeMap(exchangeData, receivedData);

        forAll(receivedData, i)
        {
            const label bpI = globalToLocal[receivedData[i]];

            vertexType_[bpI] |= LOCKED;
        }
    }
}

//- lock edge points
inline void meshSurfaceOptimizer::lockFeatureEdges()
{
    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(vertexType_, bpI)
        if( vertexType_[bpI] & (EDGE | CORNER) )
            vertexType_[bpI] |= LOCKED;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
