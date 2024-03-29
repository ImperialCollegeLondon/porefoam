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

\*---------------------------------------------------------------------------*/

#include "edgeMeshGeometryModification.H"
#include "demandDrivenData.H"
#include "dictionary.H"
#include "edgeMesh.H"

namespace Foam
{

// * * * * * * * * * * * * * * Private member functions* * * * * * * * * * * //

void edgeMeshGeometryModification::checkModification()
{
    if( meshDict_.found("anisotropicSources") )
    {
        modificationActive_ = true;

        const dictionary& anisotropicDict =
            meshDict_.subDict("anisotropicSources");

        coordinateModifierPtr_ = new coordinateModifier(anisotropicDict);
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

edgeMeshGeometryModification::edgeMeshGeometryModification
(
    const edgeMesh& em,
    const dictionary& meshDict
)
:
    edgeMesh_(em),
    meshDict_(meshDict),
    coordinateModifierPtr_(nullptr),
    modificationActive_(false)
{
    checkModification();
}

edgeMeshGeometryModification::~edgeMeshGeometryModification()
{
    deleteDemandDrivenData(coordinateModifierPtr_);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool edgeMeshGeometryModification::activeModification() const
{
    return modificationActive_;
}

const edgeMesh* edgeMeshGeometryModification::modifyGeometry() const
{
    if( !modificationActive_ )
    {
        WarningIn
        (
            "const edgeMesh* edgeMeshGeometryModification"
            "::modifyGeometry() const"
        ) << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = edgeMesh_.points();

    pointField newPts(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
        newPts[pointI] = coordinateModifierPtr_->modifiedPoint(pts[pointI]);

    const edgeMesh* newEdgeMesh = new edgeMesh(newPts, edgeMesh_.edges());

    return newEdgeMesh;
}

const edgeMesh* edgeMeshGeometryModification::
revertGeometryModification() const
{
    if( !modificationActive_ )
    {
        WarningIn
        (
            "const edgeMesh* edgeMeshGeometryModification"
            "::revertGeometryModification() const"
        ) << "Modification is not active" << endl;

        return nullptr;
    }

    const pointField& pts = edgeMesh_.points();

    pointField newPts(pts.size());

    # ifdef USE_OMP
    # pragma omp parallel for schedule(dynamic, 50)
    # endif
    forAll(pts, pointI)
    {
        newPts[pointI] =
            coordinateModifierPtr_->backwardModifiedPoint(pts[pointI]);
    }

    const edgeMesh* newEdgeMeshPtr = new edgeMesh(newPts, edgeMesh_.edges());

    return newEdgeMeshPtr;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
