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

#include "polyMeshGenModifier.H"
#include "demandDrivenData.H"
#include "labelList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void polyMeshGenModifier::removeUnusedVertices()
{
    faceListPMG& faces = mesh_.faces_;
    pointFieldPMG& points = mesh_.points_;

    boolList usePoint(points.size(), false);
    forAll(faces, faceI)
    {
        const face& f = faces[faceI];

        forAll(f, pI)
            usePoint[f[pI]] = true;
    }

    labelLongList newLabel(points.size(), -1);
    label nPoints(0);
    forAll(points, pI)
        if( usePoint[pI] )
            newLabel[pI] = nPoints++;

    //- remove unused points from the list
    forAll(newLabel, pI)
        if( (newLabel[pI] != -1) && (newLabel[pI] < pI) )
        {
            points[newLabel[pI]] = points[pI];
        }

    points.setSize(nPoints);

    forAll(faces, faceI)
    {
        face& f = faces[faceI];

        forAll(f, pI)
            f[pI] = newLabel[f[pI]];
    }

    mesh_.updatePointSubsets(newLabel);

    mesh_.clearOut();
    this->clearOut();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
