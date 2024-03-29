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

#include "wallDist.H"
#include "patchWave.H"
#include "fvMesh.H"
#include "wallPolyPatch.H"
#include "fvPatchField.H"
#include "Field.H"
#include "emptyFvPatchFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::wallDist::wallDist(const fvMesh& mesh, const bool correctWalls)
:
	volScalarField
	(
		IOobject
		(
			"y",
			mesh.time().timeName(),
			mesh
		),
		mesh,
		dimensionedScalar("y", dimLength, GREAT)
	),
	cellDistFuncs(mesh),
	correctWalls_(correctWalls),
	nUnset_(0)
{
	wallDist::correct();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::wallDist::~wallDist()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Correct for mesh geom/topo changes. Might be more intelligent in the
// future (if only small topology change)
void Foam::wallDist::correct()
{
	// AJ: make sure to pick up all patches that are specified as a wall
	const polyBoundaryMesh& bMesh = cellDistFuncs::mesh().boundaryMesh();
	labelHashSet wallPatchIDs(bMesh.size());

	forAll (bMesh, patchI)
	{
		if (bMesh[patchI].isWall())
		{
			wallPatchIDs.insert(patchI);
		}
	}

	// Calculate distance starting from wallPatch faces.
	patchWave wave(cellDistFuncs::mesh(), wallPatchIDs, correctWalls_);

	// Transfer cell values from wave into *this
	transfer(wave.distance());

	// Make near-wall distance consistent with wall distance
	// This is needed by immersed boundary walls
	// HJ, 29/May/2018
	const fvPatchList& patches = volScalarField::mesh().boundary();

	forAll (patches, patchI)
	{
		boundaryField()[patchI] = 1/patches[patchI].deltaCoeffs();
	}

	// Transfer number of unset values
	nUnset_ = wave.nUnset();
}


// ************************************************************************* //
