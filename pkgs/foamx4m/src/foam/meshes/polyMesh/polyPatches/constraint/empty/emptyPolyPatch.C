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

#include "emptyPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(emptyPolyPatch, 0);

	addToRunTimeSelectionTable(polyPatch, emptyPolyPatch, word);
	addToRunTimeSelectionTable(polyPatch, emptyPolyPatch, dictionary);
}

// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::emptyPolyPatch::emptyPolyPatch
(
	const word& name,
	const label size,
	const label start,
	const label index,
	const polyBoundaryMesh& bm
)
:
	polyPatch(name, size, start, index, bm)
{}


Foam::emptyPolyPatch::emptyPolyPatch
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyBoundaryMesh& bm
)
:
	polyPatch(name, dict, index, bm)
{}


Foam::emptyPolyPatch::emptyPolyPatch
(
	const emptyPolyPatch& pp,
	const polyBoundaryMesh& bm,
	const label index,
	const label newSize,
	const label newStart
)
:
	polyPatch(pp, bm, index, newSize, newStart)
{}


Foam::emptyPolyPatch::emptyPolyPatch
(
	const emptyPolyPatch& pp
)
:
	polyPatch(pp)
{}


Foam::emptyPolyPatch::emptyPolyPatch
(
	const emptyPolyPatch& pp,
	const polyBoundaryMesh& bm
)
:
	polyPatch(pp, bm)
{}


// ************************************************************************* //
