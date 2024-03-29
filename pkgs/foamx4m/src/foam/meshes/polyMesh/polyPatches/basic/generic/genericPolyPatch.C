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

#include "genericPolyPatch.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(genericPolyPatch, 0);

	addToRunTimeSelectionTable(polyPatch, genericPolyPatch, word);
	addToRunTimeSelectionTable(polyPatch, genericPolyPatch, dictionary);
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

Foam::genericPolyPatch::genericPolyPatch
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


Foam::genericPolyPatch::genericPolyPatch
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyBoundaryMesh& bm
)
:
	polyPatch(name, dict, index, bm),
	actualTypeName_(dict.lookup("type")),
	dict_(dict)
{}


Foam::genericPolyPatch::genericPolyPatch
(
	const genericPolyPatch& pp,
	const polyBoundaryMesh& bm,
	const label index,
	const label newSize,
	const label newStart
)
:
	polyPatch(pp, bm, index, newSize, newStart),
	actualTypeName_(pp.actualTypeName_),
	dict_(pp.dict_)
{}


Foam::genericPolyPatch::genericPolyPatch
(
	const genericPolyPatch& pp
)
:
	polyPatch(pp),
	actualTypeName_(pp.actualTypeName_),
	dict_(pp.dict_)
{}


Foam::genericPolyPatch::genericPolyPatch
(
	const genericPolyPatch& pp,
	const polyBoundaryMesh& bm
)
:
	polyPatch(pp, bm),
	actualTypeName_(pp.actualTypeName_),
	dict_(pp.dict_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::genericPolyPatch::~genericPolyPatch()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::genericPolyPatch::write(Ostream& os) const
{
	os.writeKeyword("type") << actualTypeName_ << token::END_STATEMENT << nl;
	patchIdentifier::write(os);
	os.writeKeyword("nFaces") << size() << token::END_STATEMENT << nl;
	os.writeKeyword("startFace") << start() << token::END_STATEMENT << nl;

	for
	(
		dictionary::const_iterator iter = dict_.begin();
		iter != dict_.end();
		++iter
	)
	{
		if
		(
			iter().keyword() != "type"
		 && iter().keyword() != "nFaces"
		 && iter().keyword() != "startFace"
		)
		{
			iter().write(os);
		}
	}
}


// ************************************************************************* //
