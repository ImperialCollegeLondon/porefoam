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

#include "polyPatch.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
	const word& patchType,
	const word& name,
	const label size,
	const label start,
	const label index,
	const polyBoundaryMesh& bm
)
{
	if (debug)
	{
		Info<< "polyPatch::New(const word&, const word&, const label, "
			   "const label, const label, const polyBoundaryMesh&) : "
			   "constructing polyPatch"
			<< endl;
	}

	wordConstructorTable::iterator cstrIter =
		wordConstructorTablePtr_->find(patchType);

	if (cstrIter == wordConstructorTablePtr_->end())
	{
		FatalErrorIn
		(
			"polyPatch::New(const word&, const word&, const label, "
			"const label, const label, const polyBoundaryMesh&) "
		)   << "Unknown polyPatch type " << patchType << " for patch " << name
			<< endl << endl
			<< "Valid polyPatch types are :" << endl
			<< wordConstructorTablePtr_->sortedToc()
			<< exit(FatalError);
	}

	return autoPtr<polyPatch>(cstrIter()(name, size, start, index, bm));
}


Foam::autoPtr<Foam::polyPatch> Foam::polyPatch::New
(
	const word& name,
	const dictionary& dict,
	const label index,
	const polyBoundaryMesh& bm
)
{
	if (debug)
	{
		Info<< "polyPatch::New(const word&, const dictionary&, const label, "
			   "const polyBoundaryMesh&) : constructing polyPatch"
			<< endl;
	}

	word patchType(dict.lookup("type"));

	dict.readIfPresent("geometricType", patchType);

	dictionaryConstructorTable::iterator cstrIter =
		dictionaryConstructorTablePtr_->find(patchType);

	if (cstrIter == dictionaryConstructorTablePtr_->end())
	{
		if (!disallowGenericPolyPatch)
		{
			cstrIter = dictionaryConstructorTablePtr_->find("genericPatch");
		}

		if (cstrIter == dictionaryConstructorTablePtr_->end())
		{
			FatalIOErrorIn
			(
				"polyPatch::New(const word&, const dictionary&, "
				"const label, const polyBoundaryMesh&)",
				dict
			)   << "Unknown polyPatch type " << patchType
				<< " for patch " << name
				<< endl << endl
				<< "Valid polyPatch types are :" << endl
				<< dictionaryConstructorTablePtr_->sortedToc()
				<< exit(FatalIOError);
		}
	}

	return autoPtr<polyPatch>(cstrIter()(name, dict, index, bm));
}


// ************************************************************************* //
