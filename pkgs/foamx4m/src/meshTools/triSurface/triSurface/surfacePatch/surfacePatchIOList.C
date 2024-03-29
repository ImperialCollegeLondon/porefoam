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

#include "surfacePatchIOList.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(Foam::surfacePatchIOList, 0);


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

// Construct from IOObject
Foam::surfacePatchIOList::surfacePatchIOList
(
	const IOobject& io
)
:
	surfacePatchList(),
	regIOobject(io)
{
	Foam::string functionName =
		"surfacePatchIOList::surfacePatchIOList"
		"(const IOobject& io)";


	if (readOpt() == IOobject::MUST_READ)
	{
		surfacePatchList& patches = *this;

		// read polyPatchList
		Istream& is = readStream(typeName);

		PtrList<entry> patchEntries(is);
		patches.setSize(patchEntries.size());

		label faceI = 0;

		forAll(patches, patchI)
		{
			const dictionary& dict = patchEntries[patchI].dict();

			label patchSize = readLabel(dict.lookup("nFaces"));
			label startFaceI = readLabel(dict.lookup("startFace"));

			patches[patchI] =
				surfacePatch
				(
					word(dict.lookup("geometricType")),
					patchEntries[patchI].keyword(),
					patchSize,
					startFaceI,
					patchI
				);


			if (startFaceI != faceI)
			{
				FatalErrorIn(functionName)
					<< "Patches are not ordered. Start of patch " << patchI
					<< " does not correspond to sum of preceding patches."
					<< endl
					<< "while reading " << io.objectPath()
					<< exit(FatalError);
			}

			faceI += patchSize;
		}

		// Check state of IOstream
		is.check(functionName.c_str());

		close();
	}
}

// Construct from IOObject
Foam::surfacePatchIOList::surfacePatchIOList
(
	const IOobject& io,
	const surfacePatchList& patches
)
:
	surfacePatchList(patches),
	regIOobject(io)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::surfacePatchIOList::~surfacePatchIOList()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //


// writeData member function required by regIOobject
bool Foam::surfacePatchIOList::writeData(Ostream& os) const
{
	os << *this;
	return os.good();
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

Foam::Ostream& Foam::operator<<(Ostream& os, const surfacePatchIOList& patches)
{
	os  << patches.size() << nl << token::BEGIN_LIST;

	forAll(patches, patchI)
	{
		patches[patchI].writeDict(os);
	}

	os  << token::END_LIST;

	return os;
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Friend Functions  * * * * * * * * * * * * * //

// ************************************************************************* //
