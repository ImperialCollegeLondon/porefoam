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

#include "ReadFields.H"
#include "HashSet.H"
#include "Pstream.H"
#include "IOobjectList.H"

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

// Read all fields of type. Returns names of fields read. Guarantees all
// processors to read fields in same order.
template<class GeoField, class Mesh>
Foam::wordList Foam::ReadFields
(
	const Mesh& mesh,
	const IOobjectList& objects,
	PtrList<GeoField>& fields,
	const bool syncPar
)
{
	// Search list of objects for wanted type
	IOobjectList fieldObjects(objects.lookupClass(GeoField::typeName));

	wordList masterNames(fieldObjects.names());

	if (syncPar && Pstream::parRun())
	{
		// Check that I have the same fields as the master
		const wordList localNames(masterNames);
		Pstream::scatter(masterNames);

		HashSet<word> localNamesSet(localNames);

		forAll(masterNames, i)
		{
			const word& masterFld = masterNames[i];

			HashSet<word>::iterator iter = localNamesSet.find(masterFld);

			if (iter == localNamesSet.end())
			{
				FatalErrorIn
				(
					"ReadFields<class GeoField, class Mesh>"
					"(const Mesh&, const IOobjectList&, PtrList<GeoField>&"
					", const bool)"
				)   << "Fields not synchronised across processors." << endl
					<< "Master has fields " << masterNames
					<< "  processor " << Pstream::myProcNo()
					<< " has fields " << localNames << exit(FatalError);
			}
			else
			{
				localNamesSet.erase(iter);
			}
		}

		forAllConstIter(HashSet<word>, localNamesSet, iter)
		{
			FatalErrorIn
			(
				"ReadFields<class GeoField, class Mesh>"
				"(const Mesh&, const IOobjectList&, PtrList<GeoField>&"
				", const bool)"
			)   << "Fields not synchronised across processors." << endl
				<< "Master has fields " << masterNames
				<< "  processor " << Pstream::myProcNo()
				<< " has fields " << localNames << exit(FatalError);
		}
	}


	fields.setSize(masterNames.size());

	// Make sure to read in masterNames order.

	forAll(masterNames, i)
	{
		Info<< "Reading " << GeoField::typeName << ' ' << masterNames[i]
			<< endl;

		const IOobject& io = *fieldObjects[masterNames[i]];

		fields.set
		(
			i,
			new GeoField
			(
				IOobject
				(
					io.name(),
					io.instance(),
					io.local(),
					io.db(),
					IOobject::MUST_READ,
					IOobject::AUTO_WRITE,
					io.registerObject()
				),
				mesh
			)
		);
	}
	return masterNames;
}


// ************************************************************************* //
