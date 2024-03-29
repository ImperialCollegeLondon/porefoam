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

#include "includeEtcEntry.H"
#include "dictionary.H"
#include "IFstream.H"
#include "addToMemberFunctionSelectionTable.H"
#include "stringOps.H"
#include "OSspecific.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

const Foam::word Foam::functionEntries::includeEtcEntry::typeName
(
	Foam::functionEntries::includeEtcEntry::typeName_()
);

// Don't lookup the debug switch here as the debug switch dictionary
// might include includeEtcEntry
Foam::debug::debugSwitch
Foam::functionEntries::includeEtcEntry::debug
(
	"includeEntry",
	0
);

bool Foam::functionEntries::includeEtcEntry::report(false);


namespace Foam
{
namespace functionEntries
{
	addToMemberFunctionSelectionTable
	(
		functionEntry,
		includeEtcEntry,
		execute,
		dictionaryIstream
	);

	addToMemberFunctionSelectionTable
	(
		functionEntry,
		includeEtcEntry,
		execute,
		primitiveEntryIstream
	);
}
}

// * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * * //

Foam::fileName Foam::functionEntries::includeEtcEntry::includeEtcFileName
(
	const fileName& f,
	const dictionary& dict
)
{
	fileName fName(f);

	// Substitute dictionary and environment variables.
	// Allow empty substitutions.
	stringOps::inplaceExpand(fName, dict, true, true);

	if (fName.empty() || fName.isAbsolute())
	{
		return fName;
	}
	else
	{
		// Search the etc directories for the file
		return findEtcFile(fName);
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionEntries::includeEtcEntry::execute
(
	dictionary& parentDict,
	Istream& is
)
{
	const fileName rawFName(is);
	const fileName fName
	(
		includeEtcFileName(rawFName, parentDict)
	);
	IFstream ifs(fName);

	if (ifs)
	{
		if (Foam::functionEntries::includeEtcEntry::report)
		{
			Info<< fName << endl;
		}
		parentDict.read(ifs);
		return true;
	}
	else
	{
		FatalIOErrorIn
		(
			"functionEntries::includeEtcEntry::includeEtcEntry"
			"(dictionary& parentDict, Istream&)",
			is
		)   << "Cannot open etc file "
			<< (ifs.name().size() ? ifs.name() : rawFName)
			<< " while reading dictionary " << parentDict.name()
			<< exit(FatalIOError);

		return false;
	}
}


bool Foam::functionEntries::includeEtcEntry::execute
(
	const dictionary& parentDict,
	primitiveEntry& entry,
	Istream& is
)
{
	const fileName rawFName(is);
	const fileName fName
	(
		includeEtcFileName(rawFName, parentDict)
	);
	IFstream ifs(fName);

	if (ifs)
	{
		if (Foam::functionEntries::includeEtcEntry::report)
		{
			Info<< fName << endl;
		}
		entry.read(parentDict, ifs);
		return true;
	}
	else
	{
		FatalIOErrorIn
		(
			"functionEntries::includeEtcEntry::includeEtcEntry"
			"(dictionary& parentDict, primitiveEntry&, Istream&)",
			is
		)   << "Cannot open etc file "
			<< (ifs.name().size() ? ifs.name() : rawFName)
			<< " while reading dictionary " << parentDict.name()
			<< exit(FatalIOError);

		return false;
	}
}


// ************************************************************************* //
