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

#include "objectRegistry.H"
#include "surfaceFormatsCore.H"

#include "foamTime.H"
#include "IFstream.H"
#include "OFstream.H"
#include "SortableList.H"
#include "surfMesh.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

Foam::word Foam::fileFormats::surfaceFormatsCore::nativeExt("ofs");

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

Foam::string Foam::fileFormats::surfaceFormatsCore::getLineNoComment
(
	IFstream& is
)
{
	string line;
	do
	{
		is.getLine(line);
	}
	while ((line.empty() || line[0] == '#') && is.good());

	return line;
}


#if 0
Foam::fileName Foam::fileFormats::surfaceFormatsCore::localMeshFileName
(
	const word& surfName
)
{
	const word name(surfName.size() ? surfName : surfaceRegistry::defaultName);

	return fileName
	(
		surfaceRegistry::prefix/name/surfMesh::meshSubDir
	  / name + "." + nativeExt
	);
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::findMeshInstance
(
	const Time& t,
	const word& surfName
)
{
	fileName localName = localMeshFileName(surfName);

	// Search back through the time directories list to find the time
	// closest to and lower than current time

	instantList ts = t.times();
	label instanceI;

	for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
	{
		if (ts[instanceI].value() <= t.timeOutputValue())
		{
			break;
		}
	}

	// Noting that the current directory has already been searched
	// for mesh data, start searching from the previously stored time directory

	if (instanceI >= 0)
	{
		for (label i = instanceI; i >= 0; --i)
		{
			if (isFile(t.path()/ts[i].name()/localName))
			{
				return ts[i].name();
			}
		}
	}

	return "constant";
}


Foam::fileName Foam::fileFormats::surfaceFormatsCore::findMeshFile
(
	const Time& t,
	const word& surfName
)
{
	fileName localName = localMeshFileName(surfName);

	// Search back through the time directories list to find the time
	// closest to and lower than current time

	instantList ts = t.times();
	label instanceI;

	for (instanceI = ts.size()-1; instanceI >= 0; --instanceI)
	{
		if (ts[instanceI].value() <= t.timeOutputValue())
		{
			break;
		}
	}

	// Noting that the current directory has already been searched
	// for mesh data, start searching from the previously stored time directory

	if (instanceI >= 0)
	{
		for (label i = instanceI; i >= 0; --i)
		{
			fileName testName(t.path()/ts[i].name()/localName);

			if (isFile(testName))
			{
				return testName;
			}
		}
	}

	// fallback to "constant"
	return t.path()/"constant"/localName;
}
#endif


bool Foam::fileFormats::surfaceFormatsCore::checkSupport
(
	const wordHashSet& available,
	const word& ext,
	const bool verbose,
	const word& functionName
)
{
	if (available.found(ext))
	{
		return true;
	}
	else if (verbose)
	{
		wordList toc = available.toc();
		SortableList<word> known(toc.xfer());

		Info<<"Unknown file extension for " << functionName
			<< " : " << ext << nl
			<<"Valid types: (";
		// compact output:
		forAll(known, i)
		{
			Info<<" " << known[i];
		}
		Info<<" )" << endl;
	}

	return false;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::fileFormats::surfaceFormatsCore::surfaceFormatsCore()
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::fileFormats::surfaceFormatsCore::~surfaceFormatsCore()
{}


// ************************************************************************* //
