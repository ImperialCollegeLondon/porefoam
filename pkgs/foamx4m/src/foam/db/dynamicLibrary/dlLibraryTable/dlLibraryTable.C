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

#include "dlLibraryTable.H"
#include "OSspecific.H"
#include "int.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(dlLibraryTable, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::dlLibraryTable::dlLibraryTable()
{}


Foam::dlLibraryTable::dlLibraryTable
(
	const dictionary& dict,
	const word& libsEntry
)
{
	open(dict, libsEntry);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::dlLibraryTable::~dlLibraryTable()
{
	forAllReverse(libPtrs_, i)
	{
		if (libPtrs_[i])
		{
			if (debug)
			{
				InfoInFunction
					<< "Closing " << libNames_[i]
					<< " with handle " << uintptr_t(libPtrs_[i]) << endl;
			}
			dlClose(libPtrs_[i]);
		}
	}
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::dlLibraryTable::open
(
	const fileName& functionLibName,
	const bool verbose
)
{
	if (functionLibName.size())
	{
		void* functionLibPtr = dlOpen(functionLibName, verbose);

		if (debug)
		{
			InfoInFunction
				<< "Opened " << functionLibName
				<< " resulting in handle " << uintptr_t(functionLibPtr) << endl;
		}

		if (!functionLibPtr)
		{
			if (verbose)
			{
				WarningInFunction
					<< "could not load " << functionLibName
					<< endl;
			}

			return false;
		}
		else
		{
			libPtrs_.append(functionLibPtr);
			libNames_.append(functionLibName);
			return true;
		}
	}
	else
	{
		return false;
	}
}


bool Foam::dlLibraryTable::close
(
	const fileName& functionLibName,
	const bool verbose
)
{
	label index = -1;
	forAllReverse(libNames_, i)
	{
		if (libNames_[i] == functionLibName)
		{
			index = i;
			break;
		}
	}

	if (index != -1)
	{
		if (debug)
		{
			InfoInFunction
				<< "Closing " << functionLibName
				<< " with handle " << uintptr_t(libPtrs_[index]) << endl;
		}

		bool ok = dlClose(libPtrs_[index]);

		libPtrs_[index] = nullptr;
		libNames_[index] = fileName::null;

		if (!ok)
		{
			if (verbose)
			{
				WarningInFunction
					<< "could not close " << functionLibName
					<< endl;
			}

			return false;
		}

		return true;
	}
	return false;
}


void* Foam::dlLibraryTable::findLibrary(const fileName& functionLibName)
{
	label index = -1;
	forAllReverse(libNames_, i)
	{
		if (libNames_[i] == functionLibName)
		{
			index = i;
			break;
		}
	}

	if (index != -1)
	{
		return libPtrs_[index];
	}
	return nullptr;
}


bool Foam::dlLibraryTable::open
(
	const dictionary& dict,
	const word& libsEntry
)
{
	if (dict.found(libsEntry))
	{
		fileNameList libNames(dict.lookup(libsEntry));

		bool allOpened = !libNames.empty();

		forAll(libNames, i)
		{
			allOpened = dlLibraryTable::open(libNames[i]) && allOpened;
		}

		return allOpened;
	}
	else
	{
		return false;
	}
}


// ************************************************************************* //
