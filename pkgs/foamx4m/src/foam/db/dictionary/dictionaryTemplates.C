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

#include "dictionary.H"
#include "primitiveEntry.H"

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class T>
T Foam::dictionary::lookupOrDefault
(
	const word& keyword,
	const T& deflt,
	bool recursive,
	bool patternMatch
) const
{
	const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

	if (entryPtr)
	{
		return pTraits<T>(entryPtr->stream());
	}
	else
	{
		if (writeOptionalEntries)
		{
			IOInfoIn("dictionary::lookupOrDefault", *this)
				<< "Optional entry '" << keyword << "' is not present,"
				<< " returning the default value '" << deflt << "'"
				<< endl;
		}

		return deflt;
	}
}


template<class T>
T Foam::dictionary::lookupOrAddDefault
(
	const word& keyword,
	const T& deflt,
	bool recursive,
	bool patternMatch
)
{
	const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

	if (entryPtr)
	{
		return pTraits<T>(entryPtr->stream());
	}
	else
	{
		if (writeOptionalEntries)
		{
			IOInfoIn("dictionary::lookupOrAddDefault", *this)
				<< "Optional entry '" << keyword << "' is not present,"
				<< " adding and returning the default value '" << deflt << "'"
				<< endl;
		}

		add(new primitiveEntry(keyword, deflt));
		return deflt;
	}
}


template<class T>
bool Foam::dictionary::readIfPresent
(
	const word& keyword,
	T& val,
	bool recursive,
	bool patternMatch
) const
{
	const entry* entryPtr = lookupEntryPtr(keyword, recursive, patternMatch);

	if (entryPtr)
	{
		entryPtr->stream() >> val;
		return true;
	}
	else
	{
		if (writeOptionalEntries)
		{
			IOInfoIn("dictionary::readIfPresent", *this)
				<< "Optional entry '" << keyword << "' is not present,"
				<< " the default value '" << val << "' will be used."
				<< endl;
		}

		return false;
	}
}


template<class T>
void Foam::dictionary::add(const keyType& k, const T& t, bool overwrite)
{
	add(new primitiveEntry(k, t), overwrite);
}


template<class T>
void Foam::dictionary::set(const keyType& k, const T& t)
{
	set(new primitiveEntry(k, t));
}


// ************************************************************************* //
