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

Class
	Foam::wordReListMatcher

Description
	A wrapper for matching a List of wordRe.

Note
	The constructor should remain non-explicit. This allows automatic
	conversion from UList\<wordRe\> to wordReListMatcher in search
	functions.

SourceFiles
	wordReListMatcherI.H

\*---------------------------------------------------------------------------*/

#ifndef wordReListMatcher_H
#define wordReListMatcher_H

#include "wordReList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class wordReListMatcher
{
	// Private data

		//- Reference to underlying list
		const UList<wordRe>& reList_;

public:

	// Constructors

		//- Construct from a List of wordRe
		inline wordReListMatcher(const UList<wordRe>&);


	// Member Functions

		// Access

			inline label size() const;
			inline bool  empty() const;

			//- Return underlying list of wordRe
			inline const UList<wordRe>& operator()() const;


	// Searching

			//- Return true if string matches any of the regular expressions
			//  Smart match as regular expression or as a string.
			//  Optionally specify a literal match only.
			inline bool match(const string&, bool literalMatch=false) const;

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "wordReListMatcherI.H"

#endif

// ************************************************************************* //
