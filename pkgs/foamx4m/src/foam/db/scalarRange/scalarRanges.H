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
	Foam::scalarRanges

Description
	A List of scalarRange.

SourceFiles
	scalarRanges.C

\*---------------------------------------------------------------------------*/

#ifndef scalarRanges_H
#define scalarRanges_H

#include "scalarRange.H"
#include "scalarList.H"
#include "boolList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class scalarRanges
:
	public List<scalarRange>
{
public:

	// Constructors

		//- Construct Null
		scalarRanges();

		//- Construct from Istream.
		//  The list items are comma-delimited.
		scalarRanges(Istream&);

	// Member Functions

		//- Return true if the given value is within the ranges
		bool selected(const scalar&) const;

		//- Return the set of selected entries in the given list
		//  that are within the ranges
		boolList selected(const scalarList&) const;

		//- Select a list of values that are within the ranges
		scalarList select(const scalarList&) const;

		//- Select a list of values that are within the ranges
		void inplaceSelect(List<scalar>&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
