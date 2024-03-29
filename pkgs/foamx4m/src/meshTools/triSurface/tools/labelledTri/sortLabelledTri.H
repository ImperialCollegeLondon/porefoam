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
	Foam::sortLabelledTri

Description
	Helper class which when constructed with a triSurface
	sorts the faces according to region number (or rather constructs a
	mapping).

SourceFiles
	sortLabelledTri.C

\*---------------------------------------------------------------------------*/

#ifndef sortLabelledTri_H
#define sortLabelledTri_H

#include "labelList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class sortLabelledTri;
class triSurface;


//- Hold surface and label
class surfAndLabel
{
	const triSurface* surfPtr_;

	label index_;

	// Private Classes

		//- scalar comparison function used for sorting
		class less
		{
		public:

			inline bool operator()
			(
				const surfAndLabel& one,
				const surfAndLabel& two
			) const;
		};


public:

	friend class sortLabelledTri;

	// Constructors

		//- Construct null
		surfAndLabel()
		:
			surfPtr_(nullptr),
			index_(-1)
		{}

		//- Construct from surface and index
		surfAndLabel(const triSurface& surf, const label index)
		:
			surfPtr_(&surf),
			index_(index)
		{}
};



class sortLabelledTri
:
	public List<surfAndLabel>
{

public:

	// Constructors

		//- Construct from surface, sorting the faces according to patch
		sortLabelledTri(const triSurface&);


	// Member Functions

		// Access

			//- Set the labelList to those of sorted point indices
			void indices(labelList&) const;

			//- Return the list of sorted point indices
			labelList indices() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
