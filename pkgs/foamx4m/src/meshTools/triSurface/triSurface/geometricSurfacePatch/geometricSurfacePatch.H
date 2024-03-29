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
	Foam::geometricSurfacePatch

Description
	The geometricSurfacePatch is like patchIdentifier but for surfaces.
	Holds type, name and index.

SourceFiles
	geometricSurfacePatch.C

\*---------------------------------------------------------------------------*/

#ifndef geometricSurfacePatch_H
#define geometricSurfacePatch_H

#include "word.H"
#include "label.H"
#include "typeInfo.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class dictionary;


class geometricSurfacePatch
{
	// Private data

		//- Type name of patch
		word geometricType_;

		//- Name of patch
		word name_;

		//- Index of patch in boundary
		label index_;


public:

	//- Runtime type information
	ClassName("geometricSurfacePatch");


	// Constructors

		//- Construct null
		geometricSurfacePatch();

		//- Construct from components
		geometricSurfacePatch
		(
			const word& geometricType,
			const word& name,
			const label index
		);

		//- Construct from Istream
		geometricSurfacePatch(Istream&, const label index);

		//- Construct from dictionary
		geometricSurfacePatch
		(
			const word& name,
			const dictionary& dict,
			const label index
		);


	// Member Functions

		//- Return name
		const word& name() const
		{
			return name_;
		}

		//- Return name
		word& name()
		{
			return name_;
		}

		//- Return the type of the patch
		const word& geometricType() const
		{
			return geometricType_;
		}

		//- Return the type of the patch
		word& geometricType()
		{
			return geometricType_;
		}

		//- Return the index of this patch in the boundaryMesh
		label index() const
		{
			return index_;
		}

		//- Return the index of this patch in the boundaryMesh
		label& index()
		{
			return index_;
		}

		//- Write
		void write(Ostream&) const;

		//- Write dictionary
		void writeDict(Ostream&) const;


	// Member Operators

		bool operator!=(const geometricSurfacePatch&) const;

		//- compare.
		bool operator==(const geometricSurfacePatch&) const;


	// Ostream Operator

		friend Ostream& operator<<(Ostream&, const geometricSurfacePatch&);
		friend Istream& operator>>(Istream&, geometricSurfacePatch&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
