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
	Foam::surfacePatchIOList

Description
	IOobject for a surfacePatchList

SourceFiles
	surfacePatchIOList.C

\*---------------------------------------------------------------------------*/

#ifndef surfacePatchIOList_H
#define surfacePatchIOList_H

#include "surfacePatchList.H"
#include "regIOobject.H"
#include "faceList.H"
#include "className.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
// Forward declaration of friend functions and operators
class surfacePatchIOList;

Ostream& operator<<(Ostream&, const surfacePatchIOList&);


class surfacePatchIOList
:
	public surfacePatchList,
	public regIOobject
{
	// Private data


	// Private Member Functions

		//- Disallow default bitwise copy construct
		surfacePatchIOList(const surfacePatchIOList&);

		//- Disallow default bitwise assignment
		void operator=(const surfacePatchIOList&);


public:

	//- Runtime type information
	TypeName("surfacePatchIOList");


	// Static data members

		//- Static data someStaticData

	// Constructors

		//- Construct from IOobject
		explicit surfacePatchIOList(const IOobject& io);

		//- Construct from IOobject
		surfacePatchIOList(const IOobject& io, const surfacePatchList&);

	// Destructor

		~surfacePatchIOList();


	// Member Functions

		//- writeData member function required by regIOobject
		bool writeData(Ostream&) const;


	// Member Operators

	// Friend Functions

	// Friend Operators

	// IOstream Operators

		friend Ostream& operator<<(Ostream&, const surfacePatchIOList&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
