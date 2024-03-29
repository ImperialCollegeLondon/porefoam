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
	Foam::IOMap

Description
	A Map of objects of type \<T\> with automated input and output.

SourceFiles
	IOMap.C

\*---------------------------------------------------------------------------*/

#ifndef IOMap_H
#define IOMap_H

#include "Map.H"
#include "regIOobject.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class T>
class IOMap
:
	public regIOobject,
	public Map<T>
{

public:

	//- Runtime type information
	TypeName("Map");


	// Constructors

		//- Construct from IOobject
		IOMap(const IOobject&);

		//- Construct from IOobject and size of Map
		IOMap(const IOobject&, const label);

		//- Construct from IOobject and a Map
		IOMap(const IOobject&, const Map<T>&);

		//- Construct by transferring the Map contents
		IOMap(const IOobject&, const Xfer<Map<T> >&);


	//- Destructor
	virtual ~IOMap();


	// Member functions

		bool writeData(Ostream&) const;


	// Member operators

		//- Assignment of other IOMap's entries to this IOMap
		void operator=(const IOMap<T>&);

		//- Assignment of other Map's entries to this IOMap
		void operator=(const Map<T>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "IOMap.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
