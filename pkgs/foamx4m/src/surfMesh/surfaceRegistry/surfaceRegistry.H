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
	Foam::surfaceRegistry

Description
	Wraps the normal objectRegistry with a local instance for %surfaces.

SourceFiles
	surfaceRegistry.C

\*---------------------------------------------------------------------------*/

#ifndef surfaceRegistry_H
#define surfaceRegistry_H

#include "objectRegistry.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class surfaceRegistry
:
	public objectRegistry
{

	// Private Member Functions

		//- Disallow default bitwise copy construct
		surfaceRegistry(const surfaceRegistry&);

		//- Disallow default bitwise assignment
		void operator=(const surfaceRegistry&);


public:

	//- Runtime type information
	TypeName("surfaceRegistry");

		//- The prefix to local: %surfaces
		static const word prefix;

		//- The default surface name: %default
		static word defaultName;


	// Constructors

		//- Construct for the given objectRegistry and named surface
		surfaceRegistry(const objectRegistry&, const word& surfName = "");


	// Destructor

		virtual ~surfaceRegistry();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
