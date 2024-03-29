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
	Foam::tolerances

Description
	Selector class for solution tolerances.

SourceFiles
	tolerances.C

\*---------------------------------------------------------------------------*/

#ifndef tolerances_H
#define tolerances_H

#include "IOdictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class tolerances
:
	public IOdictionary
{
	// Private data

		dictionary relaxationFactors_;
		dictionary solverTolerances_;
		dictionary solverRelativeTolerances_;


	// Private Member Functions

		//- Disallow default bitwise copy construct and assignment
		tolerances(const tolerances&);
		void operator=(const tolerances&);


public:

	// Constructors

		//- Construct from time
		tolerances(const Time& t, const fileName& dictName);


	// Member Functions

		// Access

			bool relax(const word& name) const;
			scalar relaxationFactor(const word& name) const;

			scalar solverTolerance(const word& name) const;

			bool solverRelativeTolerances() const;
			scalar solverRelativeTolerance(const word& name) const;


		// Read

			//- Read the tolerances
			bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
