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
	Foam::sigQuit

Description
	Signal handler for QUIT interupt.

	The standard interupt handler is overridden to ensure that the
	runningJob file is removed.

See Also
	Foam::JobInfo

SourceFiles
	sigQuit.C

\*---------------------------------------------------------------------------*/

#ifndef sigQuit_H
#define sigQuit_H

#include <signal.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class sigQuit
{
	// Private data

		//- Saved old signal trapping setting
		static struct sigaction oldAction_;


	// Private Member Functions

		//- Handler for caught signals
		static void sigHandler(int);


public:

	// Constructors

		//- Construct null
		sigQuit();


	//- Destructor
	~sigQuit();


	// Member functions

		//- Activate SIGQUIT signal handler
		void set(const bool verbose);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
