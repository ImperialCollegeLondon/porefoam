/*-------------------------------*- C++ -*-----------------------------------*\
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

Description
	Define the globals used in the OpenFOAM library.
	It is important that these are constructed in the appropriate order to
	avoid the use of unconstructed data in the global namespace.

	This file has the extension .Cver to trigger a Makefile rule that converts
	'VERSION\_STRING' and 'BUILD\_STRING' into the appropriate strings.

\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "foamVersion.H"

const char* const Foam::FOAMversion = WM_PROJECT_VERSION;  // compiler-flag: dep files won't recognise change of WM_PROJECT_VERSION, need cleaning to be updated, 
const char* const Foam::FOAMbuild = __DATE__;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Setup an error handler for the global new operator

#include "new.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create the nullObject singleton

#include "nullObject.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global IO streams

#include "IOstreams.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "JobInfo.H"
bool Foam::JobInfo::constructed(false);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global error definitions (initialised by construction)

#include "messageStream.C"
#include "error.C"
#include "IOerror.C"
#include "token.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Read the debug and info switches

#include "debug.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Read and set cell models

#include "globalCellModeller.C"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Create the jobInfo file in the $FOAM_JOB_DIR/runningJobs directory

#include "JobInfo.C"

// ************************************************************************* //
