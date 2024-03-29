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

\*---------------------------------------------------------------------------*/

#include "IOstreams.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Define the default IOstream versions and precision

const IOstream::versionNumber IOstream::originalVersion(0.5);
const IOstream::versionNumber IOstream::currentVersion(2.0);

Foam::debug::infoSwitch
IOstream::precision_
(
	"writePrecision",
	6
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Global IO streams

ISstream Sin(cin, "Sin");
OSstream Sout(cout, "Sout");
OSstream Serr(cerr, "Serr");
OFstream Snull("/dev/null");

prefixOSstream Pout(cout, "Pout");
prefixOSstream Perr(cerr, "Perr");

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
