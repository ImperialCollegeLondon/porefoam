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

#include "objectRegistry.H"
#include "GTSsurfaceFormat.H"

#include "addToRunTimeSelectionTable.H"
#include "addToMemberFunctionSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace fileFormats
{

// read UnsortedMeshedSurface
addNamedTemplatedToRunTimeSelectionTable
(
	UnsortedMeshedSurface,
	GTSsurfaceFormat,
	face,
	fileExtension,
	gts
);
addNamedTemplatedToRunTimeSelectionTable
(
	UnsortedMeshedSurface,
	GTSsurfaceFormat,
	triFace,
	fileExtension,
	gts
);

// write MeshedSurface
addNamedTemplatedToMemberFunctionSelectionTable
(
	MeshedSurface,
	GTSsurfaceFormat,
	face,
	write,
	fileExtension,
	gts
);
addNamedTemplatedToMemberFunctionSelectionTable
(
	MeshedSurface,
	GTSsurfaceFormat,
	triFace,
	write,
	fileExtension,
	gts
);

// write UnsortedMeshedSurface
addNamedTemplatedToMemberFunctionSelectionTable
(
	UnsortedMeshedSurface,
	GTSsurfaceFormat,
	face,
	write,
	fileExtension,
	gts
);
addNamedTemplatedToMemberFunctionSelectionTable
(
	UnsortedMeshedSurface,
	GTSsurfaceFormat,
	triFace,
	write,
	fileExtension,
	gts
);

}
}

// ************************************************************************* //
