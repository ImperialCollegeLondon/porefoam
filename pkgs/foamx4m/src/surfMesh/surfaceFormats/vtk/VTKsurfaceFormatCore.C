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
#include "VTKsurfaceFormatCore.H"
#include "clock.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::fileFormats::VTKsurfaceFormatCore::writeHeader
(
	Ostream& os,
	const pointField& pointLst
)
{
	// Write header
	os  << "# vtk DataFile Version 2.0" << nl
		<< "surface written " << clock::dateTime().c_str() << nl
		<< "ASCII" << nl
		<< nl
		<< "DATASET POLYDATA" << nl;

	// Write vertex coords
	os  << "POINTS " << pointLst.size() << " float" << nl;
	forAll(pointLst, ptI)
	{
		const point& pt = pointLst[ptI];

		os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
	}
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeTail
(
	Ostream& os,
	const UList<surfZone>& zoneLst
)
{
	label nFaces = 0;
	forAll(zoneLst, zoneI)
	{
		nFaces += zoneLst[zoneI].size();
	}

	// Print zone numbers
	os  << nl
		<< "CELL_DATA " << nFaces << nl
		<< "FIELD attributes 1" << nl
		<< "zone 1 " << nFaces << " int" << nl;


	forAll(zoneLst, zoneI)
	{
		forAll(zoneLst[zoneI], localFaceI)
		{
			if (localFaceI)
			{
				if (localFaceI % 20)
				{
					os << ' ';
				}
				else
				{
					os << nl;
				}
			}
			os  << zoneI + 1;
		}
		os  << nl;
	}
}


void Foam::fileFormats::VTKsurfaceFormatCore::writeTail
(
	Ostream& os,
	const UList<label>& zoneIds
)
{
	// Print zone numbers
	os  << nl
		<< "CELL_DATA " << zoneIds.size() << nl
		<< "FIELD attributes 1" << nl
		<< "zone 1 " << zoneIds.size() << " float" << nl;

	forAll(zoneIds, faceI)
	{
		if (faceI)
		{
			if (faceI % 20)
			{
				os << ' ';
			}
			else
			{
				os << nl;
			}
		}
		os  << zoneIds[faceI] + 1;
	}
	os  << nl;
}


// ************************************************************************* //
