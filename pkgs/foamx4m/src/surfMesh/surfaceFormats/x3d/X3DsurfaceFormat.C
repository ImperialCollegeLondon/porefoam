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

#include "X3DsurfaceFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "Ostream.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::X3DsurfaceFormat<Face>::X3DsurfaceFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::X3DsurfaceFormat<Face>::write
(
	const fileName& filename,
	const MeshedSurfaceProxy<Face>& surf
)
{
	const pointField& pointLst = surf.points();
	const List<Face>&  faceLst = surf.faces();
	const labelList& faceMap = surf.faceMap();

	// for no zones, suppress the group name
	const List<surfZone>& zones =
	(
		surf.surfZones().size() > 1
	  ? surf.surfZones()
	  : X3DsurfaceFormat::oneZone(faceLst, "")
	);

	const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

	OFstream os(filename);
	if (!os.good())
	{
		FatalErrorIn
		(
			"fileFormats::X3DsurfaceFormat::write"
			"(const fileName&, const MeshedSurfaceProxy<Face>&)"
		)
			<< "Cannot open file for writing " << filename
			<< exit(FatalError);
	}

	writeHeader(os);

	os  << "\n"
		"<Group>\n"
		" <Shape>\n";

	writeAppearance(os);


	// NOTE: we could provide an optimized IndexedTriangleSet output for
	// triangulated surfaces too

	os  <<
		"  <IndexedFaceSet coordIndex='\n";

	label faceIndex = 0;
	forAll(zones, zoneI)
	{
		const surfZone& zone = zones[zoneI];

		if (useFaceMap)
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceMap[faceIndex++]];

				forAll(f, fp)
				{
					os << f[fp] << ' ';
				}
				os << "-1\n";
			}
		}
		else
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceIndex++];

				forAll(f, fp)
				{
					os << f[fp] << ' ';
				}
				os << "-1\n";
			}
		}
	}

	os <<
		"' >\n"
		"    <Coordinate point='\n";

	// Write vertex coords
	forAll(pointLst, ptI)
	{
		const point& pt = pointLst[ptI];

		os  << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
	}

	os  <<
		"' />\n"        			   // end Coordinate
		"   </IndexedFaceSet>\n"
		"  </Shape>\n"
		" </Group>\n"
		"</X3D>\n";

}


// ************************************************************************* //
