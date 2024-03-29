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

#include "SMESHsurfaceFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "OFstream.H"
#include "Ostream.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::SMESHsurfaceFormat<Face>::SMESHsurfaceFormat()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::SMESHsurfaceFormat<Face>::write
(
	const fileName& filename,
	const MeshedSurfaceProxy<Face>& surf
)
{
	const pointField& pointLst = surf.points();
	const List<Face>&  faceLst = surf.faces();
	const labelList& faceMap = surf.faceMap();

	const List<surfZone>& zones =
	(
		surf.surfZones().size() > 1
	  ? surf.surfZones()
	  : SMESHsurfaceFormat::oneZone(faceLst)
	);

	const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);


	OFstream os(filename);
	if (!os.good())
	{
		FatalErrorIn
		(
			"fileFormats::SMESHsurfaceFormat::write"
			"(const fileName&, const MeshedSurfaceProxy<Face>&)"
		)
			<< "Cannot open file for writing " << filename
			<< exit(FatalError);
	}


	// Write header
	os  << "# tetgen .smesh file written " << clock::dateTime().c_str() << nl
		<< "# <points count=\"" << pointLst.size() << "\">" << nl
		<< pointLst.size() << " 3" << nl;	// 3: dimensions

	// Write vertex coords
	forAll(pointLst, ptI)
	{
		const point& pt = pointLst[ptI];

		os  << ptI << ' ' << pt.x() << ' ' << pt.y() << ' ' << pt.z() << nl;
	}
	os  << "# </points>" << nl
		<< nl
		<< "# <faces count=\"" << faceLst.size() << "\">" << endl;

	os  << faceLst.size() << " 1" << endl;   // one attribute: zone number


	label faceIndex = 0;
	forAll(zones, zoneI)
	{
		const surfZone& zone = zones[zoneI];

		if (useFaceMap)
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceMap[faceIndex++]];

				os << f.size();
				forAll(f, fp)
				{
					os << ' ' << f[fp];
				}
				os << ' ' << zoneI << endl;
			}
		}
		else
		{
			forAll(zones[zoneI], localFaceI)
			{
				const Face& f = faceLst[faceIndex++];

				os << f.size();
				forAll(f, fp)
				{
					os << ' ' << f[fp];
				}
				os << ' ' << zoneI << endl;
			}
		}
	}

	// write tail

	os  << "# </faces>" << nl
		<< nl
		<< "# no holes or regions:" << nl
		<< '0' << nl		// holes
		<< '0' << endl;	 // regions
}


// ************************************************************************* //
