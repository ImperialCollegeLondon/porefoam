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

#include "OBJsurfaceFormat.H"
#include "clock.H"
#include "IFstream.H"
#include "IStringStream.H"
#include "Ostream.H"
#include "OFstream.H"
#include "ListOps.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::OBJsurfaceFormat<Face>::OBJsurfaceFormat
(
	const fileName& filename
)
{
	this->read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::OBJsurfaceFormat<Face>::read
(
	const fileName& filename
)
{
	const bool mustTriangulate = this->isTri();
	this->clear();

	IFstream is(filename);
	if (!is.good())
	{
		FatalErrorIn
		(
			"fileFormats::OBJsurfaceFormat::read(const fileName&)"
		)
			<< "Cannot read file " << filename
			<< exit(FatalError);
	}

	// assume that the groups are not intermixed
	bool sorted = true;

	DynamicList<point> dynPoints;
	DynamicList<Face>  dynFaces;
	dynamicLabelList dynZones;
	DynamicList<word>  dynNames;
	dynamicLabelList dynSizes;
	HashTable<label>   lookup;

	// place faces without a group in zone0
	label zoneI = 0;
	lookup.insert("zone0", zoneI);
	dynNames.append("zone0");
	dynSizes.append(0);

	while (is.good())
	{
		string line = this->getLineNoComment(is);

		// handle continuations
		if (line[line.size()-1] == '\\')
		{
			line.substr(0, line.size()-1);
			line += this->getLineNoComment(is);
		}

		// Read first word
		IStringStream lineStream(line);
		word cmd;
		lineStream >> cmd;

		if (cmd == "v")
		{
			scalar x, y, z;
			lineStream >> x >> y >> z;
			dynPoints.append(point(x, y, z));
		}
		else if (cmd == "g")
		{
			word name;
			lineStream >> name;

			HashTable<label>::const_iterator fnd = lookup.find(name);
			if (fnd != lookup.end())
			{
				if (zoneI != fnd())
				{
					// group appeared out of order
					sorted = false;
				}
				zoneI = fnd();
			}
			else
			{
				zoneI = dynSizes.size();
				lookup.insert(name, zoneI);
				dynNames.append(name);
				dynSizes.append(0);
			}
		}
		else if (cmd == "f")
		{
			dynamicLabelList dynVertices;

			// Assume 'f' is followed by space.
			string::size_type endNum = 1;

			while (true)
			{
				string::size_type startNum =
					line.find_first_not_of(' ', endNum);

				if (startNum == string::npos)
				{
					break;
				}

				endNum = line.find(' ', startNum);

				string vertexSpec;
				if (endNum != string::npos)
				{
					vertexSpec = line.substr(startNum, endNum-startNum);
				}
				else
				{
					vertexSpec = line.substr(startNum, line.size() - startNum);
				}

				string::size_type slashPos = vertexSpec.find('/');

				label vertI = 0;
				if (slashPos != string::npos)
				{
					IStringStream intStream(vertexSpec.substr(0, slashPos));

					intStream >> vertI;
				}
				else
				{
					IStringStream intStream(vertexSpec);

					intStream >> vertI;
				}
				dynVertices.append(vertI - 1);
			}
			dynVertices.shrink();

			UList<label>& f = static_cast<UList<label>&>(dynVertices);

			if (mustTriangulate && f.size() > 3)
			{
				// simple face triangulation about f[0]
				// points may be incomplete
				for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
				{
					label fp2 = f.fcIndex(fp1);

					dynFaces.append(triFace(f[0], f[fp1], f[fp2]));
					dynZones.append(zoneI);
					dynSizes[zoneI]++;
				}
			}
			else
			{
				dynFaces.append(Face(f));
				dynZones.append(zoneI);
				dynSizes[zoneI]++;
			}
		}
	}


	// transfer to normal lists
	this->storedPoints().transfer(dynPoints);

	this->sortFacesAndStore(dynFaces.xfer(), dynZones.xfer(), sorted);

	// add zones, culling empty ones
	this->addZones(dynSizes, dynNames, true);
	return true;
}


template<class Face>
void Foam::fileFormats::OBJsurfaceFormat<Face>::write
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
	  : OBJsurfaceFormat::oneZone(faceLst, "")
	);

	const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);

	OFstream os(filename);
	if (!os.good())
	{
		FatalErrorIn
		(
			"fileFormats::OBJsurfaceFormat::write"
			"(const fileName&, const MeshedSurfaceProxy<Face>&)"
		)
			<< "Cannot open file for writing " << filename
			<< exit(FatalError);
	}


	os  << "# Wavefront OBJ file written " << clock::dateTime().c_str() << nl
		<< "o " << os.name().lessExt().name() << nl
		<< nl
		<< "# points : " << pointLst.size() << nl
		<< "# faces  : " << faceLst.size() << nl
		<< "# zones  : " << zones.size() << nl;

	// Print zone names as comment
	forAll(zones, zoneI)
	{
		os  << "#	" << zoneI << "  " << zones[zoneI].name()
			<< "  (nFaces: " << zones[zoneI].size() << ")" << nl;
	}

	os  << nl
		<< "# <points count=\"" << pointLst.size() << "\">" << nl;

	// Write vertex coords
	forAll(pointLst, ptI)
	{
		const point& pt = pointLst[ptI];

		os  << "v " << pt.x() << ' '  << pt.y() << ' '  << pt.z() << nl;
	}

	os  << "# </points>" << nl
		<< nl
		<< "# <faces count=\"" << faceLst.size() << "\">" << endl;


	label faceIndex = 0;
	forAll(zones, zoneI)
	{
		const surfZone& zone = zones[zoneI];

		if (zone.name().size())
		{
			os << "g " << zone.name() << endl;
		}

		if (useFaceMap)
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceMap[faceIndex++]];

				os << 'f';
				forAll(f, fp)
				{
					os << ' ' << f[fp] + 1;
				}
				os << endl;
			}
		}
		else
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceIndex++];

				os << 'f';
				forAll(f, fp)
				{
					os << ' ' << f[fp] + 1;
				}
				os << endl;
			}
		}
	}
	os << "# </faces>" << endl;
}


// ************************************************************************* //
