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

#include "STARCDsurfaceFormat.H"
#include "ListOps.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
inline void Foam::fileFormats::STARCDsurfaceFormat<Face>::writeShell
(
	Ostream& os,
	const Face& f,
	const label cellId,
	const label cellTableId
)
{
	os  << cellId                    // includes 1 offset
		<< ' ' << starcdShellShape_  // 3(shell) shape
		<< ' ' << f.size()
		<< ' ' << cellTableId
		<< ' ' << starcdShellType_;  // 4(shell)

	// primitives have <= 8 vertices, but prevent overrun anyhow
	// indent following lines for ease of reading
	label count = 0;
	forAll(f, fp)
	{
		if ((count % 8) == 0)
		{
			os  << nl << "  " << cellId;
		}
		os  << ' ' << f[fp] + 1;
		count++;
	}
	os  << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Face>
Foam::fileFormats::STARCDsurfaceFormat<Face>::STARCDsurfaceFormat
(
	const fileName& filename
)
{
	read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::STARCDsurfaceFormat<Face>::read
(
	const fileName& filename
)
{
	const bool mustTriangulate = this->isTri();
	this->clear();

	fileName baseName = filename.lessExt();

	// STAR-CD index of points
	labelList pointId;

	// read points from .vrt file
	readPoints
	(
		IFstream(baseName + ".vrt")(),
		this->storedPoints(),
		pointId
	);

	// Build inverse mapping (STAR-CD pointId -> index)
	Map<label> mapPointId(2*pointId.size());
	forAll(pointId, i)
	{
		mapPointId.insert(pointId[i], i);
	}
	pointId.clear();

	//
	// read .cel file
	// ~~~~~~~~~~~~~~
	IFstream is(baseName + ".cel");
	if (!is.good())
	{
		FatalErrorIn
		(
			"fileFormats::STARCDsurfaceFormat::read(const fileName&)"
		)
			<< "Cannot read file " << is.name()
			<< exit(FatalError);
	}

	readHeader(is, "PROSTAR_CELL");

	DynamicList<Face>  dynFaces;
	dynamicLabelList dynZones;
	DynamicList<word>  dynNames;
	dynamicLabelList dynSizes;
	Map<label> lookup;

	// assume the cellTableIds are not intermixed
	bool sorted = true;
	label zoneI = 0;

	label lineLabel, shapeId, nLabels, cellTableId, typeId;
	dynamicLabelList vertexLabels(64);

	while ((is >> lineLabel).good())
	{
		is >> shapeId >> nLabels >> cellTableId >> typeId;

		vertexLabels.clear();
		vertexLabels.reserve(nLabels);

		// read indices - max 8 per line
		for (label i = 0; i < nLabels; ++i)
		{
			label vrtId;
			if ((i % 8) == 0)
			{
			   is >> lineLabel;
			}
			is >> vrtId;

			// convert original vertex id to point label
			vertexLabels.append(mapPointId[vrtId]);
		}

		if (typeId == starcdShellType_)
		{
			// Convert groupID into zoneID
			Map<label>::const_iterator fnd = lookup.find(cellTableId);
			if (fnd != lookup.end())
			{
				if (zoneI != fnd())
				{
					// cellTableIds are intermixed
					sorted = false;
				}
				zoneI = fnd();
			}
			else
			{
				zoneI = dynSizes.size();
				lookup.insert(cellTableId, zoneI);
				dynNames.append(word("cellTable_") + ::Foam::name(zoneI));
				dynSizes.append(0);
			}

			SubList<label> vertices(vertexLabels, vertexLabels.size());
			if (mustTriangulate && nLabels > 3)
			{
				face f(vertices);

				faceList triFaces(f.nTriangles());
				label nTri = 0;
				f.triangles(this->points(), nTri, triFaces);

				forAll(triFaces, faceI)
				{
					// a triangular face, but not yet a triFace
					dynFaces.append
					(
					    triFace
					    (
					        static_cast<UList<label>&>(triFaces[faceI])
					    )
					);
					dynZones.append(zoneI);
					dynSizes[zoneI]++;
				}
			}
			else
			{
				dynFaces.append(Face(vertices));
				dynZones.append(zoneI);
				dynSizes[zoneI]++;
			}
		}
	}
	mapPointId.clear();

	this->sortFacesAndStore(dynFaces.xfer(), dynZones.xfer(), sorted);

	// add zones, culling empty ones
	this->addZones(dynSizes, dynNames, true);
	return true;
}


template<class Face>
void Foam::fileFormats::STARCDsurfaceFormat<Face>::write
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
	  : STARCDsurfaceFormat::oneZone(faceLst)
	);

	const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);


	fileName baseName = filename.lessExt();

	writePoints(OFstream(baseName + ".vrt")(), pointLst);
	OFstream os(baseName + ".cel");
	writeHeader(os, "CELL");

	label faceIndex = 0;
	forAll(zones, zoneI)
	{
		const surfZone& zone = zones[zoneI];

		if (useFaceMap)
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceMap[faceIndex++]];
				writeShell(os, f, faceIndex, zoneI + 1);
			}
		}
		else
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceIndex++];
				writeShell(os, f, faceIndex, zoneI + 1);
			}
		}
	}

	// write simple .inp file
	writeCase
	(
		OFstream(baseName + ".inp")(),
		pointLst,
		faceLst.size(),
		zones
	);
}


// ************************************************************************* //
