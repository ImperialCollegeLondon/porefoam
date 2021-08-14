/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2012-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2018 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/


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

#include "VTKsurfaceFormat.H"
#include "labelIOField.H"
#include "scalarIOField.H"
#include "vectorIOField.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::writeHeaderPolygons
(
	Ostream& os,
	const UList<Face>& faceLst
)
{
	label nNodes = 0;

	forAll(faceLst, faceI)
	{
		nNodes += faceLst[faceI].size();
	}

	os  << nl
		<< "POLYGONS " << faceLst.size() << ' '
		<< faceLst.size() + nNodes << nl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

namespace Foam
{
template<class T>
void readBlock
(
    Istream& inFile,
    const label n,
    List<T>& lst
)
{
    lst.setSize(n);
    forAll(lst, i)
    {
        inFile >> lst[i];
    }
}

void readField
(
    ISstream& inFile,
    objRegistry& obj,
    const word& arrayName,
    const word& dataType,
    const label size
) 
{

	if(dataType=="VTK_INT"
		||dataType=="VTK_UINT"
		||dataType=="VTK_LONG"
		||dataType=="VTK_ULONG"
		||dataType=="int"
		||dataType=="VTK_ID")
	{
		labelField& fieldVals = dbset(obj, arrayName, new labelField(size));
		readBlock(inFile, fieldVals.size(), fieldVals);
	}


	else if(dataType=="VTK_FLOAT" || dataType=="VTK_DOUBLE")
	{
		scalarField& fieldVals = dbset(obj, arrayName, new scalarField(size));
		readBlock(inFile, fieldVals.size(), fieldVals);
	}
	else if(dataType=="VTK_STRING")
	{
		//if (debug)
		//{
			//Info<< "Reading strings:" << size << endl;
		//}
		autoPtr<stringList> fieldVals
		(
			new stringList
			(
				size
			)
		);
		// Consume current line.
		inFile.getLine(fieldVals()[0]);

		// Read without parsing
		forAll(fieldVals(), i)
		{
			inFile.getLine(fieldVals()[i]);
		}
		//regIOobject::store(fieldVals);
	}
	else
	{
		IOWarningInFunction(inFile)
			<< "Unhandled type " << dataType << endl
			<< "Skipping " << size
			<< " words." << endl;
		scalarField fieldVals;
		readBlock(inFile, size, fieldVals);
	}

}


wordList readFieldArray
(
    ISstream& inFile,
    objRegistry& obj,
    const label wantedSize
) 
{
    DynamicList<word> fields;

    word dataName(inFile);
    //if (debug)
    //{
        //Info<< "dataName:" << dataName << endl;
    //}
    label numArrays(readLabel(inFile));
    //if (debug)
    //{
        //Pout<< "numArrays:" << numArrays << endl;
    //}
    for (label i = 0; i < numArrays; i++)
    {
        word arrayName(inFile);
        label numComp(readLabel(inFile));
        label numTuples(readLabel(inFile));
        word dataType(inFile);

        //if (debug)
        //{
            //Info<< "Reading field " << arrayName
                //<< " of " << numTuples << " tuples of rank " << numComp << endl;
        //}

        if (wantedSize != -1 && numTuples != wantedSize)
        {
            FatalIOErrorInFunction(inFile)
                << "Expected " << wantedSize << " tuples but only have "
                << numTuples << exit(FatalIOError);
        }

        readField
        (
            inFile,
            obj,
            arrayName,
            dataType,
            numTuples*numComp
        );
        fields.append(arrayName);
    }
    return fields.shrink();
}

}

template<class Face>
Foam::fileFormats::VTKsurfaceFormat<Face>::VTKsurfaceFormat
(
    const fileName& filename
)
{
    read(filename);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Face>
bool Foam::fileFormats::VTKsurfaceFormat<Face>::read
(
	const fileName& filename
)
{
	const bool mustTriangulate = this->isTri();
	this->clear();

	IFstream is(filename);
	if (!is.good())
	{
        FatalErrorInFunction
			<< "Cannot read file " << filename
			<< exit(FatalError);
	}

    // assume that the groups are not intermixed
    bool sorted = true;


	//pointField points;
	DynamicList<point>  dynPoints;
	DynamicList<Face>  faces;

	//dynamicLabelList zones;
	dynamicLabelList zoneSizes;


	bool debug=true;

    //- Header
    string header_;

    //- Title
    string title_;

    //- DataType
    string dataType_;


    is.getLine(header_);
										if (debug)  Info<< "Header   : " << header_ << endl; 
    is.getLine(title_);
									if (debug)  Info<< "Title    : " << title_ << endl; 
    is.getLine(dataType_);
									if (debug)  Info<< "dataType : " << dataType_ << endl; 

										if (dataType_ == "BINARY")
										{
											FatalIOErrorInFunction(is)
												<< "Binary reading not supported " << exit(FatalIOError);
										}

    string readMode = "None";
    label wantedSize = -1;


    // Temporary storage for vertices of cells.
    labelList cellVerts;

    while (is.good())
    {
        word tag(is);

        if (!is.good())
        {
            break;
        }

        if (debug)
        {
            Info<< "line:" << is.lineNumber()
                << " tag:" << tag << endl;
        }

        if (tag == "DATASET")
        {
            word geomType(is);
            if (debug)
            {
                Info<< "geomType : " << geomType << endl;
            }
            //readMode = parseModeNames[geomType];
            wantedSize = -1;
        }
        else if (tag == "POINTS")
        {
            label nPoints(readLabel(is));
            dynPoints.setSize(nPoints);    ///3);
            if (debug)
            {
                Info<< "Reading " << nPoints << " numbers representing "
                    << dynPoints.size() << " coordinates." << endl;
            }

            word primitiveTag(is);
            if (primitiveTag != "float" && primitiveTag != "double")
            {
                FatalIOErrorInFunction(is)
                    << "Expected 'float' entry but found "
                    << primitiveTag
                    << exit(FatalIOError);
            }
            forAll(dynPoints, i)
            {
                is >> dynPoints[i].x() >> dynPoints[i].y() >> dynPoints[i].z();
            }
        }
        else if (tag == "CELLS")
        {
            label nCells(readLabel(is));
            label nNumbers(readLabel(is));
            if (debug)
            {
                Info<< "Reading " << nCells << " cells or faces." << endl;
            }
            readBlock(is, nNumbers, cellVerts);
        }
        else if (tag == "CELL_TYPES")
        {
            label nCellTypes(readLabel(is));

            labelList cellTypes;
            readBlock(is, nCellTypes, cellTypes);

            if (cellTypes.size() > 0 && cellVerts.size() == 0)
            {
                FatalIOErrorInFunction(is)
                    << "Found " << cellTypes.size()
                    << " cellTypes but no cells."
                    << exit(FatalIOError);
            }

            //extractCells(is, cellTypes, cellVerts);
            //cellVerts.clear();
        }
        else if (tag == "LINES")
        {
            label nLines(readLabel(is));
            label nNumbers(readLabel(is));
            //if (debug)
            {
                Info<< "Reading " << nLines << " lines." << endl;
            }
            labelList lineVerts;
            readBlock(is, nNumbers, lineVerts);
			labelListList lines_;

            label lineI = lines_.size();
            lines_.setSize(lineI+nLines);
            //lineMap_.setSize(lines_.size());

            label elemI = 0;
            for (label i = 0; i < nLines; i++)
            {
                //lineMap_[lineI] = lineI;
                labelList& f = lines_[lineI];
                f.setSize(lineVerts[elemI++]);
                forAll(f, fp)
                {
                    f[fp] = lineVerts[elemI++];
                }
                lineI++;
            }
        }
        else if (tag == "POLYGONS")
        {
            // If in polydata mode

            label nFaces(readLabel(is));
            label nNumbers(readLabel(is));
            if (debug)
            {
                Info<< "Reading " << nFaces << " faces." << endl;
            }
            labelList faceVerts;
            readBlock(is, nNumbers, faceVerts);

            label facei = faces.size();
            faces.setSize(facei+nFaces);
            //faceMap_.setSize(faces.size());

            label elemI = 0;
            for (label i = 0; i < nFaces; i++)
            {
                //faceMap_[facei] = facei;
                Face& f = faces[facei];
                f.setSize(faceVerts[elemI++]);
                forAll(f, fp)
                {
                    f[fp] = faceVerts[elemI++];
                }
                facei++;
            }
        }
        else if (tag == "POINT_DATA")
        {
            // 'POINT_DATA 24'
            readMode = "POINT_DATA";
            wantedSize = dynPoints.size();

            label nPoints(readLabel(is));
															if (nPoints != wantedSize)
															{
																FatalIOErrorInFunction(is)
																	<< "Reading POINT_DATA : expected " << wantedSize
																	<< " but read " << nPoints << exit(FatalIOError);
															}
        }
        else if (tag == "CELL_DATA")
        {
            readMode = "CELL_DATA";
            wantedSize = faces.size();//+cells_.size()+lines_.size();

            label nCells(readLabel(is));
															if (nCells != wantedSize)
															{
																FatalIOErrorInFunction(is)
																	<< "Reading CELL_DATA : expected "
																	<< wantedSize
																	<< " but read " << nCells << exit(FatalIOError);
															}
        }
        else if (tag == "FIELD")
        {
            // wantedSize already set according to type we expected to read.
            readFieldArray(is,selectRegistry(readMode), wantedSize);//, 
        }
        else if (tag == "SCALARS")
        {
            string line;
            is.getLine(line);
            IStringStream is(line);
            word dataName(is);
            word dataType(is);
            //label numComp(readLabel(is));

																	if (debug)
																	{
																		Info<< "Reading scalar " << dataName
																			<< " of type " << dataType
																			<< " from lookup table" << endl;
																	}

																word lookupTableTag(is);
																if (lookupTableTag != "LOOKUP_TABLE")
																{
																	FatalIOErrorInFunction(is)
																		<< "Expected tag LOOKUP_TABLE but read "
																		<< lookupTableTag
																		<< exit(FatalIOError);
																}

            word lookupTableName(is);

            readField
            (
                is,
                selectRegistry(readMode),//
                dataName,
                dataType,
                wantedSize//*numComp
            );
        }
        else if (tag == "VECTORS" || tag == "NORMALS")
        {
            // 'NORMALS Normals float'
            string line;
            is.getLine(line);
            IStringStream is(line);
            word dataName(is);
            word dataType(is);
																if (debug)
																{
																	Info<< "Reading vector " << dataName
																		<< " of type " << dataType << endl;
																}


            readField
            (
                is,
                selectRegistry(readMode),
                dataName,
                dataType,
                3*wantedSize
            );

            //if
            //( dataType == "VTK_FLOAT" || dataType == "VTK_DOUBLE")
            //{
                //objRegistry::iterator iter = reg.find(dataName);
                //scalarField s(*dynamic_cast<const scalarField*>(iter()));
                //reg.erase(iter);
                //autoPtr<vectorIOField> fieldVals
                //(
                    //new vectorIOField
                    //(
                        //IOobject
                        //(
                            //dataName,
                            //"",
                            //reg
                        //),
                        //s.size()/3
                    //)
                //);

                //label elemI = 0;
                //forAll(fieldVals(), i)
                //{
                    //fieldVals()[i].x() = s[elemI++];
                    //fieldVals()[i].y() = s[elemI++];
                    //fieldVals()[i].z() = s[elemI++];
                //}
                //regIOobject::store(fieldVals);
            //}
        }
        else if (tag == "TEXTURE_COORDINATES")
        {
            // 'TEXTURE_COORDINATES TCoords 2 float'
            string line;
            is.getLine(line);
            IStringStream is(line);
            word dataName(is);          //"Tcoords"
            label dim(readLabel(is));
            word dataType(is);

            if (debug)
            {
                Info<< "Reading texture coords " << dataName
                    << " dimension " << dim
                    << " of type " << dataType << endl;
            }

            scalarField coords(dim*dynPoints.size());
            readBlock(is, coords.size(), coords);
        }
        else if (tag == "TRIANGLE_STRIPS")
        {
            label nStrips(readLabel(is));
            label nNumbers(readLabel(is));
            if (debug)
            {
                Info<< "Reading " << nStrips << " triangle strips." << endl;
            }
            labelList faceVerts;
            readBlock(is, nNumbers, faceVerts);

            // Count number of triangles
            label elemI = 0;
            label nTris = 0;
            for (label i = 0; i < nStrips; i++)
            {
                label nVerts = faceVerts[elemI++];
                nTris += nVerts-2;
                elemI += nVerts;
            }


            // Store
            label facei = faces.size();
            faces.setSize(facei+nTris);
            //faceMap_.setSize(faces.size());
            elemI = 0;
            for (label i = 0; i < nStrips; i++)
            {
                label nVerts = faceVerts[elemI++];
                label nTris = nVerts-2;

                // Read first triangle
                //faceMap_[facei] = facei;
                Face& f = faces[facei++];
                f.setSize(3);
                f[0] = faceVerts[elemI++];
                f[1] = faceVerts[elemI++];
                f[2] = faceVerts[elemI++];
                for (label triI = 1; triI < nTris; triI++)
                {
                    //faceMap_[facei] = facei;
                    Face& f = faces[facei++];
                    f.setSize(3);
                    f[0] = faceVerts[elemI-1];
                    f[1] = faceVerts[elemI-2];
                    f[2] = faceVerts[elemI++];
                }
            }
        }
        else if (tag == "METADATA")
        {
            word infoTag(is);
            if (infoTag != "INFORMATION")
            {
                FatalIOErrorInFunction(is)                << "Unsupported tag "                << infoTag << exit(FatalIOError);
            }
            label nInfo(readLabel(is));
            if (debug)
            {
                Info<< "Consuming " << nInfo << " metadata information."              << endl;
            }
            string line;
            // Consume rest of line
            is.getLine(line);
            for (label i = 0; i < 2*nInfo; i++)
            {
                is.getLine(line);
            }
        }
        else
        {
            FatalIOErrorInFunction(is)                << "Unsupported tag "                << tag << exit(FatalIOError);
        }
    }

	//AQ copied from of-v1606+
    // Assume all faces in zone0 unless a region field is present
    labelList zones(faces.size(), 0); 
    Info<<"\ncellData_: ";
    for(auto&d:cellData_) Info<<d.first<<endl;
    if (cellData_.find("zone")!=cellData_.end())
    {
        const labelField& region = *static_cast<labelField*>(dbget(cellData_,"zone"));
        forAll(region, i)
        {
            zones[i] = label(region[i]);
        }
        Info<<"\n nZons: "<<max(zones)+1<<endl;
    }
    else if (cellData_.find("STLSolidLabeling")!=cellData_.end())
    {
        const labelField& region = *static_cast<labelField*>(dbget(cellData_,"STLSolidLabeling"));
        forAll(region, i)
        {
            zones[i] = label(region[i]);
        }
    }

    // Create zone names
    const label nZones = max(zones)+1;
    wordList zoneNames(nZones);
    forAll(zoneNames, i)
    {
        zoneNames[i] = "zone" + Foam::name(i);
    }





	// Read faces - ignore optional zone information     OFF format
	// use a DynamicList for possible on-the-fly triangulation

    label nTri = 0;
    if (mustTriangulate)
    {
        forAll(faces, faceI)
        {
            nTri += faces[faceI].size()-2;
        }
 
	 DynamicList<Face>  dynTriFs;
        dynamicLabelList dynZones(nTri);
        forAll(faces, faceI)
	 {


            const face& f = faces[faceI];

		if (f.size() > 3)
		{
			for (label fp1 = 1; fp1 < f.size() - 1; fp1++)
			{
				label fp2 = f.fcIndex(fp1);

				dynTriFs.append(triFace(f[0], f[fp1], f[fp2]));
                dynZones.append(zones[faceI]);
			}
		}
		else
		{
			dynTriFs.append(Face(f));
		}
	 }
        // Count
        labelList zoneSizes(nZones, 0);
        forAll(dynZones, triI)
        {
            zoneSizes[dynZones[triI]]++;
        }
	 
		this->storedPoints().transfer(dynPoints);		// transfer to normal lists


		//this->sortFacesAndStore(dynTriFs.xfer(), dynTriZones.xfer(), sorted);

		this->addZones(zoneSizes, zoneNames, true);			// add zones, culling empty ones
	}
	else
	{
		this->storedPoints().transfer(dynPoints);		// transfer to normal lists

        // Count
        labelList zoneSizes(nZones, 0);
        forAll(zones, faceI)
        {
            zoneSizes[zones[faceI]]++;
        }

		this->sortFacesAndStore(faces.xfer(), zones.xfer(), sorted);

		this->addZones(zoneSizes, zoneNames, true);		// add zones, culling empty ones
	}



	return true;
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
	const fileName& filename,
	const MeshedSurfaceProxy<Face>& surf
)
{
	
	OFstream os(filename);
	if (!os.good())
	{
		FatalErrorIn
		(
			"fileFormats::VTKsurfaceFormat::write"
			"(const fileName&, const MeshedSurfaceProxy<Face>&)"
		)
			<< "Cannot open file for writing " << filename
			<< exit(FatalError);
	}
	write(os,surf);
}
template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
	Ostream& os,
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
	  : VTKsurfaceFormat::oneZone(faceLst)
	);

	const bool useFaceMap = (surf.useFaceMap() && zones.size() > 1);


	writeHeader(os, pointLst);
	writeHeaderPolygons(os, faceLst);

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
				os << ' ' << nl;
			}
		}
		else
		{
			forAll(zone, localFaceI)
			{
				const Face& f = faceLst[faceIndex++];

				os << f.size();
				forAll(f, fp)
				{
					os << ' ' << f[fp];
				}
				os << ' ' << nl;
			}
		}
	}

	writeTail(os, zones);
}


template<class Face>
void Foam::fileFormats::VTKsurfaceFormat<Face>::write
(
	const fileName& filename,
	const UnsortedMeshedSurface<Face>& surf
)
{
	OFstream os(filename);
	if (!os.good())
	{
        FatalErrorInFunction
			<< "Cannot open file for writing " << filename
			<< exit(FatalError);
	}


	const List<Face>& faceLst = surf.faces();

	writeHeader(os, surf.points());
	writeHeaderPolygons(os, faceLst);

	forAll(faceLst, faceI)
	{
		const Face& f = faceLst[faceI];

		os << f.size();
		forAll(f, fp)
		{
			os << ' ' << f[fp];
		}
		os << ' ' << nl;
	}

	writeTail(os, surf.zoneIds());
}


// ************************************************************************* //
