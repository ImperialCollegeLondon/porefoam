/*-------------------------------------------------------------------------*\
 Copyright (C) 2010-2020  Ali Qaseminejad Raeini

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/

//! Description:
//!   volume-preserving surface smoothing, retired due to incompatibility
//!   with some OpenFOAm versions

#include "MeshedSurfaces.H"

#include "argList.H"
#include "OFstream.H"
#include "boundBox.H"

using namespace Foam;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:
int appendUnique(DynamicList<label>& dynList, label value)
{
	bool append=true;
	forAll(dynList,i)
	{
		if(dynList[i]==value) append=false;
	}
	if(append)
	{
		dynList.append(value);
		return 1;
	}
	else
	{
		return 0;
	}
}


int main(int argc, char *argv[])
{

#   include "setRootCase.H"
#   include "createTime.H"

	IOdictionary meshingDict
	(	IOobject
		( "meshingDict", runTime.system(), runTime,
			IOobject::MUST_READ, IOobject::NO_WRITE
	)	);
	dictionary smoothingDict(meshingDict.subDict("surfaceSmoothing"));

	word surfFileName(smoothingDict.lookup("inputSurface"));



	scalar relax(readScalar(smoothingDict.lookup("relaxationFactor")));
	if ((relax < 0) || (relax > 1))
	{
		FatalErrorIn(args.executable()) << "Illegal relaxation factor "
		    << relax << endl
		    << "0: no change   1: move vertices to average of neighbours"
		    << exit(FatalError);
	}
	label kernelRadius(readLabel(smoothingDict.lookup("kernelRadius")));
	label iters(readLabel(smoothingDict.lookup("nIterations")));

	Info<< "Relax:" << relax << endl;
	Info<< "kernel radius:" << kernelRadius << endl;
	Info<< "Iters:" << iters << endl;


	Info<< "Reading surface from " << surfFileName << " ..." << endl;

	meshedSurface surf1(surfFileName);

	Info<< "nFaces    : " << surf1.size() << endl;
	Info<< "nVertices     : " << surf1.nPoints() << endl;
	Info<< "Bounding Box : " << boundBox(surf1.localPoints()) << endl;


	#define labelLoop face

	const pointField & points1=surf1.points();
	const List<face> & faces=surf1.faces();
	DynamicList<DynamicList<label> > pointPointsTmp(points1.size());
	List<labelLoop> pointPoints(points1.size());
	forAll(faces,faceI)
	{
		const face & f=faces[faceI];
		forAll(f, pI)
		{
		   pointPointsTmp[ f[pI] ].append(f.nextLabel(pI));
		}
	}
	forAll(pointPoints,i)
	{
		pointPoints[i].setSize(pointPoints[i].size());
		pointPoints[i]=face(pointPointsTmp[i]);

	}

	//const labelList& meshPoints = surf1.meshPoints();

	Info<<"........"<<endl<<endl ;



	pointField newPoints(points1);
	for(label iter = 0; iter < iters; iter++)
	{
		Info<<"\nSmoothing: "<<iter<<"  ";
		#include "./gauss.H"
	}
	Info<<":/"<<endl;


	surf1.movePoints(newPoints);


	fileName outFileName(smoothingDict.lookup("outputSurface"));

	Info<< "Writing surface to " << outFileName << " ..." << endl;

	//Info<<"triangulating the surface before write"<<endl;
	//surf1.triangulate();

	surf1.write(outFileName);

	Info << "End\n" << endl;

	return 0;
}


// ************************************************************************* //
