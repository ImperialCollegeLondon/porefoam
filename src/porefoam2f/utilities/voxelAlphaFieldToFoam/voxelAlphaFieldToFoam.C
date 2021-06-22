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
//!   map images, used to initialize Openfoam fields


    #include <fstream>
    #include <iostream>
    #include <vector>

    #include <assert.h>

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"

#include "OFstream.H"
 
#include "voxelImage.h"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:





int main(int argc, char *argv[])
{
	
    argList::validArgs.append("alpha1Header");
    argList::validOptions.insert("invertAlpha","");
    argList::validOptions.insert("nGrowAlpha", "label");
    //argList::addOption
    //(
        //"alpha1",
        //" create alpha1 "
    //);


	
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"



    
    word alpha1Header(args.argRead<word>(1));
   Info<<endl<<"alpha1Header:    "<<alpha1Header<<endl;



	//word headerName(meshingDict.lookup("headerName"));
	//word inputName(meshingDict.lookup("inputName"));

	//unsigned int n[3];
	//double  xmin[3];
	voxelImage vximage(alpha1Header); //dummy

	int3 n=vximage.size3();
	dbl3 xmin=vximage.X0(); //xmin*=1e-6;
	dbl3 dx=vximage.dx();   //dx*=1e-6;

	if (!n[0]) Info<<"\nError: no image read\n"<<endl;
	if (args.optionFound("invertAlpha"))	vximage.threshold101(1,255);


	runTime.setTime(timeDirs[timeDirs.size()-1], 0);


	 


	


	const fvBoundaryMesh& boundary = mesh.boundary();


	
	volScalarField alpha1
	(
		IOobject
		(
			"alpha1",
			runTime.timeName(),
			mesh,
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		mesh
	);	
		

		const vectorField & C =	mesh.C().internalField();
		double sumAlpha=0.0, sumWalpha=1e-64;
		if (vximage.nx())
		{
		  forAll(C,c)
		  {
			int i=(C[c][0]-xmin[0])/dx[0]*0.999999999999;
			int j=(C[c][1]-xmin[1])/dx[1]*0.999999999999;
			int k=(C[c][2]-xmin[2])/dx[2]*0.999999999999;
			alpha1[c]=vximage(i,j,k);
			sumAlpha += alpha1[c];
			sumWalpha += 1.0;
		  }
		}
		Info<<" AvgAlpha: "<<sumAlpha/sumWalpha<<endl;

		forAll(boundary, bi)
			alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();


		if (args.optionFound("nGrowAlpha"))
		{
			label nGrowAlpha(0);
			args.optionLookup("nGrowAlpha")() >> nGrowAlpha ;
			Info<<" nGrowAlpha  "<<nGrowAlpha<<endl;
			alpha1.correctBoundaryConditions();
			alpha1=fvc::average(linearInterpolate(alpha1));
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-90.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-150.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-60.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-180.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-30.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-210.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-120.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-120.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();
			alpha1=min(max(fvc::average(linearInterpolate(alpha1))*(241.0)-120.0,0.0),1.0);
					forAll(boundary, bi)	alpha1.boundaryFieldRef()[bi]==alpha1.boundaryField()[bi].patchInternalField();

			if(nGrowAlpha<0)
			{
			
				for (label itr=0;itr<-nGrowAlpha;++itr)
				{
					alpha1.correctBoundaryConditions();
					alpha1=1.0-min(fvc::average(linearInterpolate(1.0-alpha1))*(241.0),1.0);				
				}
			}
			else if(nGrowAlpha>0)
			{				
				for (label itr=0;itr<nGrowAlpha;++itr)
				{
					alpha1.correctBoundaryConditions();
					alpha1=min(fvc::average(linearInterpolate(alpha1))*(241.0),1.0);				
				}
			}

			alpha1.correctBoundaryConditions();
		}
		else
			Info<<"Warning use -nGrowAlpha 0 to apply some smoothing"<<endl;

		OFstream AOF(alpha1.time().timeName()+"/alpha1");
		alpha1.writeHeader(AOF);
		alpha1.writeData(AOF);

    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
