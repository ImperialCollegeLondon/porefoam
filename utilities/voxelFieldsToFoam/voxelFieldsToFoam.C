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
//!   voxelFieldToFoam: map images, used to initialize Openfoam fields


#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>
#include "voxelImage.h"

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"

#include "OFstream.H"

#ifdef FOAMX
# include "foamTime.H"
#else
# include "Time.H"
#endif

#include <sys/stat.h> //mkdir

using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

template<typename T1,typename T2> void setEq(T1& v1,const T2& v2) { v1 = v2; }
template<typename T1> void setEq(T1& v1,const float3& v2) { v1[0]=v2[0];   v1[1]=v2[1]; v1[2]=v2[2]; }

template<typename ImgT, typename Type>
void convertToFoamWrite
(
	const word& alpha1Header ,
	const word& name ,
	Time& runTime,
	const fvMesh& mesh
)
{
	voxelImageT<ImgT> vximage(alpha1Header); //dummy

	int3 n=vximage.size3();
	dbl3 xmin=vximage.X0();
	dbl3 dx=vximage.dx();   //dx*=1e-6;

	if (!n[0]) Info<<"\nError: no image read\n"<<endl;



	const fvBoundaryMesh& boundary = mesh.boundary();



	GeometricField<Type, fvPatchField, volMesh> Vfield(
		IOobject ( name, "0", mesh, IOobject::MUST_READ, IOobject::AUTO_WRITE ),  mesh);


	const vectorField & C =	mesh.C().internalField();
	Type sumAlpha = Type();
	double sumWalpha=1e-64;
	if (vximage.nx())
	{
	  forAll(C,c)
	  {
		int i=(C[c][0]-xmin[0])/dx[0]*0.999999999999;
		int j=(C[c][1]-xmin[1])/dx[1]*0.999999999999;
		int k=(C[c][2]-xmin[2])/dx[2]*0.999999999999;
		setEq(Vfield[c],vximage(i,j,k));
		sumAlpha += Vfield[c];
		sumWalpha += 1.;
	  }
	}
	Info<<" AvgAlpha: "<<sumAlpha/sumWalpha<<endl;

	forAll(boundary, bi)
		Vfield.boundaryFieldRef()[bi]==Vfield.boundaryField()[bi].patchInternalField();


	OFstream AOF(Vfield.time().timeName()+"/"+name);
	Vfield.writeHeader(AOF);
	Vfield.writeData(AOF);
}



int main(int argc, char *argv[])
{

	argList::validArgs.append("fileList");
	//argList::validOptions.insert("fileList", "label");



	#include "setRootCase.H"
	#include "createTime.H"


	wordList fileList(args.argRead<wordList>(1));
	Info<<endl<<"fileList:    "<<fileList<<endl;

	#include "createMesh.H"

	std::string dummy; double tim=0.;
	scalar timO=tim;
	forAll(fileList,ii)
	{
		std::string basNam=fileList[ii];
		Info<<endl<<"basNam:    "<<basNam<<endl;
		std::stringstream ss;  ss.str(replaceFromTo(basNam,"_", " "));  ss>>dummy>>tim;  tim *= 4.267e-7; //time conversion s/fileSecondNumber
		scalar dt = tim-timO;   timO = tim;

		runTime.setDeltaT(dt);
		runTime++;
		Info<<endl<<"time: "<<runTime.timeName()<<"   timO: "<<timO<<"   dt: "<<dt<<endl;

		::mkdir(runTime.timeName().c_str(), 0733);

		convertToFoamWrite<float,Foam::scalar>(basNam+"_alpha.mhd", "alpha1", runTime,mesh);
		convertToFoamWrite<float,Foam::scalar>(basNam+"_pc.mhd", "pc", runTime,mesh);
		#ifdef VMMLIB__VECTOR__HPP
		convertToFoamWrite<float3,Foam::vector>(basNam+"_U.mhd", "U", runTime,mesh);
		#else
		Info<<basNam+"_U.mhd (float3) conversion to openfoam not supported"<<endl;
		#endif
		//basNam=""; ss>>basNam;
	}


	Info<< "end" << endl;

	return 0;
}


// ************************************************************************* //
