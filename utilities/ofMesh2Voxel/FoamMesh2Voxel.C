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
//!   This app converts an unstructured openfoam mesh into a 3D image.


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
#include "voxelImageI.h"

using namespace Foam;



int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"


    IOdictionary meshingDict
    (
        IOobject
        (
            "meshingDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    word headerName(meshingDict.lookupOrDefault("headerName",word("")));
    word outputFormat(meshingDict.lookupOrDefault("outputFormat",word("binary")));


	int3 nnn;
	dbl3 xmin, dx;
    voxelImage rock(headerName);
	surfaceVectorField Cf =  mesh.Cf();
	surfaceVectorField Sf =  mesh.Sf();

    if (rock.nx()==0)
    {
		if(headerName.size()) Info<<"\n\nWarning: ignoring header file from system/meshingDict:"<<headerName<<endl<<endl;
		headerName="";
		Info <<"computing mesh extents"<<endl;
		dx[0]=std::pow(gAverage(mesh.V()),1./3.); dx[1]=dx[0]; dx[2]=dx[1];
		rock.dxCh()=dx;

		boundBox box(mesh.points());
		rock.X0Ch()[0]=box.min()[0];  rock.X0Ch()[1]=box.min()[1]; rock.X0Ch()[2]=box.min()[2];
		//rock.X0Ch()[1]-=dx[1]; rock.X0Ch()[2]-=dx[2];
		nnn[0]=((box.max()[0]-box.min()[0]) + 0.6*dx[0])/dx[0];
		nnn[1]=((box.max()[1]-box.min()[1]) + 0.6*dx[1])/dx[1];
		nnn[2]=((box.max()[2]-box.min()[2]) + 0.6*dx[2])/dx[2];

		Info<<"X0: "<<rock.X0()[0]<<" "<<rock.X0()[1]<<" "<<rock.X0()[2]<<endl;
		Info<<"dx:   "<<rock.dx()[0]<<"  "<<rock.dx()[1]<<"  "<<rock.dx()[2]<<endl;
		Info<<"nnn:    "<<nnn[0]<<" "<<nnn[1]<<" "<<nnn[2]<<endl;


		const dbl3 X0 = rock.X0();
		const dbl3 dx = rock.dx();
		label patchI = mesh.boundaryMesh().findPatchID("Left");
		if (patchI >= 0)
		{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[0];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[0]>.99 && mag(avgPos-X0[0]) < 2.*dx[0])
				rock.X0Ch()[0] = avgPos;
		}

		patchI = mesh.boundaryMesh().findPatchID("Right");
		if (patchI >= 0)
		{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[0];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[0]>.99 && mag(avgPos-X0[0]-nnn[0]*dx[0]) < 2.*dx[0])
				rock.dxCh()[0] = (avgPos-X0[0])/nnn[0];
		}

		patchI = mesh.boundaryMesh().findPatchID("Bottom");
		if (patchI >= 1)
		{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[1];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[1]>.99 && mag(avgPos-X0[1]) < 2.*dx[1])
				rock.X0Ch()[1] = avgPos;
		}
		patchI = mesh.boundaryMesh().findPatchID("Top");
		if (patchI >= 1)
		{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[1];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[1]>.99 && mag(avgPos-X0[1]-nnn[1]*dx[1]) < 2.*dx[1])
				rock.dxCh()[1] = (avgPos-X0[1])/nnn[1];
		}

		patchI = mesh.boundaryMesh().findPatchID("Back");
		if (patchI >= 2)
		{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[2];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[2]>.99 && mag(avgPos-X0[2]) < 2.*dx[2])
				rock.X0Ch()[2] = avgPos;
		}
		patchI = mesh.boundaryMesh().findPatchID("Front");
		if (patchI >= 0)
		{
			double avgPos = gAverage(Cf.boundaryField()[patchI])[2];
			vector avgN = gSum(Cf.boundaryField()[patchI]); avgN/=mag(avgN);
			if (avgN[2]>.99 && mag(avgPos-X0[2]-nnn[2]*dx[2]) < 2.*dx[2])
				rock.dxCh()[2] = (avgPos-X0[2])/nnn[2];
		}


		Info<<"->X0:   "<<rock.X0()[0]<<"  "<<rock.X0()[1]<<"  "<<rock.X0()[2]<<endl;
		Info<<"->dx:   "<<rock.dx()[0]<<"  "<<rock.dx()[1]<<"  "<<rock.dx()[2]<<endl;
		double dxAvg = (rock.dxCh()[0]+rock.dxCh()[1]+rock.dxCh()[2])/3.;

		rock.dxCh()[0]=dxAvg;  rock.dxCh()[1]=dxAvg;  rock.dxCh()[2]=dxAvg;
		Info<<"->dx:   "<<rock.dx()[0]<<"  "<<rock.dx()[1]<<"  "<<rock.dx()[2]<<endl;

		//rock.X0Ch()*=1e6;
		//rock.dxCh()*=1e6;
	}
	else
    {
		Info <<"mesh size read from header"<<endl;
		nnn=rock.size3();
		rock.reset(0,0,0,0);

	}













	xmin=rock.X0(); dx=rock.dx();
	//dx*=1e-6;
	//xmin*=1e-6;

	Info<<"xmin: "<<xmin[0]<<" "<<xmin[1]<<" "<<xmin[2]<<endl;
	Info<<"dx:   "<<dx[0]<<" "<<dx[1]<<" "<<dx[2]<<endl;
	Info<<"nnn:    "<<nnn[0]<<" "<<nnn[1]<<" "<<nnn[2]<<endl;



	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);




	//voxelField<double> pVoxel(nnn[0],nnn[1],nnn[2],0.);
	rock.reset(nnn[0],nnn[1],nnn[2],1);


	const vectorField & C =	mesh.C();
	const scalarField & V =	mesh.V();

	Info<<"mesh size: "<<C.size()<<endl;

	forAll(C,c)
	{
		int i=(C[c][0]-xmin[0])/dx[0];
		int j=(C[c][1]-xmin[1])/dx[1];
		int k=(C[c][2]-xmin[2])/dx[2];

		if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
		else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
	}

	//if (headerName.empty())
	{
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]+dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]-dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]-dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]-dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]+dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]-dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]-dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]-dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]+dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]-dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]+dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]+dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]+dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]-dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]+dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]+dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]+dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]-dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}
		forAll(C,c)
		{
			double dxi = std::pow(V[c],1./3.);
			int i=(C[c][0]-xmin[0]+dxi/4.)/dx[0];
			int j=(C[c][1]-xmin[1]+dxi/4.)/dx[1];
			int k=(C[c][2]-xmin[2]+dxi/4.)/dx[2];
			if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
			//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
		}

	}
	rock.FaceMedian06(2,4);
	if (headerName.empty())
	forAll(Cf,c)
	{

		int i=(Cf[c][0]-xmin[0])/dx[0];
		int j=(Cf[c][1]-xmin[1])/dx[1];
		int k=(Cf[c][2]-xmin[2])/dx[2];

		if (i>=0 && j>=0 && k>=0 && i<nnn[0] && j<nnn[1] && k<nnn[2])	rock(i,j,k)=0;
		//else							Info<<"Error: ijk: "<<i<<" "<<j<<" "<<k<<endl;
	}
	rock.FaceMedian06(1,5);
	rock.FaceMedian06(1,5);
	rock.FaceMedian06(1,5);

/*		const fvBoundaryMesh& patches = mesh.boundary();
		forAll(patches, patchI)
		{
			Info<<patches[patchI].name()<< mesh.Cf().boundaryField()[patchI].size()<<endl;

			const Field<vector> & Cfp =  mesh.Cf().boundaryField()[patchI];
			const Field<vector> & Sfp = mesh.Sf().boundaryField()[patchI];


			forAll(Cfp, c)
			{

				int i=(Cfp[c][0]-xmin[0]+dx[0]*0.1)/dx[0];
				int j=(Cfp[c][1]-xmin[1]+dx[1]*0.1)/dx[1];
				int k=(Cfp[c][2]-xmin[2]+dx[2]*0.1)/dx[2];

			  //if( i<nnn[0] && j<nnn[1] && k<nnn[2] )
			  {
					vector nf=(Sfp[c]/mag(Sfp[c]));
				if      (mag(nf[0])>0.95  &&  i<nnn[0] && j<nnn[1] && k<nnn[2])
					rock(i,j,k)=0;
				else if (mag(nf[1])>0.95  &&  i<nnn[0] && j<nnn[1] && k<nnn[2])
					rock(i,j,k)=0;
				else if (mag(nf[2])>0.95  &&  i<nnn[0] && j<nnn[1] && k<nnn[2])
					rock(i,j,k)=0;
				else
					Info<<c<<"  face not in x nor y nor z direction, or not on boundaries, nf = "<<nf<<"   i"<<i<<" j"<<j<<" k"<<k<<" v"<<rock(i,j,k)<<endl;
			  }

			}
		}*/


	//rock.crop(1,nnn[0]-2, 1,nnn[1]-2, 1,nnn[2]-2,1,1);

	rock.replacexLayer(0,1);
	rock.replacexLayer(nnn[0]-1,nnn[0]-2);

	rock.replaceyLayer(0,1);
	rock.replaceyLayer(nnn[1]-1,nnn[1]-2);

	rock.replacezLayer(0,1);
	rock.replacezLayer(nnn[2]-1,nnn[2]-2);

	if (outputFormat=="binary")
	{
		 rock.write(runTime.path()+"/vxlImage"+imgExt());
	}
	else
	{
		 rock.write(runTime.path()+"/vxlImage.dat");
    }



    Info<< "\nend" << endl;

    return 0;
}
