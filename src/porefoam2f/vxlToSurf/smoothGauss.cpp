// Volume-preserving Gaussian smoothing

// Copyright (C) 2018  IC
// Copyright (C) 2020  Ali Qaseminejad Raeini 

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// For any queries, contact Ali Q. Raeini: email aliqasemi@gmail.com




#include "InputFile.h"
#include "surfUtils.h"

using namespace std;
int  smoothSurf(InputFile& meshingDict, facePieceList& facezsZ, piece<point> pointsAll)
{
	Info<<"\n... smoothing surface ..."<<endl<<endl ;

	//word surfFileName(meshingDict.getOr(std::string("Chomic251Cyl.vtk"),"inputSurface"));
	//fileName outFileName(meshingDict.getOr(std::string("surfSmooth.vtk"),"outputSurface"));
	int solidIndex(meshingDict.getOr(1024,"solidIndex"));
	const int OWInd(meshingDict.getOr(512,"oilWaterIndex"));

	scalar relaxParIntrf(meshingDict.getOr(0.1,"relaxParIntrf"));
	int kernelRadius(meshingDict.getOr(6,"kernelRadiusGaussVP"));
	int nIters(meshingDict.getOr(8,"nIterationsGaussVP"));
	scalar relax(meshingDict.getOr(0.2,"relaxFactorGaussVP"));  ensure(0<=relax && relax<=1, "Illegal relaxation factor:\n 0: no change, 1: move vertices to average of neighbours",2);
	scalar relaxCL(meshingDict.getOr(0.1,"relaxFactorGaussVPCL"));
	Info<< "Relax:" << relax << endl;
	Info<< "Relax CL:" << relaxCL << endl;
	Info<< "kernel radius:" << kernelRadius << endl;
	Info<< "Iters:" << nIters << endl;

	ensure(relaxCL<=0.00000001,"relaxCL does not produce good results",0); // 

	labelList pMarks(pointsAll.size(),0);
	for(const auto& facezs:facezsZ) for(const auto& fac:facezs) for(const auto& pI:fac) // this only guaranties up to  3faces contact lines work ok
	{  int zn=fac.zone+1, pm=pMarks[pI]&255;  pMarks[pI] |= ((pm^zn)*(pm!=0))<<8 | zn;  } //syncCalcPMarks

	for(int iter = 0; iter < nIters; iter++)
	{///    Volume-preserving Gaussian smoothing
	 if(iter==nIters-1) writeSurfaceFiles(facezsZ,pointsAll, "dumpSurfSmooth.vtk");
	  //for_(facezsZ, ifcs) 
	 for (int ifcs=facezsZ.size()-1; ifcs>=0; --ifcs) if(facezsZ[ifcs].size())
	 {const auto& facezs=facezsZ[ifcs]; 
	vars<point> newPoints(pointsAll);



		//const labelListList pFaces  = pointFaces(newPoints.size(),faces);
		//const labelListList pNeips = getPointPoints(newPoints.size(),facezs);

		//dbls pWeights(pNeips.size(),1.0);
		vectorField pNw(pointsAll.size(),dbl3(0,0,0));
		dbls pAreas;
	  {// weits


		Info<< "\npMarks"<<ifcs<<": " << min(pMarks)<<" - "<< max(pMarks)<< "  "<<endl ;

		///. point-normal vectors
		 //vars<vectorField>  Cf=faceCentres(facezsZ,newPoints);

		 //vectorField pNs(pNeips.size(),dbl3(0,0,0));
		
		 for(const auto& fac:facezs)   //for(const auto& facezs:facezsZ) 
		 {
			dbl3 cf = centre(fac,newPoints);
			for_(fac, pI)
			{   int e0=fac[pI], e1=fac[(pI+1)%4];
				dbl3 Ce=0.5*(newPoints[e0] + newPoints[e1]); 

				 dbl3 pNE=0.5*((newPoints[e0]-cf) ^ (Ce-cf)); 
				  pNw[e0] += pNE;
				 //if( pMarks[e0]==fac.zone ) {    pNw[e0] += pNE;     }// pNs[e0] += pNE;  
				 //else if( fac.zone != OWInd )        pNw[e0] += pNE;
				 //else                                pNs[e0] += pNE;

				 pNE=0.5*((newPoints[e1]-cf) ^ (Ce-cf));
				  pNw[e1] -= pNE; 
				 //if( pMarks[e1]==fac.zone ) {    pNw[e1] -= pNE;    } // pNs[e1] -= pNE;  
				 //else if( fac.zone != OWInd )        pNw[e1] -= pNE;
				 //else                                pNs[e1] -= pNE;//2.0*pNE-mag(pNE)*fNorms[faceI];
			}
		 }
		 pAreas=(mag(pNw));
		 Info<<"A: ["<<min(pAreas)<< "-"<<max(pAreas)<<"]="<<pAreas.avg(); cout.flush();

		 //pNs/=(mag(pNs)+1.0e-18);
		 pNw/=(mag(pNw)+1.0e-18);


	  }



	 double avgA=pAreas.avg();
	 vectorField fNorms(facezs.size());
	 {
	 for_(facezs, fI)  fNorms[fI]=0.01*normal(facezs[fI],newPoints);
	 labelListList pFacess = pointFaces(newPoints.size(), facezs);
	 for_(pFacess, pI) {
		const auto& pfacs=pFacess[pI];
		for_(pfacs, fI)
		{
			if(pMarks[pI]==facezs[fI].zone+1)
			{
				 fNorms[fI]+=pNw[pI];
			}
			else
				 pAreas[pI]=(0.9*avgA+0.1*pAreas[pI]);
		}
	 }//for_(pFaces, pI) 
	 fNorms/=(mag(fNorms)+1.0e-32);
	 }

	 { //smoothing



	  const vectorField previousPoints(newPoints);
	  Info<< iter   <<"  "; cout.flush();
	  
	  for (int iKern=0; iKern<kernelRadius; ++iKern)
	  {
		vectorField displac(newPoints.size(),dbl3(0.,0.,0.));
		dbls        sumWeis(newPoints.size(),1.0e-64);
		 double wt;
		 for_(facezs, fI)   //for(const auto& facezs:facezsZ) 
		 {	const auto& fac=facezs[fI];
			for_(fac, pI)
			{	int e0=fac[pI], e1=fac[(pI+1)%4];
				dbl3 delp=newPoints[e1]-newPoints[e0];
				if(pMarks[e1]>pMarks[e0])
				{
					wt=0.9;       displac[e0]+=wt*delp;	sumWeis[e0]+=wt;
					dbl3 CL=fNorms[fI]^pNw[e1]; CL/=mag(CL)+1.0e-18;
					wt=0.5;  displac[e1]-=wt*((delp&CL)*CL);	sumWeis[e1]+=wt*0.3;
				} else if(pMarks[e0]>pMarks[e1])
				{
					dbl3 CL=fNorms[fI]^pNw[e0]; CL/=mag(CL)+1.0e-18;
					wt=0.5;  displac[e0]+=wt*((delp&CL)*CL);	sumWeis[e0]+=wt*0.3;
					wt=0.9;       displac[e1]-=wt*delp;   sumWeis[e1]+=wt;
				} else {
					wt=1.0;       displac[e0]+=wt*delp;	sumWeis[e0]+=wt*(1.0-0.5*(pMarks[e0]!=fac.zone+1));
					wt=1.0;       displac[e1]-=wt*delp;	sumWeis[e1]+=wt*(1.0-0.5*(pMarks[e1]!=fac.zone+1));;
				}
			}
		 }
		 displac/=sumWeis;
		 displac*=relax*0.5;
		 newPoints += displac;
	 }


		vectorField dispAvgO(newPoints-previousPoints);
		vectorField dispAvg(newPoints-previousPoints);


		Info<< "  disp:  "<<mag(dispAvg).avg(); cout.flush();

		Info <<"  VolPresev dAvg:"; cout.flush();

	  for (int iKern=0; iKern<kernelRadius+2; ++iKern)
	  {
		const vectorField dispAvgTmp(dispAvg);
		//vectorField displac(newPoints.size(),dbl3(0.,0.,0.));
		//dbls        sumWeis(newPoints.size(),1.0e-64);
		dbls        sumWeis(pAreas.size(),1.0e-64);
		dispAvg *= sumWeis;
		 for_(facezs, fI)   //for(const auto& facezs:facezsZ) 
		 {	const auto& fac=facezs[fI];
			double wt;
			for_(fac, pI)
			{	int e0=fac[pI], e1=fac[(pI+1)%4];
				//dbl3 delp=dispAvgTmp[e1]-dispAvgTmp[e0];
				if(pMarks[e1]>pMarks[e0])
				{
					wt=0.01*pAreas[e1];  dispAvg[e0]+=wt*((dispAvgTmp[e1]&fNorms[fI])*fNorms[fI]);	sumWeis[e0]+=wt;
					wt=0.9*pAreas[e1];  dispAvg[e1]+=wt*((dispAvgTmp[e0]&fNorms[fI])*fNorms[fI]);	sumWeis[e1]+=wt;
				} else if(pMarks[e0]>pMarks[e1])
				{
					wt=0.9*pAreas[e0];  dispAvg[e0]+=wt*((dispAvgTmp[e1]&fNorms[fI])*fNorms[fI]);	sumWeis[e0]+=wt;
					wt=0.01*pAreas[e0];  dispAvg[e1]+=wt*((dispAvgTmp[e0]&fNorms[fI])*fNorms[fI]);	sumWeis[e1]+=wt;
				} else {
					wt=pAreas[e1];  dispAvg[e0]+=wt*dispAvgTmp[e1];	sumWeis[e0]+=wt;
					wt=pAreas[e0];  dispAvg[e1]+=wt*dispAvgTmp[e0];	sumWeis[e1]+=wt;
				}
			}
		 }
		 dispAvg/=sumWeis;
		 //dispAvg*=pWeights;
		 //displac*=relax*0.5;
		 //newPoints += displac;
	 }

		//dispAvg   -= 0.3*(dispAvg-(dispAvg& pNw)*pNw);	 // to allow larger weight below to avoid sliding near boundary	//dispAvg -= 0.3*(dispAvg-(dispAvg& pNs)*pNs);
		newPoints -= 1.07*dispAvg;     ///.  TODO change 0.5 0.5,  sum should be 1



		pointsAll = newPoints;





	  if(iter==nIters-1) {
			std::vector<std::vector<face> > facess = getOrderedFaces(facezsZ, ifcs);
			ints pValidIs(pointsAll.size(),-1); int pInd=-1;
			for(auto&faces:facess) for(auto&fac:faces) for(auto& pi:fac) if(pValidIs[pi]==-1) pValidIs[pi]=++pInd;
			for(auto&faces:facess) for(auto&fac:faces) for(auto& pi:fac) pi=pValidIs[pi];

			fstream vtk("dumpSurfSmooth"+_s(ifcs)+".vtk",ios::app);
			vtk<< "POINT_DATA " << pInd+1 << '\n'
			<< "FIELD attributes 3\n";

			dbls  pmrks(pInd+1);   for_i(pValidIs) if(pValidIs[i]!=-1) { pmrks[pValidIs[i]]=pMarks[i]; }
			vtk<< "\npmrks 1 " << pmrks.size() << " int\n" ;
			for_(pmrks,ii) vtk<<pmrks[ii]<<'\n';

			dbl3s dAvgs(pInd+1);   for_i(pValidIs) if(pValidIs[i]!=-1) { dAvgs[pValidIs[i]]=dispAvg[i]; }
			vtk<< "\ndAvgs 3 " << dAvgs.size() << " float\n" ;
			for_(dAvgs,ii) vtk<<dAvgs[ii]<<'\n';

			dbl3s dAvgO(pInd+1);   for_i(pValidIs) if(pValidIs[i]!=-1) { dAvgO[pValidIs[i]]=dispAvgO[i]; }
			vtk<< "\ndAvgO 3 " << dAvgO.size() << " float\n" ;
			for_(dAvgO,ii) vtk<<dAvgO[ii]<<'\n';
		}

		dispAvg=(newPoints-previousPoints);
		Info<< "  max:  "<<max(mag(dispAvg)); cout.flush();
		Info<< "  avgMag:  "<<mag(dispAvg).avg(); cout.flush();
		Info<< "  avg:  "<<dispAvg.avg(); cout.flush();
		Info<<endl;

	  Info<<"Gasss smoothed  :/"<<endl<<endl;

	  }
	}//===================================================================================



  }

	Info<<mag(pointsAll).avg()<< "  :avgMagPoints    "<<endl;
	Info<<max(mag(pointsAll))<< "  :maxMagPoints    "; cout.flush();


    Info << "End\n" << endl;

    return 0;
}


// ************************************************************************* //
