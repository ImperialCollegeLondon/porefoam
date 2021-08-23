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
//!   smooths mesh to improve its quality, as reported by checkMeshh


#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "transformGeometricField.H"
#include "IOobjectList.H"


#include "syncTools.H"

using namespace Foam;




// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	timeSelector::addOptions();

	argList::validArgs.append("relax");
	argList::validArgs.append("nIter");
	argList::validArgs.append("strictOnBoundary");
	argList::validArgs.append("boundaryAttractionFactor");

#   include "setRootCase.H"
#   include "createTime.H"

	scalar relax = args.argRead<scalar>(1);

	label nIter = args.argRead<label>(2);
	label strictOnBoundary = args.argRead<label>(3);
	scalar boundaryAttractionFactor = args.argRead<scalar>(4);

	Info<<"  relax:"<<relax <<"  nIter:"<<nIter  
	    <<"  strictOnBoundary:"<<strictOnBoundary  <<"  boundaryAttractionFactor:"<<boundaryAttractionFactor <<endl;

	pointIOField points (
		IOobject (
			 "points",
			 runTime.findInstance(polyMesh::meshSubDir, "points"),
			 polyMesh::meshSubDir,
			 runTime,
			 IOobject::MUST_READ,
			 IOobject::NO_WRITE,
			 false
		)
	);





#   include "createMesh.H"




	scalarField pointWeis(points.size(), 1.);
	labelList   pointOnBondry(points.size(), 0);
	{
		
		int nPoints = points.size();
		const labelListList&  pointfaces = mesh.pointFaces ();
		const labelListList&  pointPoints = mesh.pointPoints ();
		const label        nInteFaces = mesh.nInternalFaces();

		//OMPFor()
		for(int pI=0; pI<nPoints; ++pI)  {
			const auto& facs=pointfaces[pI];
			forAll(facs,j)  {
				label fJ=facs[j];
				if(fJ>=nInteFaces) {
					pointOnBondry[pI] = 1;  pointWeis[pI]=boundaryAttractionFactor; } ///SYNC11
			}
		}

		for(int i=0;i<3;++i){// smooth pointWeis  the results are sensitive to nIter=3 used here
			scalarField pointWeisTmp = pointWeis;
			//OMPFor()
			for(int pI=0; pI<nPoints; ++pI)
			{
				scalar PwBSum(1.*pointWeis[pI]);
				scalar wPwBSum(1.);
				const labelList& 	neiPoints = pointPoints[pI];
				forAll(neiPoints, ppI) 
				{
				 label pJ = neiPoints[ppI];
				 //if(pointOnBondry[pJ])
				 {
					wPwBSum += 1.; 
					PwBSum += pointWeisTmp[pJ]; 
				 }
				}
				//if(pointOnBondry[pI]) 	pointWeis[pI]=0.9*PwBSum/wPwBSum;///SYNC11
				//else	   				
					pointWeis[pI]=PwBSum/wPwBSum;
			}
		}

	}






	for (int i=0; i<nIter;++i)
	{

		const labelList& 	    owns = mesh.owner();
		const labelList& 	    neis = mesh.neighbour();
		const labelListList&  pointfaces = mesh.pointFaces ();
		const labelListList&  pointPoints = mesh.pointPoints ();
		const pointField      pointsOrig = points;
		const vectorField&    cellCentres = mesh.cellCentres();    
		vectorField           faceNormals = mesh.faceAreas();  
		const scalarField     magAfs = mag(faceNormals)+1e-32; 
		faceNormals/=magAfs;
		vectorField           pNorms(points.size(), vector::zero);

		const scalar       VavgInv = 1./gAverage(mesh.cellVolumes()); 
		const scalarField  Vols = max(mesh.cellVolumes(),0.1/VavgInv);
		const vectorField& Cfs = mesh.faceCentres(); 
		const label        nInteFaces = mesh.nInternalFaces();
		vectorField        avgPCellCntrs(points.size(), vector::zero);
		vectorField        avgPPoints(points.size(), vector::zero);

		vectorField avgPFCentres(points.size(), vector::zero);
		
		int nPoints = points.size();

		//OMPFor()
		for(int pI=0; pI<nPoints; ++pI)
		{
			const auto& facs=pointfaces[pI];
			vector avgFcntr(points[pI]*0.001);  double sumWfc=0.001;
			vector pointNrm(vector::zero); int sumWnr=0;
			forAll(facs,j)  {
				label fJ=facs[j];
				if(fJ>=nInteFaces) {
					pointNrm+= faceNormals[fJ]; ++sumWnr;}
				
				scalar wei = max(0.2, VavgInv*(fJ<nInteFaces ?  
												(Vols[owns[fJ]]+Vols[neis[fJ]]) : 6.*Vols[owns[fJ]])); 
				avgFcntr += wei*Cfs[fJ]; 
				sumWfc+=wei;
			}
			if(sumWnr) {
				avgFcntr=(points[pI]*0.001);  sumWfc=0.001;
				forAll(facs,j)
				{
					label fJ=facs[j];
					if(fJ>=nInteFaces) {
						scalar wei = max(0.2, 6.*Vols[owns[fJ]]); 
						avgFcntr += wei*Cfs[fJ]; 
						sumWfc+= wei; }
				}
				vector disp=avgFcntr/sumWfc-points[pI];
				pointNrm/=sumWnr;
				pNorms[pI]=pointNrm;
				avgPFCentres[pI]=points[pI]+(1+relax)*disp-(disp&pointNrm)*pointNrm;
			}
			else 
				avgPFCentres[pI]=avgFcntr/sumWfc;
		}

		//OMPFor()
		for(int pI=0; pI<nPoints; ++pI)
		{
			const labelList& 	pointCells = mesh.pointCells(pI);
			vector CCSum(vector::zero);
			scalar WSum(1e-64);
			  forAll(pointCells, pcI)
			  {
					label cI = pointCells[pcI];
					scalar w=(Vols[cI]);
					//w*=1e15*w;
					CCSum += w*cellCentres[cI]; 
					WSum+=w;
			  }
			  avgPCellCntrs[pI]=CCSum/WSum;
		}


		{ //avgPPoints

			//OMPFor()
			for(int pI=0; pI<nPoints; ++pI)
			{
				const labelList& 	neiPoints = pointPoints[pI];
				vector CCSum(vector::zero);
				scalar WSum(1e-64);
				label  pIOnB=pointOnBondry[pI];
				forAll(neiPoints, ppI)
				{
					label pII = neiPoints[ppI];
					scalar w=pointWeis[pII];//1.;//(pointsOrig[pI]-pointsOrig[pII]);
					if(pIOnB && pointOnBondry[pII]!=pIOnB)  w = 1e-64;

					CCSum += w*pointsOrig[pII]; 
					WSum+=w;
				}
				avgPPoints[pI]=CCSum/WSum;
			}
		}



		 //syncTools::syncPointList ( mesh, avgPCellCntrs, plusEqOp<point>(), vector::zero   );
		 //syncTools::syncPointList ( mesh, avgPPoints, plusEqOp<point>(),  vector::zero   );

		//OMPFor()
		for(int pI=0; pI<nPoints; ++pI)
		{
			if (!pointOnBondry[pI])
			  points[pI]= relax*(0.4*avgPCellCntrs[pI] + 0.3*avgPPoints[pI] +0.3*avgPFCentres[pI]) 
			        + (1-relax)*pointsOrig[pI];
			else if(!(strictOnBoundary&4))
			{
				vector disp=relax*(0.4*avgPCellCntrs[pI] + 0.3*avgPPoints[pI] +0.3*avgPFCentres[pI] - pointsOrig[pI]);

				disp=disp-(disp&pNorms[pI])*pNorms[pI];
				
				if (!(strictOnBoundary&2)){
						disp+=0.2*relax*(avgPFCentres[pI] - pointsOrig[pI]);
				}

				if (strictOnBoundary&1) {
					forAll(pointfaces[pI],j)
					{
						if (pointfaces[pI][j]>=nInteFaces)
						{
							vector norm = faceNormals[pointfaces[pI][j]];
							disp=disp-(disp&norm)*norm;
						}
					}
					forAll(pointfaces[pI],j)
					{
						if (pointfaces[pI][j]>=nInteFaces)
						{
							vector norm = faceNormals[pointfaces[pI][j]];
							disp=disp-(disp&norm)*norm;
						}
					}
				}



				points[pI]+=disp;
			}
		 }







		vectorField  displacements(points-pointsOrig);
		//Info<< " -> "<<max(mag(displacements)); cout.flush();

		for(int iter3=0;iter3<3;iter3++)
		{
			  vectorField dispTmp(displacements);
			  forAll(pointOnBondry, pI)
			  {
					if (pointOnBondry[pI])
					{
						vector avgDispV(displacements[pI]);
						scalar sumWeights(1.0000000001);

						const labelList& neiPoints = pointPoints[pI];
						forAll(neiPoints, neipI)
						{
							 scalar weight=0.1+(pNorms[neiPoints[neipI]] & pNorms[pI]);
							 weight*=weight; weight*=weight;
							 avgDispV += weight * (dispTmp[neiPoints[neipI]] & pNorms[neiPoints[neipI]])*pNorms[neiPoints[neipI]];
							 sumWeights += weight;
						}
						displacements[pI]=avgDispV/sumWeights;
					}
				}
		 }
		 
		forAll(pointOnBondry, pI)
			if (pointOnBondry[pI])    points[pI] -= 1.001*(displacements[pI]& pNorms[pI]) * pNorms[pI];



		Info<< " -> "<<max(mag(displacements))<<"  "; cout.flush();











		mesh.clearOut();
	}



	// Set the precision of the points data to 10
	IOstream::defaultPrecision(max(10u, IOstream::defaultPrecision()));

	Info<< "Writing points into directory " << points.path() << nl << endl;
	points.write();
	instantList timeDirs = timeSelector::select0(runTime, args);


	Info<< "\nEnd.\n" << endl;

	return 0;
}


// ************************************************************************* //
