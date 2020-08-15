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

#   include "setRootCase.H"
#   include "createTime.H"

	scalar relax = args.argRead<scalar>(1);

	label nIter = args.argRead<label>(2);
	bool strictOnBoundary = args.argRead<bool>(3);

	Info<<"  relax:"<<relax <<"  nIter:"<<nIter <<endl;
	//tensor T = rotationTensor(n1, n2);

	pointIOField points
	(
		IOobject
		(
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











	for (int i=0; i<nIter;++i)
	{

		const labelList& 	owners = mesh.owner();
		const labelList& 	neibrs = mesh.neighbour();
		const labelListList& 	pointfaces = mesh.pointFaces ();
		const labelListList& 	pointPoints = mesh.pointPoints ();
		const pointField pointsOrig = points;
		const vectorField& cellCentres = mesh.cellCentres();    
		vectorField faceNormals = mesh.faceAreas();  
		const scalarField magAfs = mag(faceNormals)+1.0e-32; 
		faceNormals/=magAfs;
		vectorField pNorms(points.size(), vector::zero);

		const scalarField& Vols = mesh.cellVolumes(); 
		const scalar& VavgInv = 1.0/gAverage(Vols); 
		const vectorField& Cfs = mesh.faceCentres(); 
		const label nInteFaces = mesh.nInternalFaces();
		vectorField avgPCellCntrs(points.size(), vector::zero);
		vectorField avgPPoints(points.size(), vector::zero);

		scalarField pointWeis(points.size(), 1.0);
		labelList pointOnBoundary(points.size(), 1);
		vectorField avgPFCentres(points.size(), vector::zero);
		{
			scalarField nPointFaces(points.size(), 1.0e-64);
			forAll(Cfs, fI)
			{
				const face& 	f = mesh.faces()[fI];
				vector CCSum(Cfs[fI]);
				scalar wei = VavgInv*(fI<nInteFaces ?  //Warning 1.0e15 should be replaced by 1/AvgVolCell
												(Vols[owners[fI]]+Vols[neibrs[fI]]) : 
												6.0*Vols[owners[fI]]);//*(magAfs[fI])
				forAll(f, fpI)
				{
						label pI = f[fpI];
						avgPFCentres[pI] += wei*CCSum; 
						nPointFaces[pI]+=wei;
						pointOnBoundary[pI]=fI<nInteFaces ? 0 : 1; // WAS ALWAYS 1
						pointWeis[pI]=fI<nInteFaces ? 1.0 : 27.0;
				}
			}
			avgPFCentres/=nPointFaces;
		}


		forAll(pointsOrig, pI)
		{
			  const labelList& 	pointCells = mesh.pointCells(pI);
			vector CCSum(vector::zero);
			scalar WSum(1.0e-64);
			  forAll(pointCells, pcI)
			  {
					label cI = pointCells[pcI];
					scalar w=(Vols[cI]);
					//w*=1.0e15*w;
					CCSum += w*cellCentres[cI]; 
					WSum+=w;
			  }
			  avgPCellCntrs[pI]=CCSum/WSum;
		}

		{
			scalarField pointWeisTmp = pointWeis;
			forAll(pointWeisTmp, pI)
			{
				scalar PwBSum(0.0);
				scalar wPwBSum(0.0);
				const labelList& 	neiPoints = pointPoints[pI];
				forAll(neiPoints, ppI) 
				{
				 label pII = neiPoints[ppI];
				 if(pointOnBoundary[pII])
				 {
					wPwBSum += 1.0; 
					PwBSum += pointWeisTmp[pII]; 
				 }
				}
				if(!pointOnBoundary[pI])	pointWeis[pI]=(PwBSum+2.0*pointWeis[pI])/(2.0+wPwBSum);
			}
		}

		{
			forAll(pointsOrig, pI)
			{
				const labelList& 	neiPoints = pointPoints[pI];
				vector CCSum(vector::zero);
				scalar WSum(1.0e-64);
				label  pIOnB=pointOnBoundary[pI];
				forAll(neiPoints, ppI)
				{
					label pII = neiPoints[ppI];
					scalar w=pointWeis[pII];//1.0;//(pointsOrig[pI]-pointsOrig[pII]);
					if(pIOnB && pointOnBoundary[pII]!=pIOnB)  w = 1.0e-64;

					CCSum += w*pointsOrig[pII]; 
					WSum+=w;
				}
				avgPPoints[pI]=CCSum/WSum;
			}
		}



		 //syncTools::syncPointList
		 //( mesh, avgPCellCntrs,
			  //plusEqOp<point>(),  // combine op
			  //vector::zero        // null value
		 //);
		 //syncTools::syncPointList
		 //( mesh, avgPPoints,
				//plusEqOp<point>(),  // combine op
			  //vector::zero        // null value
		 //);

		 forAll(points, pI)
		 {
			if (!pointOnBoundary[pI])
			  points[pI]= relax*(0.4*avgPCellCntrs[pI] + 0.3*avgPPoints[pI] +0.3*avgPFCentres[pI]) 
			        + (1-relax)*pointsOrig[pI];
			else
			{
				vector disp=relax*(0.5*avgPCellCntrs[pI] + 0.2*avgPPoints[pI] +0.3*avgPFCentres[pI] - pointsOrig[pI])
				;


				vector normbfs(0.0,0.0,0.0);
				forAll(pointfaces[pI],j)
				{
					if (pointfaces[pI][j]>=nInteFaces)
					{
						normbfs += faceNormals[pointfaces[pI][j]];
					}
				}
				normbfs/=mag(normbfs)+1.0e-28;
				pNorms[pI] = normbfs;
				disp=disp-(disp&normbfs)*normbfs;
				
				disp+=0.1*relax*((avgPPoints[pI] - pointsOrig[pI])&normbfs)*normbfs;


				if (strictOnBoundary)
				{
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

		for(int iter3=0;iter3<2;iter3++)
		{
			  vectorField dispTmp(displacements);
			  forAll(pointOnBoundary, pI)
			  {
					if (pointOnBoundary[pI])
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
		 
		forAll(pointOnBoundary, pI)
			if (pointOnBoundary[pI])    points[pI] -= 1.0625*(displacements[pI]& pNorms[pI]) * pNorms[pI];



		Info<< " -> "<<max(mag(displacements)); cout.flush();











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
