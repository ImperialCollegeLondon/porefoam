/*-------------------------------------------------------------------------*\
This code is part of poreFOAM, a suite of codes written using OpenFOAM
for direct simulation of flow at the pore scale. 	
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.



The code has been developed by Ali Qaseminejad Raeini as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
 
 Description:
	post-processing code Calculates energy losses for control 
	volumes given in postprossingDict file in $case/system directory. 
	
	Note: Parallel case post-processing is added to the code, no need to 
	reconstract parallel cases anymore.
\*-------------------------------------------------------------------------*/


#include <fstream>
#include <string>
#include <assert.h>
#include <vector>
#include "fvCFD.H"
 
 
#include "argList.H"
#include "interpolationCellPoint.H"
#include "primitivePatchInterpolation.H"
//#include "sampledPlane.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "twoPhaseMixture.H"
#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"

 #include "voxelImage.h"
 #include "AverageData.h"
 #include "myFVC.H"
 #define _USE_MPI_
 #include <mpi.h>



using namespace Foam;





int main(int argc, char *argv[])
{


    timeSelector::addOptions();
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"
 
	runTime.setTime(timeDirs[0], 0);


	#include  "createFields.H"


	runTime.setTime(timeDirs[0], 0);
	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);
	runTime.setTime(timeDirs[0], 0);

	//if (readData && 0)
	//{
		//tzData tzDatas(timeDirs.size(),CVBounds1.size()-1);
		//tzDatas.read();
		//tzDatas.write("readDataWriteNOTIMPLEMENTED");
		//return(0);
	//}




	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI);
		Info<< endl<<timeI<< "    Time = " << runTime.timeName() << "''''''VVVVVVVV'''' ";		
			volScalarField pc
			(
				IOobject
				(
					"pc",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);
 
			pd==volScalarField 
			(
				IOobject
				(
					"pd",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

			alpha1== volScalarField 
			(
				IOobject
				(
					"alpha1",
					runTime.timeName(),
					mesh,
					IOobject::MUST_READ,
					IOobject::NO_WRITE
				),
				mesh
			);

		  //alphaContactAngleFvPatchScalarField::reset = true;
		  //alpha1.correctBoundaryConditions();
			U== volVectorField
			(	IOobject
				(	"U",
					runTime.timeName(),
					mesh,
				  IOobject::READ_IF_PRESENT,
				  IOobject::NO_WRITE
				),
				  fvc::reconstruct
				  (surfaceScalarField 
					(///. to delete, chnage back to must_read
					  IOobject
					  (
						  "phi",
						  runTime.timeName(),
						  mesh,
						  IOobject::READ_IF_PRESENT,
						  IOobject::NO_WRITE
					  ),
					 (linearInterpolate(U) & mesh.Sf())*0.0
					)
				  )
			);

				phi ==  surfaceScalarField 
				  (
					  IOobject
					  (
						  "phi",
						  runTime.timeName(),
						  mesh,
						  IOobject::READ_IF_PRESENT,
						  IOobject::NO_WRITE
					  ),
					  linearInterpolate(U) & mesh.Sf()
				  );

//TODO test if necassary:
				#include "calcPhis.H"

				volVectorField gPc
				(	IOobject
					(	"gPc",
						runTime.timeName(),
						mesh,
						IOobject::READ_IF_PRESENT,
						IOobject::NO_WRITE
					),
					0.0*fvc::grad(pc)
				);	

		  	  //surfaceScalarField sgPc   
			  //(
				  //IOobject
				  //(
					  //"sgPc",
					  //runTime.timeName(),
					  //mesh,
					  //IOobject::READ_IF_PRESENT,
					  //IOobject::NO_WRITE
				  //),
				  //linearInterpolate(gPc) & mesh.Sf()
			  //);

			surfaceScalarField muEff
			(
				"muEff",
				twoPhaseProperties.muf()
			);

			gradP=fvc::reconstruct
					(
						(
						  fvc::snGrad(pd)
						) * mesh.magSf()
					);

		#define 	curtailBADOFSET(a,b) (min (max((a),(b)),(1.-(b))))


		#include "calc_grads.H"
	}

    Info<< "end" << endl<< endl<< endl<< endl;


    return 0;
}



