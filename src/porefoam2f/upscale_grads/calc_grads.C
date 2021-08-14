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
//!   post-processing code Calculates energy losses for control 
//!   volumes given in postprossingDict file in $case/system directory. 
//!   
//!   Note: Parallel case post-processing is added to the code, no need to 
//!   reconstract parallel cases anymore.


#include <fstream>
#include <string>
#include <assert.h>
#include <vector>
#include "fvCFD.H"
 
 
#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "twoPhaseMixture.H"
#include "interfaceProps.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "interpolationCellPoint.H"
#include "primitivePatchInterpolation.H"
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
					IOobject::MUST_READ
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
					IOobject::MUST_READ
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
				  IOobject::READ_IF_PRESENT
				),
				  fvc::reconstruct
				  (surfaceScalarField 
					(///. to delete, chnage back to must_read
					  IOobject
					  (
						  "phi",
						  runTime.timeName(),
						  mesh,
						  IOobject::READ_IF_PRESENT
					  ),
					 (linearInterpolate(U) & mesh.Sf())*0.
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
						  IOobject::READ_IF_PRESENT
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
						IOobject::READ_IF_PRESENT
					),
					0.*fvc::grad(pc)
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



