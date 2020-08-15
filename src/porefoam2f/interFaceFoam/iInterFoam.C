/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of porefoam.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your (at your (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

	-----------------------------------------------------------------


The code has been developed by Mosayeb Shams and Ali Qaseminejad Raeini
, under the supervision of Branko Bijeljic and Martin Blunt. 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

	For further information please contact:
	Mosayeb Shams:   m.shams14@imperial.ac.uk
	Ali Q Raeini:    a.q.raeini@imperial.ac.uk
	Branko Bijeljic: b.bijeljic@imperial.ac.uk
	Martin J Blunt:  m.blunt@imperial.ac.uk

Description:
    

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "interfaceProperties.H"
#include "twoPhaseMixture.H"
//#include "turbulenceModel.H"
#include "pimpleControl.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "pressureDirectedInletVelocityFvPatchVectorField.H"

#define _POSTPROCESS_
#ifdef _POSTPROCESS_
 #include "voxelImage.h"
 #include "AverageData.h"
 #include "myFVC.H"

 #define _USE_MPI_
 #include <mpi.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define ifMonitor  if( runTime.timeIndex()%10== 0 ) 
#define curtailBADOFSET(a,b) (min (max((a),(b)),(1.-(b))))

int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	//pimpleControl pimple(mesh);
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#ifdef _POSTPROCESS_
	  #include "createCVs.H"
	#endif
	#include "createTimeControls.H"
	#include "correctPhi.H"
	#include "CourantNo.H"
	//#include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;


	scalar tOld= runTime.elapsedCpuTime() ;


	while (runTime.run())
	{
		scalar timeStepExecutionTime= runTime.elapsedCpuTime() ;
		#include "readPIMPLEControls.H"
		#include "readTimeControls.H"
		#include "CourantNo.H"
		#include "alphaCourantNo.H"
		#include "setDeltaT.H"

		runTime++;

		Info<< "Time = " << runTime.timeName() << nl << endl;



		dimensionedScalar totalDeltaT = runTime.deltaT();


		int oCorr=nOuterCorr-1;
		subCycleTime alphaSubCycle(runTime, 2);
		{
			tOld= runTime.elapsedCpuTime() ;
			interface.correct(0.0);
			sgPc = interface.sgPc();
			ifMonitor  {Info<< "  ExeTime pc = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;}
			twoPhaseProperties.correct();

			#include "UEqn.H"
			for (int corr=0; corr<nCorr; corr++)// --- PISO loop
			{
				#include "pEqn.H"
			}
			#include "alphaEqn.H"
		}
		alphaSubCycle++;



		ifMonitor {Info<< "ExeTime t half cycle = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;}

		for (oCorr=0; oCorr<nOuterCorr; oCorr++)
		{

			{
				tOld= runTime.elapsedCpuTime() ;
				interface.correct(0.5);
				ifMonitor  {Info<< oCorr <<"ExeTime pc = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;}

				twoPhaseProperties.correct();
				if (oCorr==nOuterCorr-1)  sgPc = interface.sgPc();
				else                     sgPc = 0.6*sgPc+0.4*interface.sgPc() ;
			}


			tOld= runTime.elapsedCpuTime() ;

			#include "UEqn.H"
			for (int corr=0; corr<nCorr; corr++)// --- PISO loop
			{
				#include "pEqn.H"
			}
			#include "alphaEqn.H"
			ifMonitor   { Info<< "ExeTime pd-u = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;}

		}



		alphaSubCycle++;
		alphaSubCycle.endSubCycle();


		#include "continuityErrs.H"






		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< "  timeStepExecutionTime  = " << runTime.elapsedCpuTime()-timeStepExecutionTime << " s"
			<< endl;

		ifMonitor
		{
			Info<< "\n\n Time = " << runTime.timeName()  << " s"
				<< "   S_a = "
				<< alpha1.weightedAverage(mesh.V()).value()
				<< "   Pc = " << (average(interface.pc()*alpha1)/(mag(average(alpha1))+1e-16)-average(interface.pc()*(1-alpha1))/(mag(average(1-alpha1))+1e-16)).value() << " Pa "
				<< "   Uavg = " << mag(average(U)).value() << " m/s"
				<< "   Umax = " << max(mag(U)).value() << " m/s"
				<< "   D_alpha = " << max(alpha1).value() - min(alpha1).value()
				<< "   D_p = " << (max(pd)-min(pd)).value() << " Pa"
				<< nl<< nl << endl;
		}

		#ifdef _POSTPROCESS_
		if( runTime.timeIndex()%50== 0 )
		{
			const volScalarField& pc=interface.pc();
			const volVectorField& gPc=interface.gPc();
			#include "calc_grads.H"
		}
		#endif



		scalar maxSimTime( readScalar(runTime.controlDict().lookup("maxExecutionTime")) );
		scalar maxS1 (  runTime.controlDict().lookupOrDefault("maxS1", scalar(1.1)) );
		scalar minS1  (  runTime.controlDict().lookupOrDefault("minS1",scalar(-0.1)) );
		scalar aAvg = average(alpha1).value();
		if (runTime.elapsedClockTime()>maxSimTime || aAvg>=maxS1 || aAvg<=minS1)
		{
			Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
				<< "  average(alpha) = " << aAvg << "  "
				<< "  maxExecutionTime ("<<maxSimTime<<"s) reached,   " 
				<< "  or maxS1 ("<<maxS1<<") reached,   " 
				<< "  or minS1 ("<<minS1<<") reached, ending simulation. " 
				<< "  or S1 < " << minS1 << "  or " << maxS1 << "<  S1"
				<< endl;

			runTime.writeNow();
			return 0;
		}

		#include "write.H"
	}

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
