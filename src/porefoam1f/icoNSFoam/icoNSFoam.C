/*-------------------------------------------------------------------------*\
 Direct single-phase flow solver

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



#define SINGLE_PHASE
#define ifMonitor(TEN)  if (runTime.timeIndex()%TEN==0) 

#include "fvCFD.H"


//#include "singlePhaseTransportModel.H"
//#include "turbulentTransportModel.H"


#include "fixedFluxPressureFvPatchScalarField.H"

#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	#include "createMesh.H"
	if (!mesh.cells().size()) {Info<<"Error: no cells in (processor) mesh"<<endl; exit(-1);}
	pimpleControl pimple(mesh);
	#include "initContinuityErrs.H"
	#include "createFields.H"
	#include "createTimeControls.H"
	#include "correctPhi.H"
	#include "CourantNo.H"

	//#include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\nStarting time loop\n" << endl;
#define 	curtail(a,b) (min (max((a),(b)),(1.-(b))))

	ifMonitor(10)
	{
		Info<< "\n         Umax = " << max(mag(U)).value() << " m/s  "
		<< "Uavg = " << mag(average(U)).value() << " m/s"
		<< "   DP = " << (max(p)-min(p)).value() << " Pa"
		<< nl<< nl << endl;
	}

	Info <<"min(p): "<<min(p)<<"  max(p): "<<max(p)<<endl;



	scalar pRelaxF=0.1;




	while (runTime.run())
	{
		scalar timeStepExecutionTime= runTime.elapsedCpuTime() ;
		//#include "readPIMPLEControls.H"
		#include "readTimeControls.H"
		#include "CourantNo.H"
		#include "setDeltaT.H"

		runTime++;

		Info << nl<< "Time = " << runTime.timeName() << endl;



		dimensionedScalar sgPc("sgPc", dimPressure/dimLength*dimArea, 0.);
		while (pimple.loop())
		{
			//rhoPhi = rho1*phi;


			{
				volScalarField pOldCopy=p; 
				p=0.30*pOldOld+0.34*pOld+0.36*p;
				pOldOld=pOld;
				pOld=pOldCopy;
			}

			scalar tOldPU= runTime.elapsedCpuTime() ;



			#include "UEqn.H"

			ifMonitor(10)   { Info<< "ExeTime - tOldPU = " << runTime.elapsedCpuTime()-tOldPU << " s"	<< endl;}

			while (pimple.correct())// --- PISO loop
			{
				#include "pEqn.H"
			}

		}

		#include "continuityErrs.H"






		Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
			<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
			<< "  timeStepExecutionTime  = " << runTime.elapsedCpuTime()-timeStepExecutionTime << " s"
			<< endl;

		ifMonitor(10)
		{
			scalar maxU = max(mag(U)).value();
			scalar avgU = average(mag(U)).value();

			Info<< "\n         maxMagU = " <<maxU << " m/s  "
			<< "avgMagU = " << avgU << " m/s"
			<< "   DP = " << (max(p)-min(p)).value() << " Pa"
			<< nl<< nl << endl;
			
			scalar delUx10 =  mag(avgU - oldAvgU10);
			if (delUx10<refDelUx10*max(avgU,oldAvgU10) && oldDelUx10<refDelUx10*max(avgU,oldAvgU10))
			{
				Info<< "converged ! " 	<< endl;
				
				runTime.writeAndEnd();
			}
			Info<<"! convergence: "<<delUx10<<"<"<<refDelUx10*max(avgU,oldAvgU10) <<" && "
			    << oldDelUx10<<"<"<<refDelUx10*max(avgU,oldAvgU10)<<nl<<endl;
			oldDelUx10=delUx10;
			oldAvgU10=avgU;
		}




		scalar maxSimTime( readScalar(runTime.controlDict().lookup("maxExecutionTime")) );

		if (runTime.elapsedClockTime()>maxSimTime)
		{
			Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s\n"
				<< maxSimTime/24/60/60<<"  days passed!, ending simulation. " 
				<< endl;

			runTime.writeAndEnd();
		}

		#include "write.H"
	}

	if(!std::ifstream(runTime.timeName()+"/p").good())
	{
		Info<<"Error Pressure is not written, trying again:";
		U.write();
		phi.write();
		p.write();
	}

	mesh.clearPoints();// just to print the message
	mesh.clearNonOrtho();// just to print the message

    Info<< "End\n" << endl;

    return 0;
}


//----------------------------------------------------------------------
