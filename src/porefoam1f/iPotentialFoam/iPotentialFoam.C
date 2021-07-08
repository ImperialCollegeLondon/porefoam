/*-------------------------------------------------------------------------*\
 Initialize flow field, speeds up convergence in some cases

 Copyright (C) 2012-2020  Ali Qaseminejad Raeini 

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
#define ifMonitor  if( runTime.timeIndex()%10== 0 ) 

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
	//#include "correctPhi.H" // GAMGPCG:  Solving for p:  solution singularity !
	#include "CourantNo.H"

	//#include "setInitialDeltaT.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	runTime.run();


	Info<< "\nStarting time loop\n" << endl;
#define 	curtail(a,b) (min (max((a),(b)),(1.-(b))))


	Info <<"min(p): "<<min(p)<<"max(p): "<<max(p)<<endl;
	scalar pRelaxF=0.1;




	Info << "cBC "<< cBC<<endl<<endl<<endl;

	scalar tOld= runTime.elapsedCpuTime() ;


/*{
	surfaceScalarField muEff
	(	IOobject
		( "muEff", runTime.timeName(), mesh
		), mesh, rho*nu
	);
	tmp<fvVectorMatrix> UEqn 
	(
		fvm::div(rho*phi, U)
	  - fvm::laplacian(muEff, U)
	);

	volScalarField rAU = 1./UEqn().A();
	surfaceScalarField rAUf = fvc::interpolate(rAU);
	setSnGrad<fixedFluxPressureFvPatchScalarField> ( p.boundaryFieldRef(),
		( phi.boundaryField()) / (mesh.magSf().boundaryField()*rAUf.boundaryField())	);
}*/

{
	fvScalarMatrix pEqn
	(
		fvm::laplacian(p) 
	);
	pEqn.setReference(pRefCell, pRefValue);
	pEqn.solve(mesh.solutionDict().solver("pcorr"));
}




{
	fvScalarMatrix pEqn
	(
		fvm::laplacian(1e6/((fvc::grad(p)&fvc::grad(p))+0.001*average(fvc::grad(p)&fvc::grad(p))),p)
	);
	pEqn.setReference(pRefCell, pRefValue);
	pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
}

	
	Info<< "\n         Umax = " << max(mag(U)).value() << " m/s  "
	<< "Uavg = " << mag(average(U)).value() << " m/s"
	<< "   DP = " << (max(p)-min(p)).value() << " Pa"
	<< nl<< nl << endl;



{
	
	surfaceScalarField muEff
	(	IOobject ( "muEff", runTime.timeName(), mesh ),
		mesh, rho*nu
		//+ fvc::interpolate(rho*turbulence->nut()) //caution turbulance is disabled
	);



	//#include  "correctmuEff.H"
	solve
	(
		fvm::laplacian(muEff, U) - fvm::div(rho*phi, U) /*  - fvm::div(rho*phi, U) is just to avoid OpenFOAM bug with solver selection*/
		==
		fvc::grad(p)
	);	
	phi = (fvc::interpolate(U) & mesh.Sf());

	solve
	(
		fvm::laplacian(muEff, U) - fvm::div(rho*phi, U)
		==
		fvc::grad(p)
	);	
	phi = (fvc::interpolate(U) & mesh.Sf());
	
	Info<< "\n         Umax = " << max(mag(U)).value() << " m/s  "
	<< "Uavg = " << mag(average(U)).value() << " m/s"
	<< "   DP = " << (max(p)-min(p)).value() << " Pa"
	<< nl<< nl << endl;

	
	//#include  "correctmuEff.H"
	solve
	(
		fvm::laplacian(muEff, U) - fvm::div(rho*phi, U)
		==
		fvc::grad(p)
	);	
	phi = (fvc::interpolate(U) & mesh.Sf());
}	


	Info<< "\n         Umax = " << max(mag(U)).value() << " m/s  "
	<< "Uavg = " << mag(average(U)).value() << " m/s"
	<< "   DP = " << (max(p)-min(p)).value() << " Pa"
	<< nl<< nl << endl;

	Info<<endl<<endl<<endl<<max(fvc::div(phi))<<endl<<endl;

for (int ii=1;ii<10;++ii)
{
	surfaceScalarField muEff
	(	IOobject ( "muEff", runTime.timeName(), mesh ),
		mesh, rho*nu
	);
	#include  "correctmuEff.H"
	
	tmp<fvVectorMatrix> UEqn 
	(
		fvm::div(rho*phi, U) - fvm::laplacian(muEff, U)
	);
	scalar rlx=min((ii-0.9)/8.,1.);
	volScalarField rAU = 1./UEqn().A();
	surfaceScalarField rAUf = fvc::interpolate(rAU);
	//rAUf = rlx*rAUf+(1.-rlx)*fvc::interpolate(fvc::average(rAUf));
	U = rAU*(UEqn().H());
	phi = (fvc::interpolate(U) & mesh.Sf());
	UEqn.clear();
	//setSnGrad<fixedFluxPressureFvPatchScalarField> ( p.boundaryFieldRef(),
	//	( phi.boundaryField()  // - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
	//	) / (mesh.magSf().boundaryField()*rAUf.boundaryField())	);

	for(int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
	{ 
		fvScalarMatrix pEqn( fvm::laplacian(rAUf, p) 
			  == rlx*fvc::div(phi) + (1.-rlx)*fvc::laplacian(rAUf,p) );
		pEqn.setReference(pRefCell, pRefValue);
		if (nonOrth==pimple.nNonOrthCorr())
			pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
		else
			pEqn.solve(mesh.solutionDict().solver(p.name()));

		if (nonOrth == pimple.nNonOrthCorr())
			phi -= pEqn.flux();
	}
	U -= rAU*(fvc::grad(p));
	U.correctBoundaryConditions();




	Info<< "\n         Umax = " << max(mag(U)).value() << " m/s  "
	<< "Uavg = " << mag(average(U)).value() << " m/s"
	<< "   DP = " << (max(p)-min(p)).value() << " Pa"
	<< nl<< nl << endl;
}


	Info<<endl<<max(fvc::div(phi))<<endl;

	for (int ii=1;ii<=2;++ii)
	{
		surfaceScalarField muEff
		(	IOobject ( "muEff", runTime.timeName(), mesh ),
			mesh, rho*nu	//+ fvc::interpolate(rho*turbulence->nut()) //caution turbulance is disabled
		);
		#include  "correctmuEff.H"

		tmp<fvVectorMatrix> UEqn 
		(
			fvm::div(rho*phi, U) - fvm::laplacian(muEff, U)
		);

		volScalarField rAU = 1./UEqn().A();
		surfaceScalarField rAUf = fvc::interpolate(rAU);
		U = rAU*(UEqn().H());
		phi = (fvc::interpolate(U) & mesh.Sf());
		UEqn.clear();

		//setSnGrad<fixedFluxPressureFvPatchScalarField> ( p.boundaryFieldRef(),
		//	( phi.boundaryField()  // - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
		//	) / (mesh.magSf().boundaryField()*rAUf.boundaryField())	);

		for(int nonOrth=0; nonOrth<=pimple.nNonOrthCorr(); nonOrth++)
		{
			fvScalarMatrix pEqn( fvm::laplacian(rAUf, p) == fvc::div(phi) );
			pEqn.setReference(pRefCell, pRefValue);
			if (nonOrth==pimple.nNonOrthCorr())
				pEqn.solve(mesh.solutionDict().solver(p.name() + "Final"));
			else
				pEqn.solve(mesh.solutionDict().solver(p.name()));

			if (nonOrth == pimple.nNonOrthCorr())
				phi -= pEqn.flux();
		}
		U -= rAU*(fvc::grad(p));
		U.correctBoundaryConditions();

	
		Info<< "\n    Umax = " << max(mag(U)).value() << " m/s  "
			<< "Uavg = " << mag(average(U)).value() << " m/s"
			<< "   DP = " << (max(p)-min(p)).value() << " Pa"
			<< nl<< nl << endl;
	}

	Info<<endl<<max(fvc::div(phi))<<endl;
	




if(max(p).value()<100000000.)
{

	U.write();
	phi.write();

	p.write();

}else Info<<"error too high p_max, exiting"<<endl;





Info<< "End\n" << endl;

return 0;
}


// ************************************************************************* //
