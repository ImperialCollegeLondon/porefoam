
	Info << cBC<<endl<<endl<<endl;
	label nUCorr
	(
		readLabel(mesh.solutionDict().subDict("PIMPLE").lookup("nUCorr"))
	);
	dimensionedScalar beta ( mesh.solutionDict().subDict("PIMPLE").lookup("beta")  );
	
	scalar tOld= runTime.elapsedCpuTime() ;

        Info<<max(fvc::laplacian(pd))<<"  "<<min(fvc::laplacian(pd))<<endl;

{	
        fvScalarMatrix pEqn
        (
            fvm::laplacian(pd) 
        );
        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve(mesh.solver(pd.name() + "Final"));
}
{	
        fvScalarMatrix pEqn
        (
            fvm::laplacian(pd) 
        );
        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve(mesh.solver(pd.name() + "Final"));
}


    
{
    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
        //+ fvc::interpolate(rho*turbulence->nut()) //caution turbulance is disabled
    );
    solve
	(
		fvm::laplacian(muEff, U)
		==
		fvc::grad(pd)
	);
	
	
	
	#include  "correctmuEff.H"
	solve
	(
		fvm::laplacian(muEff, U)
		==
		fvc::grad(pd)
	);
}
{
	
    surfaceScalarField muEff
    (
        "muEff",
        twoPhaseProperties.muf()
        //+ fvc::interpolate(rho*turbulence->nut()) //caution turbulance is disabled
    );
	
	#include  "correctmuEff.H"
	solve
	(
		fvm::laplacian(muEff, U)
		==
		fvc::grad(pd)
	);	
	U*=0.5;
	phi = 0.5*(fvc::interpolate(U) & mesh.Sf());

	//rhoPhi = rho1*phi;
	
	Info<<max(fvc::div(phi))<<endl;
	
}


