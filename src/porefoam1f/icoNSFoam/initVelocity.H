
	Info << cBC<<endl<<endl<<endl;
	label nUCorr
	(
		readLabel(mesh.solutionDict().subDict("PIMPLE").lookup("nUCorr"))
	);
	dimensionedScalar beta ( mesh.solutionDict().subDict("PIMPLE").lookup("beta")  );
	
	scalar tOld= runTime.elapsedCpuTime() ;

        Info<<max(fvc::laplacian(p))<<min(fvc::laplacian(p))<<endl;

{	
        fvScalarMatrix pEqn
        (
            fvm::laplacian(p) 
        );
        pEqn.setReference(pRefCell, pRefValue);
        pEqn.solve(mesh.solver(p.name() + "Final"));
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
		fvc::grad(p)
	);	
	U=1.3*fvc::average(fvc::interpolate(U));
	p=fvc::average(fvc::interpolate(p));
	U.correctBoundaryConditions();
	p.correctBoundaryConditions();
	
	phi = 0.5*(fvc::interpolate(U) & mesh.Sf());

	//rhoPhi = rho1*phi;
	
	
	
	Info<<max(fvc::div(phi))<<endl;
	
}


