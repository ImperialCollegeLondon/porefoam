	surfaceScalarField muEff
	(	IOobject
		(	"muEff",  runTime.timeName(),  mesh,
		IOobject::NO_READ,  IOobject::NO_WRITE
		),  mesh,  rho*nu
	);
	#include  "correctmuEff.H"



	//explicit U correction

	scalar tOld= runTime.elapsedCpuTime() ;

	pRelaxF=pRelaxF+0.2*(1.0-pRelaxF);
	//p.relax(pRelaxF);
	for (int corr=0; corr<nUCorr; corr++)
	{
		tmp<fvVectorMatrix> UEqn
		(
			fvm::ddt(rho, U)
		  + fvm::div(rho*phi, U)
		  - fvm::laplacian(muEff, U)
		);
		UEqn().relax(pRelaxF*((corr+0.5)/(nUCorr+1)));

		volScalarField rAU = 1.0/UEqn().A();
		//U = rAU*UEqn().H()- rAU*fvc::grad(p-beta*fvc::div(phi));
		U = (1-0.5/(nUCorr+1))*U + (0.5/(nUCorr+1))*(rAU*UEqn().H()- rAU*fvc::grad(p-beta*fvc::div(phi)));
		U.correctBoundaryConditions();
		
		U = (1-0.5/(nUCorr+1))*U + (0.5/(nUCorr+1))*(rAU*UEqn().H()- rAU*fvc::grad(p-beta*fvc::div(phi)));
		U.correctBoundaryConditions();
	}

	tmp<fvVectorMatrix> UEqn 
    (
        fvm::ddt(rho, U)
      + fvm::div(rho*phi, U)
      - fvm::laplacian(muEff, U)
    );

	UEqn.ref().relax();
		//p.relax();   
		
        //solve
        //(
            //// UEqn.ref()==fvc::reconstruct(fvc::snGrad( p-beta*fvc::div(phi) )* mesh.magSf())
            //UEqn()==-fvc::grad(p)
        //);
	 U = 1.0/UEqn().A()*(UEqn().H()-fvc::grad(p));
	 U.correctBoundaryConditions();
	 ifMonitor	Info<< "ExeTime solve U = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;


	  muEff.clear();
