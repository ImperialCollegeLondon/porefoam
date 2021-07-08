
	surfaceScalarField muEff(IOobject("muEff", runTime.timeName(),  mesh),  mesh,  rho*nu);
	#include  "correctmuEff.H"



	//Aitken U correction
	if(AitkenEvery){
		scalar UAk = average(mag(U)).value();
		ifMonitor(AitkenEvery) {
			scalar delU = UAk-UAkOOld;  scalar del2U = delU-(UAkOOld-UAkOOOOld);
			scalar dUAkF= mag(del2U)>1e-18 ? -delU*delU/del2U : 0.9*delU;  
			dUAkF = min(max(-0.1,dUAkF/max(1e-16,UAk)),0.1);
			U *= 1+dUAkF;
			Info<<"\ndUAitkenF:"<<dUAkF<<nl<<endl; 
		}
		UAkOOOOld=UAkOOOld; UAkOOOld=UAkOOld;  UAkOOld=UAkOld;  UAkOld=UAk;
	}


	//explicit U correction

	scalar tOld= runTime.elapsedCpuTime() ;

	pRelaxF=pRelaxF+0.2*(1.-pRelaxF);
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

		volScalarField rAU = 1./UEqn().A();
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
	 U = 1./UEqn().A()*(UEqn().H()-fvc::grad(p));
	 U.correctBoundaryConditions();
	 ifMonitor(10)	Info<< "ExeTime solve U = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;


	  muEff.clear();
