    muEff=twoPhaseProperties.muf();


#include  "correctmuEff.H"


	fvVectorMatrix UEqn
    (
        fvm::ddt(rho, U)
      + fvm::div(rhoPhi, U)
      - fvm::laplacian(muEff, U)
      - (fvc::grad(U) & fvc::grad( surfaceScalarField("muEffOrig",twoPhaseProperties.muf()) ))
    //- fvc::div(muEff*(fvc::interpolate(dev(fvc::grad(U))) & mesh.Sf()))
    );
	if (oCorr == nOuterCorr-1 )
	{
		UEqn.relax(1);
	}
	else
	{
		UEqn.relax();
		pd.relax();   
	}

	if(filterPd)
	{
		volScalarField rAU = 1./UEqn.A();

		//scalar AverageMagPdBeforeFilter=average(mag(pd)).value();
		pFilter= pBoundBypSmoothFactor*(fvc::average(fvc::interpolate(pd)));
		pFilter.correctBoundaryConditions();
		pd = min(max(pd, -mag(pFilter)),mag(pFilter));

		pd.correctBoundaryConditions();
		//scalar AverageMagpAfterFilter=average(mag(pd)).value();

		//Info<< "AverageMagPdBeforeFilter: "<<AverageMagPdBeforeFilter<<"   AverageMagpAfterFilter:" <<AverageMagpAfterFilter<<endl;
	}

    if (0)//pimple.momentumPredictor())
    {
		scalar tOld= runTime.elapsedCpuTime() ;

        solve
        (
            UEqn
         ==
			fvc::reconstruct
						(
							sgPc-
							(
								ghf*fvc::snGrad(rho)
							  + fvc::snGrad(pd)
							) * mesh.magSf()
						)
        );
		Info<< "ExeTime solve U = " << runTime.elapsedCpuTime()-tOld << " s"	<< endl;
		

		
    }
    else
	{

		volVectorField gP = fvc::reconstruct
							(
								(
									ghf*fvc::snGrad(rho)
								    + fvc::snGrad( pd-beta*fvc::div(phi) )
								    
								) * mesh.magSf()
								- sgPc
							);


		//U =0.5*U + 0.5/UEqn.A()*( UEqn.H() -gP );
		U =1./UEqn.A()*( UEqn.H() - gP );

	}
	U.correctBoundaryConditions();

















	if(filterU)
	{
		//scalar AverageMagUBeforeFilter=average(mag(U)).value();
		UFilter= UBoundByUSmoothFactor*(fvc::average(fvc::interpolate(U,"limitedScheme")*Internalfaces1)/fvc::average(Internalfaces1));
		UFilter.correctBoundaryConditions();
		forAll(U,i)
		{
			U[i][0] = min (max(U[i][0], -mag(UFilter[i][0])),mag(UFilter[i][0]));
			U[i][1] = min (max(U[i][1], -mag(UFilter[i][1])),mag(UFilter[i][1]));
			U[i][2] = min (max(U[i][2], -mag(UFilter[i][2])),mag(UFilter[i][2]));
		}

		UFilter*= dimensionedScalar(1./UBoundByUSmoothFactor);
		UFilter.correctBoundaryConditions();

		 U.correctBoundaryConditions();
		//scalar AverageMagUAfterFilter=average(mag(U)).value();

		//Info<< "AverageMagUBeforeFilter: "<<AverageMagUBeforeFilter<<"   AverageMagUAfterFilter:" <<AverageMagUAfterFilter<<endl;
	}
