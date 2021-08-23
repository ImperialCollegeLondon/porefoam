
	volScalarField rAU = 1./UEqn.A();
	surfaceScalarField rAUf = fvc::interpolate(rAU);

	U = rAU*(UEqn.H());



	surfaceScalarField phiHbyA
	(
		  "phiU",
		  (fvc::interpolate(U) & mesh.Sf())
		+ fvc::ddtPhiCorr(rAU, rho, U, phi)
	);
	adjustPhi(phiHbyA, U, pd);

	phi = phiHbyA + ( sgPc - ghf*fvc::snGrad(rho)*mesh.magSf() )*rAUf;

	SET_REF_(pRefCell,pRefValue,Orig,alpha1);

	for(int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
	{
		fvScalarMatrix pdEqn( fvm::laplacian(rAUf, pd) == fvc::div(phi) );

		pdEqn.setReference(pRefCell, pRefValue);

		if (oCorr==nOuterCorr-1  && corr==nCorr-1  && nonOrth==nNonOrthCorr)
			pdEqn.solve(mesh.solutionDict().solver(pd.name() + "Final"));
		else
			pdEqn.solve(mesh.solutionDict().solver(pd.name()));

		if (nonOrth == nNonOrthCorr)
			phi -= pdEqn.flux();
	}

	///. U -= rAU*(fvc::grad(pd));
	U += rAU*fvc::reconstruct((phi - phiHbyA)/rAUf);
	//U = fvc::reconstruct(phi);
	U.correctBoundaryConditions();

