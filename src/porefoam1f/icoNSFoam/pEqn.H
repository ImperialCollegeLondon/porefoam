
	volScalarField rAU = 1./UEqn().A();

	U = rAU*(UEqn().H());
  
	if (pimple.nCorrPISO() <= 1)
	{
		UEqn.clear();
	}

	phi = (fvc::interpolate(U) & mesh.Sf())
		 + fvc::ddtPhiCorr(rho*rAU,U, phi);

	adjustPhi(phi, U, p);


	// Update the fixedFluxPressure BCs to ensure flux consistency
	//setSnGrad<fixedFluxPressureFvPatchScalarField> ( p.boundaryFieldRef(),
	//	( phi.boundaryField()  // - fvOptions.relative(mesh.Sf().boundaryField() & U.boundaryField())
	//	) / (mesh.magSf().boundaryField()*rAUf.boundaryField())	);

    while (pimple.correctNonOrthogonal())
	{
		fvScalarMatrix pEqn( fvm::laplacian(rAU, p) == fvc::div(phi) );

		pEqn.setReference(pRefCell, pRefValue);

		pEqn.solve
		(
			mesh.solutionDict().solver(p.select(pimple.finalInnerIter()))
		);

		if (pimple.finalNonOrthogonalIter())
			phi -= pEqn.flux();
	}


	U -= rAU*(fvc::grad(p));
	//U += rAU*fvc::reconstruct((phi - phiU)/rAUf);
	U.correctBoundaryConditions();

