
//vector FF(0., 0., 0.);
{
	Info<< "Reading field p\n" << endl;

	//volScalarField& V=p;//("V",p);

	Info<< "Reading transportProperties\n" << endl;



	Info<< "Reading diffusivity K\n" << endl;



	//dimensionedScalar fFact = transportProperties.lookupOrDefault("fF",dimensionedScalar("fF", dimless, 1.));

	volScalarField kEf
	(	IOobject
		(	"kEf",	"0",	mesh,
			IOobject::READ_IF_PRESENT,	IOobject::NO_WRITE
		),
		1./porosity/porosity,"zeroGradient"
	);
	kEf.correctBoundaryConditions();
	p.correctBoundaryConditions();
	Info<<"avg kEf: "<<average(kEf).value()<<",  ["<<min(kEf).value()<<" - "<<max(kEf).value()<<"]"<<endl;
	Info<<"avg p init: "<<average(p).value()<<",  ["<<min(p).value()<<" - "<<max(p).value()<<"]"<<endl;

	pimpleControl pimple(mesh);

	// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

	Info<< "\ncalculating formation factor \n" << endl;

	for(int iter = 0; iter<4; ++iter)
	{
 
		//runTime++;

		while (pimple.loop())
		{

			Info<< "Time = " << runTime.timeName() << nl << endl;

			while (pimple.correctNonOrthogonal())
			{
				solve
				(
					fvm::laplacian(kEf, p)
				);
			}

			Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
				<< "  ClockTime = " << runTime.elapsedClockTime() << " s"
				<< nl << endl;
		}
	}
	
	
	phi=fvc::interpolate(kEf*dimensionedScalar("one", dimLength*dimLength, 1.)/(mu))*fvc::snGrad(p)*mesh.magSf();
	Info<<"avg phi psi: "<<average(phi).value()<<",  ["<<min(phi).value()<<" - "<<max(phi).value()<<"]"<<endl;

	Info<< ".\n" << endl;

	{
		Info<<"Writting (Electrical) potential, psi"<<endl;
		volScalarField psi
		(	IOobject
			(	"psi",	"0",	mesh,
				IOobject::NO_READ,	IOobject::AUTO_WRITE
			),
			p
		);	psi.write();
	}


	{
		vector dV;
		label iLeft = mesh.boundaryMesh().findPatchID(LeftPs[iDir]);
		if (iLeft < 0)	 { Info	<< "Unable to find  patch " << LeftPs[iDir] << nl	<< endl;	}
		label iRight = mesh.boundaryMesh().findPatchID(RightPs[iDir]);
		if (iRight < 0) { Info	<< "Unable to find  patch " << RightPs[iDir] << nl	<< endl;	}
	
		Info<<"kEf_ iRight: "<<gSum(kEf.boundaryField()[iRight])<<"  "<<endl;
		Info<<"kEf_ iLeft: "<<gSum(kEf.boundaryField()[iLeft])<<"  "<<endl;
		Info<<"psi_ iRight: "<<gSum(p.boundaryField()[iRight])<<"  "<<endl;
		Info<<"psi_ iLeft: "<<gSum(p.boundaryField()[iLeft])<<"  "<<endl;

		if( iLeft>=0 && iRight>=0)
		{

			scalar fluxIn=gSum(phi.boundaryField()[iLeft]+1e-29);
			scalar fluxOut=gSum(phi.boundaryField()[iRight]+1e-29);

			scalar PLeft=gSum(p.boundaryField()[iLeft]*(phi.boundaryField()[iLeft]))/(fluxIn+1e-129);
			scalar PRight=gSum(p.boundaryField()[iRight]*(phi.boundaryField()[iRight]))/(fluxOut+1e-129);
			dV[iDir]=mag(PLeft-PRight);

			fluxIn=gSum(phi.boundaryField()[iLeft]);
		Info<<"fluxIn: "<<fluxIn<<"  @dir"<<iDir<<endl;
		Info<<"dV[iDir]: "<<dV[iDir]<<" = "<<PLeft<<" - "<<PRight<<iDir<<endl;
		Info<<"A/L: "<<A[iDir]/L[iDir]<<endl;
		Info<<"mu.value(): "<<mu.value()<<endl;

			FF[iDir]=dV[iDir]*A[iDir]/L[iDir]/(mag(fluxIn)+1e-126)/(mu.value());
		Info<<"FF[iDir]: "<<FF[iDir]<<endl;
		}
	}
}
