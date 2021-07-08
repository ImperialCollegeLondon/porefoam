
{

	volVectorField gradP=	fvc::grad(p); 

	volVectorField::Boundary& Ubf = U.boundaryField();

	const fvBoundaryMesh& boundary = mesh.boundary();
	forAll(boundary, bi)
	{
		if (isA<fixedValueFvPatchVectorField>(Ubf[bi]))
		{ 

			//muEff.boundaryField()[bi]+=cBC*mag(Ubf[bi].patchInternalField())*mag(boundary[bi].delta());
			muEff.boundaryField()[bi]+=
			0.001*cBC* //1.0526315789474;
			(mag(gradP.boundaryField()[bi].patchInternalField())) 
				* (boundary[bi].delta()&boundary[bi].delta())/
			(mag(Ubf[bi].patchInternalField())+(0.001*gAverage(mag(Ubf[bi].patchInternalField()))+1e-15));

			muEff.boundaryField()[bi]*=1./boundary[bi].deltaCoeffs()/
				(1./boundary[bi].deltaCoeffs()+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

		}
	}
}
