
{

	volVectorField gradP=	fvc::grad(p); 

    auto& Ubf = U.boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, bi)
    {
        if (isA<fixedValueFvPatchVectorField>(Ubf[bi]))
        { 
			//- delta() is not cached, not good for memory nor efficiency, use 1/deltaCoeffs() instead

			 muEff.boundaryFieldRef()[bi]+=
			 0.001*cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[bi].patchInternalField()))/
			 ( boundary[bi].deltaCoeffs()*boundary[bi].deltaCoeffs()*( mag(Ubf[bi].patchInternalField())+(0.001*gAverage(mag(Ubf[bi].patchInternalField()))+1e-15)));

			 muEff.boundaryFieldRef()[bi]*=1./boundary[bi].deltaCoeffs()/
						(1./boundary[bi].deltaCoeffs()+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

		}
	}
}
