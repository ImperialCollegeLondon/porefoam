
{

	volVectorField gradP=	fvc::grad(p); 

    auto& Ubf = U.boundaryFieldRef();

    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
        if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
        { 
			//- delta() is not cached, not good for memory nor efficiency, use 1/deltaCoeffs() instead

			 muEff.boundaryFieldRef()[patchi]+=
			 0.001*cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[patchi].patchInternalField()))/
			 ( boundary[patchi].deltaCoeffs()*boundary[patchi].deltaCoeffs()*( mag(Ubf[patchi].patchInternalField())+(0.001*gAverage(mag(Ubf[patchi].patchInternalField()))+1.0e-15)));

			 muEff.boundaryFieldRef()[patchi]*=1.0/boundary[patchi].deltaCoeffs()/
						(1.0/boundary[patchi].deltaCoeffs()+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

		}
	}
}
