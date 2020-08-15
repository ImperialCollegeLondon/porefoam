
{

	volVectorField gradP=	fvc::grad(p); 

    volVectorField::Boundary& Ubf = U.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
        if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
        { 

			 //muEff.boundaryField()[patchi]+=cBC*mag(Ubf[patchi].patchInternalField())*mag(boundary[patchi].delta());
			 muEff.boundaryField()[patchi]+=
			 0.001*cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[patchi].patchInternalField())) 
              * (boundary[patchi].delta()&boundary[patchi].delta())/
			 (mag(Ubf[patchi].patchInternalField())+(0.001*gAverage(mag(Ubf[patchi].patchInternalField()))+1.0e-15));

			 muEff.boundaryField()[patchi]*=1.0/boundary[patchi].deltaCoeffs()/
						(1.0/boundary[patchi].deltaCoeffs()+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

		}
	}
}
