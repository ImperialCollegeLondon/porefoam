
{

	volVectorField gradP=	fvc::grad(p); 


    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
		auto& Upi = U.boundaryFieldRef()[patchi];
		if (isA<fixedValueFvPatchVectorField>(Upi))
		{ 
			//- delta() is not cached, no good for memory nor efficiency, use 1/deltaCoeffs() instead
			auto& delCoefs = boundary[patchi].deltaCoeffs();

			 muEff.boundaryFieldRef()[patchi]+=
			 min(0.001*cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[patchi].patchInternalField()))/
			 ( delCoefs*delCoefs*( mag(Upi.patchInternalField())+(0.001*gAverage(mag(Upi.patchInternalField()))+1.0e-9))), 1.0*mu.value());

			 muEff.boundaryFieldRef()[patchi]*=1.0/delCoefs/(1.0/delCoefs+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

		}
	}
}
