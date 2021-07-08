
{

	volVectorField gradP=	fvc::grad(p); 


    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, bi)
    {
		auto& Upi = U.boundaryFieldRef()[bi];
		if (isA<fixedValueFvPatchVectorField>(Upi))
		{ 
			//- delta() is not cached, no good for memory nor efficiency, use 1/deltaCoeffs() instead
			auto& delCoefs = boundary[bi].deltaCoeffs();

			 muEff.boundaryFieldRef()[bi]+=
			 min(0.001*cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[bi].patchInternalField()))/
			 ( delCoefs*delCoefs*( mag(Upi.patchInternalField())+(0.001*gAverage(mag(Upi.patchInternalField()))+1e-9))), 1.*mu.value());

			 muEff.boundaryFieldRef()[bi]*=1./delCoefs/(1./delCoefs+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

		}
	}
}
