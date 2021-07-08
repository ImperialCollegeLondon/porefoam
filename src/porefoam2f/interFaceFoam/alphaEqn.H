{
	#define curtailBADOFSET(a,b) (min (max((a),(b)),(1.-(b))))

	word alphaScheme("div(phi,alpha)");
	word alpharScheme("div(phirb,alpha)");

	scalar cAlpha( readScalar ( mesh.solutionDict().subDict("PIMPLE").lookup("cAlpha") ) );
	const polyBoundaryMesh& boundary = mesh.boundaryMesh();


	volScalarField alphaSm
	(
		//IOobject( "alphaSm",  runTime.timeName(),  mesh ),
		(1./0.94)*(curtailBADOFSET(alpha1,0.03)-0.03)
		//,"zeroGradient"
	);
	alphaSm=0.7*alphaSm+0.3*fvc::average(fvc::interpolate(alpha1));
			forAll(boundary, bi)if(!boundary[bi].coupled()) { alphaSm.boundaryField()[bi]==alphaSm.boundaryField()[bi].patchInternalField();	}

	surfaceVectorField gradAlphaf = fvc::interpolate(fvc::grad(alphaSm)); 
									//forAll(boundary, bi)if(!boundary[bi].coupled()) { gradAlphaf.boundaryField()[bi]*=0.;	}
	gradAlphaf = gradAlphaf+0.1*fvc::interpolate(fvc::average(gradAlphaf)); 
	gradAlphaf = gradAlphaf+0.1*fvc::interpolate(fvc::average(gradAlphaf)); 
	surfaceScalarField nHatf = (gradAlphaf& mesh.Sf())/(mag(gradAlphaf)+interface.deltaN()); // nHatfv/=mag(nHatfv)+1e-12;
	

	surfaceScalarField phic = mag(phi/mesh.magSf());
	phic = phic+fvc::interpolate(fvc::average(phic*interface.Internalfaces1())/interface.AvgInternFaces1(),"localMax"); 
	phic = phic+fvc::interpolate(fvc::average(phic*interface.Internalfaces1())/interface.AvgInternFaces1(),"localMax"); 
	//phic = 0.5*phic+0.5*fvc::interpolate(fvc::average(phic),"localMax"); 
	//phic =  cAlpha*phic*(1.+interface.BInternalfs()) ;
	//phic = min(cAlpha*phic*(1.+interface.BInternalfs()), 1.01*max(phic));
				forAll(boundary, bi)if(!boundary[bi].coupled()){  phic.boundaryField()[bi]*=0.;	}
	surfaceScalarField isS = mag(fvc::snGrad(curtailBADOFSET(alphaSm,0.48)));
	surfaceScalarField dNotInterf=1.-0.9*isS/(isS+interface.deltaN());
	surfaceScalarField phir = min(cAlpha*dNotInterf*phic, 1.*max(phic))*nHatf;


	volScalarField alphar =  2.*curtailBADOFSET(alpha1,0.)-alpha1;

	///. to test the effect of 1.00 ;  ?-> 0.80
	//scalar cPcpEps=0.2;
	//volScalarField 
	//alphaSpread=(1./(cPcpEps*2.))*(alpha1-curtailBADOFSET(alpha1,cPcpEps))+(0.5*1.00);//+(1.00-0.20)*alpha1;
	//alphaSpread.correctBoundaryConditions();

	for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
	{
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + fvc::flux
            (
                -fvc::flux(-phir, sqrt(mag(scalar(1) - alphar)), alpharScheme),
                sqrt(mag(alphar)),
                alpharScheme
            );

		#ifndef  POST_PROCCESS
        MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);


		alphaContactAngleFvPatchScalarField::reset = true;
		alpha1.correctBoundaryConditions();
		alphaContactAngleFvPatchScalarField::reset = false;
        #endif //POST_PROCCESS

        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
		rho = alpha1*rho1 + (scalar(1) - alpha1)*rho2;//+0.*(alpha1*(1.-alpha1))*(rho1+rho2);

	}

}












