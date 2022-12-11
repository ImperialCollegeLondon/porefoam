{
    word alphaScheme("div(phi,alpha)");
    word alpharScheme("div(phirb,alpha)");

    //surfaceScalarField phic = mag(phi/mesh.magSf());
    //phic = min(interface.cAlpha()*phic+0.25*max(phic), max(phic));
    //surfaceScalarField phir = phic*interface.nHatf()*(1.-mag(deltaSs)/(mag(deltaSs)+deltaN));

   scalar cAlpha( readScalar ( mesh.solutionDict().subDict("PIMPLE").lookup("cAlpha") ) );


    surfaceScalarField phic =   mag(fvc::interpolate(U));
    phic = min(cAlpha*phic, max(phic));
    surfaceScalarField phir = phic*interface.nHatf();





//#define curtail(a,b) (min (max((a),(b)),(1.-b)))

    for (int aCorr=0; aCorr<nAlphaCorr; aCorr++)
    {
        surfaceScalarField phiAlpha =
            fvc::flux
            (
                phi,
                alpha1,
                alphaScheme
            )
          + 0.75*fvc::flux
            (
                -fvc::flux(-phir, sign(scalar(1) - alpha1)*sqrt(mag(scalar(1) - alpha1)), alpharScheme),
                sign(alpha1)*sqrt(mag(alpha1)),
                alpharScheme
            )
          + 0.25*fvc::flux
            (
                -fvc::flux(-phir, (scalar(1) - alpha1), alpharScheme),
                alpha1,
                alpharScheme
            );


        //MULES::explicitSolve(alpha1, phi, phiAlpha, 1, 0);


		alphaContactAngleFvPatchScalarField::reset = true;
		alpha1.correctBoundaryConditions();
		alphaContactAngleFvPatchScalarField::reset = false;

        rhoPhi = phiAlpha*(rho1 - rho2) + phi*rho2;
		rho == alpha1*rho1 + (scalar(1) - alpha1)*rho2;

    }

    //if( runTime.timeIndex()%10== 0 )
    //{
        //Info<< "Liquid phase volume fraction = "
            //<< alpha1.weightedAverage(mesh.V()).value()
            //<< "  Min(alpha1) = " << min(alpha1).value()
            //<< "  Max(alpha1) = " << max(alpha1).value()
            //<< endl;
    //}
}
