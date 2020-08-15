
    scalar maxAlphaCo
    (
        readScalar(runTime.controlDict().lookup("maxAlphaCo"))
    );

    scalar alphaCoNum = 0.0;

    if (mesh.nInternalFaces())
    {
        //surfaceScalarField alphaf = fvc::interpolate(alpha1);
    //
        //surfaceScalarField SfUfbyDelta =
            //pos(alphaf - 0.01)*pos(0.99 - alphaf)
           //*mesh.surfaceInterpolation::deltaCoeffs();//*mag(phi);

        alphaCoNum = max(interface.deltaS()/twoPhaseProperties.muf())
            .value()*runTime.deltaT().value()*interface.sigma().value();

    }

    Info<< "Interface Courant Number "   << " max: " << alphaCoNum << endl;




// ************************************************************************* //
