surfaceVectorField interfaceSf
  (
      IOobject
      (
          "interfaceSf",
          alpha1_.time().timeName(),
          alpha1_.mesh(),
          IOobject::NO_READ,
          IOobject::AUTO_WRITE
      ),
      alpha1_.mesh(),
      dimensionedVector("interfaceSf", dimArea, vector::zero)
  );



	// calculate interface normal vectors, interfaceSf
    forAll(magDelS,faceI)
    {
        if (magDelS[faceI] > CONTRAST_FACTOR)
        {
            // defining the right hand rule coeficient, such that interfaceSf_ always points from 0(Fluid1) to 1(Fluid2)
            scalar rightHandRuleCoef =deltaS_[faceI]/magDelS[faceI];

            // return the points of faceI
            const face& facePoints = faces[faceI];

            pointField 	xInterface(facePoints.size());

            forAll(facePoints,pointI)
            {
                label ip1(facePoints[pointI]);
                xInterface[pointI] = points[ip1] + distPointInterface_[ip1];
            }

            Foam::label nPoints = xInterface.size();

            if (nPoints == 3)
            {
                interfaceSf[faceI] = rightHandRuleCoef*triPointRef
                (
                    xInterface[0],
                    xInterface[1],
                    xInterface[2]
                ).normal();
            }
            else
            {
                point centrePoint = point::zero;

                for(Foam::label pointI=0;pointI<nPoints;++pointI)
                {
                    centrePoint += xInterface[pointI];
                }

                centrePoint /= nPoints;

                for (Foam::label pointI=0; pointI<nPoints; ++pointI)
                {
                    const point& currentPoint = xInterface[pointI];
                    const point& nextPoint = xInterface[(pointI + 1) % nPoints];

                    interfaceSf[faceI] += rightHandRuleCoef*triPointRef
                    (
                        currentPoint,
                        nextPoint,
                        centrePoint
                    ).normal();
                }
            }
        }
    }

		forAll(patches, patchI)
		{
			const scalarField & pMagDelS = magDelS.boundaryField()[patchI];
			const scalarField & pSnGradAlphaSharp = deltaS_.boundaryField()[patchI];
			vectorField & pInterfaceSf = interfaceSf.boundaryField()[patchI];
			if (patches[patchI].coupled())
			{
				forAll(pMagDelS, pfI)
				{
				  if (pMagDelS[pfI] > CONTRAST_FACTOR)
				  {
						// defining the right hand rule coeficient, such that interfaceSf_ always points from 0(Fluid1) to 1(Fluid2)
						scalar rightHandRuleCoef = pSnGradAlphaSharp[pfI]/pMagDelS[pfI];

						// return the points of faceI
						const face& facePoints = faces[pfI+patches[patchI].patch().start()];

						pointField 	xInterface(facePoints.size());

						forAll(facePoints,pointI)
						{
							 label ip1(facePoints[pointI]);
							 xInterface[pointI] = points[ip1] + distPointInterface_[ip1];
						}

						Foam::label nPoints = xInterface.size();

						if (nPoints == 3)
						{
							 pInterfaceSf[pfI] = rightHandRuleCoef*triPointRef
							 (
								  xInterface[0],
								  xInterface[1],
								  xInterface[2]
							 ).normal();
						}
						else
						{
							 point centrePoint = point::zero;

							 for(Foam::label pointI=0;pointI<nPoints;++pointI)
							 {
								  centrePoint += xInterface[pointI];
							 }

							 centrePoint /= nPoints;

							 for (Foam::label pointI=0; pointI<nPoints; ++pointI)
							 {
								  const point& currentPoint = xInterface[pointI];
								  const point& nextPoint = xInterface[(pointI + 1) % nPoints];

								  pInterfaceSf[pfI] += rightHandRuleCoef*triPointRef
								  (
										currentPoint,
										nextPoint,
										centrePoint
								  ).normal();
							 }
						}
					}
				}

			}
			else		 {			 }
		}

	//interfaceSf=(interfaceSf&curvatureNSf)*curvatureNSf;
	//interfaceSf=mesh.magSf()/(curvatureNSf&mesh.Sf()+1.0e-36)*mesh.Sf();


