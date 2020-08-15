
vector viscE(0.0,0.0,0.0);
vector gPdE(0.0,0.0,0.0);
vector rhoPhiUE(0.0,0.0,0.0);


    
{
//scalar viscE(0.0);
//scalar viscE_1(0.0);

    const fvBoundaryMesh& boundary = mesh.boundary();

    forAll(boundary, patchi)
    {//if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
			clip.boundaryField()[patchi]=0.0;
			clip.boundaryField()[patchi]==0.0;
	}

	surfaceScalarField BSurafce=fvcDelta(clip);
	
	volVectorField::Boundary& Ubf = U.boundaryField();

    forAll(boundary, patchi) if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
    {		


			scalarField muEffbp= muEff.boundaryField()[patchi];

			scalarField clipbp=    BSurafce.boundaryField()[patchi];

			vectorField bGradUp=Ubf[patchi].snGrad();//(Ubf[patchi].patchInternalField()-Ubf[patchi])*boundary[patchi].deltaCoeffs();
			bGradUp=bGradUp-boundary[patchi].Sf()*(bGradUp&boundary[patchi].Sf())/boundary[patchi].magSf()/boundary[patchi].magSf();
			viscE-=gSum(          clipbp*muEffbp*bGradUp	*boundary[patchi].magSf() );


			gPdE+=gSum(clipbp*p.boundaryField()[patchi]	*boundary[patchi].Sf() );



	rhoPhiUE +=gSum( (clipbp*phi.boundaryField()[patchi]*fvc::interpolate(rho*U)->boundaryField()[patchi]) );


	}



	surfaceVectorField gIBC= fvc::snGrad(U);
	gIBC=gIBC-mesh.Sf()*(gIBC&mesh.Sf())/mesh.magSf()/mesh.magSf();

	viscE-=gSum( (BSurafce*muEff*gIBC	*mesh.magSf())->internalField() );


	gPdE+=gSum( ( BSurafce*fvc::interpolate(p)	*mesh.Sf())->internalField() );	

	rhoPhiUE +=gSum( (BSurafce*phi*fvc::interpolate(rho*U))->internalField() );

}

