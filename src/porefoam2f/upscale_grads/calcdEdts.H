
vector viscF(0.0,0.0,0.0);
vector viscF1(0.0,0.0,0.0);
vector gPd(0.0,0.0,0.0);
vector gPd1(0.0,0.0,0.0);    
vector GPc(0.0,0.0,0.0);
vector GPc1(0.0,0.0,0.0);    
vector rhoPhiU(0.0,0.0,0.0);
vector rhoPhiU1(0.0,0.0,0.0);
vector viscInterf_1(0.0,0.0,0.0);


    
{
//scalar viscE(0.0);
//scalar viscE_1(0.0);

    const fvBoundaryMesh& boundary = mesh.boundary();

	//volScalarField clipAlpha=alphaSh*clip;
    //forAll(boundary, patchi)
    //{//if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
			//clip.boundaryField()[patchi]=0.0;
			//clip.boundaryField()[patchi]==0.0;
			//clipAlpha.boundaryField()[patchi]=0.0;
			//clipAlpha.boundaryField()[patchi]==0.0;
	//}

	surfaceScalarField BSurafce=fvcDelta(clip);
	surfaceScalarField AlphaInternalBoundary=fvcDelta(clipAlpha);
	
	volVectorField::Boundary& Ubf = U.boundaryField();

    forAll(boundary, patchi) if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
    {		


			scalarField muEffbp= muEff.boundaryField()[patchi];

			scalarField clipbp=    BSurafce.boundaryField()[patchi];
			scalarField alpha1bp=AlphaInternalBoundary.boundaryField()[patchi];

			vectorField bGradUp=Ubf[patchi].snGrad();//(Ubf[patchi].patchInternalField()-Ubf[patchi])*boundary[patchi].deltaCoeffs();
			bGradUp=bGradUp-boundary[patchi].Sf()*(bGradUp&boundary[patchi].Sf())/boundary[patchi].magSf()/boundary[patchi].magSf();
			viscF-=gSum(          clipbp*muEffbp*bGradUp
			*boundary[patchi].magSf() );
			viscF1-=gSum( alpha1bp*muEffbp*bGradUp
			*boundary[patchi].magSf() );


			//viscE+=gSum(mag(      clipbp*muEffbp*bGradUp
			//*boundary[patchi].magSf() ));
			//viscE_1+=gSum(mag( alpha1bp*muEffbp*bGradUp
			//*boundary[patchi].magSf() ));


			gPd+=gSum(          clipbp*pd.boundaryField()[patchi]
			*boundary[patchi].Sf() );
			gPd1+=gSum( alpha1bp*pd.boundaryField()[patchi]
			*boundary[patchi].Sf() );


	rhoPhiU +=gSum( (clipbp*phi.boundaryField()[patchi]*fvc::interpolate(rho*U)->boundaryField()[patchi]) );
	rhoPhiU1+=gSum( (alpha1bp*phi.boundaryField()[patchi]*fvc::interpolate(rho*U)->boundaryField()[patchi]) );


	}



	surfaceVectorField gIBC= fvc::snGrad(U);
	gIBC=gIBC-mesh.Sf()*(gIBC&mesh.Sf())/mesh.magSf()/mesh.magSf();

	viscF-=gSum( (BSurafce*muEff*gIBC
	*mesh.magSf())->internalField() );
	viscF1-=gSum( (AlphaInternalBoundary*muEff*gIBC
	*mesh.magSf())->internalField() );

	//viscE+=gSum(mag( (BSurafce*muEff*gIBC
	//*mesh.magSf())->internalField() ));
	//viscE_1+=gSum(mag( (AlphaInternalBoundary*muEff*gIBC
	//*mesh.magSf())->internalField() ));

	viscInterf_1=-gSum( (AlphaInternalBoundary*muEff*gIBC
	*mesh.magSf())->internalField() );


	gPd+=gSum( (              BSurafce*fvc::interpolate(pd)
	*mesh.Sf())->internalField() );	
	gPd1+=gSum( (AlphaInternalBoundary*fvc::interpolate(pd)
	*mesh.Sf())->internalField() );	

	rhoPhiU +=gSum( (BSurafce*phi*fvc::interpolate(rho*U))->internalField() );

	rhoPhiU1+=gSum( (AlphaInternalBoundary*phi*fvc::interpolate(rho*U))->internalField() );



}

