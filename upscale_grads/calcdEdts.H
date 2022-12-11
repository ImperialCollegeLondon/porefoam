
vector viscF(0.,0.,0.);
vector viscF1(0.,0.,0.);
vector gPd(0.,0.,0.);
vector gPd1(0.,0.,0.);
vector GPc(0.,0.,0.);
vector GPc1(0.,0.,0.);
vector rhoPhiU(0.,0.,0.);
vector rhoPhiU1(0.,0.,0.);
vector viscInterf_1(0.,0.,0.);



{
//scalar viscE(0.);
//scalar viscE_1(0.);

    const fvBoundaryMesh& boundary = mesh.boundary();

	//volScalarField clipAlpha=alphaSh*clip;
    //forAll(boundary, bi)
    //{//if (isA<fixedValueFvPatchVectorField>(Ubf[bi]))
			//clip.boundaryFieldRef()[bi]=0.;
			//clip.boundaryFieldRef()[bi]==0.;
			//clipAlpha.boundaryFieldRef()[bi]=0.;
			//clipAlpha.boundaryFieldRef()[bi]==0.;
	//}

	surfaceScalarField BSurafce=fvcDelta(clip);
	surfaceScalarField AlphaInsidBoundary=fvcDelta(clipAlpha);

	auto& Ubf = U.boundaryField();

    forAll(boundary, bi) if (isA<fixedValueFvPatchVectorField>(Ubf[bi]))
    {


			scalarField muEffbp= muEff.boundaryField()[bi];

			scalarField clipbp=    BSurafce.boundaryField()[bi];
			scalarField alpha1bp=AlphaInsidBoundary.boundaryField()[bi];

			vectorField bGradUp=Ubf[bi].snGrad();//(Ubf[bi].patchInternalField()-Ubf[bi])*boundary[bi].deltaCoeffs();
			bGradUp=bGradUp-boundary[bi].Sf()*(bGradUp&boundary[bi].Sf())/boundary[bi].magSf()/boundary[bi].magSf();
			viscF-=gSum(          clipbp*muEffbp*bGradUp
			*boundary[bi].magSf() );
			viscF1-=gSum( alpha1bp*muEffbp*bGradUp
			*boundary[bi].magSf() );


			//viscE+=gSum(mag(      clipbp*muEffbp*bGradUp
			//*boundary[bi].magSf() ));
			//viscE_1+=gSum(mag( alpha1bp*muEffbp*bGradUp
			//*boundary[bi].magSf() ));


			gPd+=gSum( clipbp*pd.boundaryField()[bi]
			*boundary[bi].Sf() );
			gPd1+=gSum( alpha1bp*pd.boundaryField()[bi]
			*boundary[bi].Sf() );


	rhoPhiU +=gSum( (clipbp*phi.boundaryField()[bi]*fvc::interpolate(rho*U)->boundaryField()[bi]) );
	rhoPhiU1+=gSum( (alpha1bp*phi.boundaryField()[bi]*fvc::interpolate(rho*U)->boundaryField()[bi]) );


	}



	surfaceVectorField gIBC= fvc::snGrad(U);
	gIBC=gIBC-mesh.Sf()*(gIBC&mesh.Sf())/mesh.magSf()/mesh.magSf();

	viscF-=gSum( (BSurafce*muEff*gIBC
	*mesh.magSf())->internalField() );
	viscF1-=gSum( (AlphaInsidBoundary*muEff*gIBC
	*mesh.magSf())->internalField() );

	//viscE+=gSum(mag( (BSurafce*muEff*gIBC
	//*mesh.magSf())->internalField() ));
	//viscE_1+=gSum(mag( (AlphaInsidBoundary*muEff*gIBC
	//*mesh.magSf())->internalField() ));

	viscInterf_1=-gSum( (AlphaInsidBoundary*muEff*gIBC
	*mesh.magSf())->internalField() );

	gPd+=gSum( (BSurafce*fvc::interpolate(pd)*mesh.Sf())->internalField() );
	gPd1+=gSum( (AlphaInsidBoundary*fvc::interpolate(pd)*mesh.Sf())->internalField() );

	rhoPhiU +=gSum( (BSurafce*phi*fvc::interpolate(rho*U))->internalField() );

	rhoPhiU1+=gSum( (AlphaInsidBoundary*phi*fvc::interpolate(rho*U))->internalField() );



}
