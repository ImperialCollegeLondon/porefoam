

{

	gradP=	fvc::grad(pd); 

    auto& Ubf = U.boundaryField();

    const fvBoundaryMesh& bndrys = mesh.boundary();
    forAll(bndrys, bi)
    {
        if (isA<fixedValueFvPatchVectorField>(Ubf[bi]))
        { 
//				Info<<bndrys[bi].name()<<" "<<gAverage(muEff.boundaryField()[bi])<<  " - ";

			 //muEff.boundaryField()[bi]+=cBC*mag(Ubf[bi].patchInternalField())*mag(bndrys[bi].delta());
			 muEff.boundaryFieldRef()[bi]*=1+
			 cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[bi].patchInternalField()
			 -interface.gPc().boundaryField()[bi].patchInternalField()
			 )) * (bndrys[bi].delta()&bndrys[bi].delta())/
			 (mag(Ubf[bi].patchInternalField())+1e-16)/(1000.*muEff.boundaryField()[bi]);

//				Info<<bndrys[bi].name()<<" "<<gAverage(muEff.boundaryField()[bi])<<  " ,       ";
			 muEff.boundaryFieldRef()[bi]*=1./bndrys[bi].deltaCoeffs()/
						(1./bndrys[bi].deltaCoeffs()+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

						
		}


		//alpha1.correctBoundaryConditions();
		const volScalarField::Boundary& abf = alpha1.boundaryField();
        
        if (isA<fixedValueFvPatchVectorField>(Ubf[bi]) && isA<alphaContactAngleFvPatchScalarField>(abf[bi]))
        { 

			primitivePatchInterpolation pinterpolator(mesh.boundaryMesh()[bi]);

			scalarField apS=(1./(1-2.*0.49))*(min (max(abf[bi].patchInternalField(),0.49),(1.-0.49))-0.49) ; 
			scalarField a12S=apS;
			int intNSlip=nSlip+1e-12;
			for (int ii=0;ii<intNSlip;ii++)
			{
				apS = pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(apS) );//+0.5*alpha2sqrtp*gradAlphab[bi];
				apS-=min(max(apS,cSSlip),(1.-cSSlip))-cSSlip;
				apS*=1./(2.*cSSlip);
				apS =min(max(apS,0.),1.);
				//Info<<intNSlip<<" - "<<nSlip<<" "<<0<<endl;
			}
			a12S=apS*(1-apS);
			if ((nSlip-intNSlip)>1e-12)
			{
			scalarField apSTmp=pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(apS) );//+0.5*alpha2sqrtp*gradAlphab[bi];
			apSTmp-=min(max(apSTmp,cSSlip),(1.-cSSlip))-cSSlip;
			apSTmp*=1./(2.*cSSlip);
			a12S=(nSlip-intNSlip)*apSTmp*(1-apSTmp) + (1-(nSlip-intNSlip))*a12S;

			}
				//Info<<min(apS)<<" "<<max(apS)<<" "<<average((1./bndrys[bi].deltaCoeffs())/
						//(1./bndrys[bi].deltaCoeffs()+lambdaS*4.*a12S))<<endl;
				//Info<<min(apS)<<" "<<max(apS)<<" "<<average(a12S)<<endl;
//alpha1.boundaryField()[bi]=a12S;
//alpha1.boundaryField()[bi]==a12S;
			muEff.boundaryFieldRef()[bi]*=(1./bndrys[bi].deltaCoeffs())/
						(1./bndrys[bi].deltaCoeffs()+lambdaS*4.*a12S);//(1.-cSlip*4.*(apS*(1-apS)));
		}
		
		

	}
}
