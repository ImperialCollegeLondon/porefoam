

{

	gradP=	fvc::grad(pd); 

    volVectorField::Boundary& Ubf = U.boundaryField();

    const fvBoundaryMesh& boundary = mesh.boundary();
    forAll(boundary, patchi)
    {
        if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]))
        { 
//				Info<<boundary[patchi].name()<<" "<<gAverage(muEff.boundaryField()[patchi])<<  " - ";

			 //muEff.boundaryField()[patchi]+=cBC*mag(Ubf[patchi].patchInternalField())*mag(boundary[patchi].delta());
			 muEff.boundaryField()[patchi]*=1+
			 cBC* //1.0526315789474;
			 (mag(gradP.boundaryField()[patchi].patchInternalField()
			 -interface.gPc().boundaryField()[patchi].patchInternalField()
			 )) * (boundary[patchi].delta()&boundary[patchi].delta())/
			 (mag(Ubf[patchi].patchInternalField())+1e-16)/(1000.*muEff.boundaryField()[patchi]);

//				Info<<boundary[patchi].name()<<" "<<gAverage(muEff.boundaryField()[patchi])<<  " ,       ";
			 muEff.boundaryField()[patchi]*=1.0/boundary[patchi].deltaCoeffs()/
						(1.0/boundary[patchi].deltaCoeffs()+lambda);//(1.-cSlip*4.*(apS*(1-apS)));

						
		}


		//alpha1.correctBoundaryConditions();
		const volScalarField::Boundary& abf = alpha1.boundaryField();
        
        if (isA<fixedValueFvPatchVectorField>(Ubf[patchi]) && isA<alphaContactAngleFvPatchScalarField>(abf[patchi]))
        { 

			primitivePatchInterpolation pinterpolator(mesh.boundaryMesh()[patchi]);

			scalarField apS=(1./(1-2.*0.49))*(min (max(abf[patchi].patchInternalField(),0.49),(1.-0.49))-0.49) ; 
			scalarField a12S=apS;
			int intNSlip=nSlip+1e-12;
			for (int ii=0;ii<intNSlip;ii++)
			{
				apS = pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(apS) );//+0.5*alpha2sqrtp*gradAlphab[patchi];
				apS-=min(max(apS,cSSlip),(1.-cSSlip))-cSSlip;
				apS*=1./(2.0*cSSlip);
				apS =min(max(apS,0.0),1.0);
				//Info<<intNSlip<<" - "<<nSlip<<" "<<0<<endl;
			}
			a12S=apS*(1-apS);
			if ((nSlip-intNSlip)>1.e-12)
			{
			scalarField apSTmp=pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(apS) );//+0.5*alpha2sqrtp*gradAlphab[patchi];
			apSTmp-=min(max(apSTmp,cSSlip),(1.-cSSlip))-cSSlip;
			apSTmp*=1./(2.0*cSSlip);
			a12S=(nSlip-intNSlip)*apSTmp*(1-apSTmp) + (1-(nSlip-intNSlip))*a12S;

			}
				//Info<<min(apS)<<" "<<max(apS)<<" "<<average((1./boundary[patchi].deltaCoeffs())/
						//(1./boundary[patchi].deltaCoeffs()+lambdaS*4.*a12S))<<endl;
				//Info<<min(apS)<<" "<<max(apS)<<" "<<average(a12S)<<endl;
//alpha1.boundaryField()[patchi]=a12S;
//alpha1.boundaryField()[patchi]==a12S;
			muEff.boundaryField()[patchi]*=(1./boundary[patchi].deltaCoeffs())/
						(1./boundary[patchi].deltaCoeffs()+lambdaS*4.*a12S);//(1.-cSlip*4.*(apS*(1-apS)));
		}
		
		

	}
}
