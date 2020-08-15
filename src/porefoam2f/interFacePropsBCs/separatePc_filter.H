
#define SET_REF_(rfc, rfv, sfx, _alpha1_) \
	for(rfc = rfc##sfx; rfc<_alpha1_.size()-5; ++rfc) if ((_alpha1_[rfc]<0.01) && (magGradAlpha[rfc]<100.0*deltaN_.value()) )  break;   \
	rfv=0.995*pc_[rfc]+0.005*rfv##sfx 


{ ///. pc computation and filtering

	fcCorrectTangent_=   readScalar( pimple.lookup("fcCorrectTangent") );
	fcCorrectTangentRelax_=   readScalar( pimple.lookup("fcCorrectTangentRelax") );
	fcdFilter_=   readScalar( pimple.lookup("fcdFilter") );
	nPcNonOrthCorr_=   readLabel( pimple.lookup("nPcNonOrthCorr") );

    pcRefCellOrig_= readLabel( pimple.lookup("pcRefCell") );
    pcRefValueOrig_= readScalar( pimple.lookup("pcRefValue") );

	volVectorField gradAlpha = fvc::grad(alpha1f_); 
	volScalarField magGradAlpha = (mag(gradAlpha) + 0.0001*deltaN_); 
	volVectorField nS = gradAlpha/magGradAlpha;
	surfaceVectorField nHatfv = fvc::interpolate(gradAlpha)/fvc::interpolate(magGradAlpha); nHatfv/=mag(nHatfv)+1.0e-15;
	//correctContactAngle(nHatfv.boundaryField(),gradAlpha.boundaryField(),nS.boundaryField(),alpha1S_.boundaryField(), a1a2.boundaryField());
	nHatf_ = nHatfv & Sf;

	surfaceScalarField dInterf=(mag(deltaS_)/(mag(deltaS_)+1.0e-4*deltaN_));
	//surfaceScalarField dInterf=/*0.001*(mag(deltaS_)/(mag(deltaS_)+deltaN_))+0.999**/(mag(snGradaS)/(mag(snGradaS)+deltaN_));
	sgPcErr_*=fcCorrectTangentRelax_*dInterf;
	sgPc_=(stf -sgPcErr_);

	SET_REF_(pcRefCell,pcRefValue,Orig_,alpha1_);

    for(int nonOrth=0; nonOrth<=nPcNonOrthCorr_; nonOrth++)
    {

        fvScalarMatrix pcEqn
        (
            fvm::laplacian( pc_) == fvc::div(sgPc_)
        );

        pcEqn.setReference(pcRefCell, pcRefValue);

        pcEqn.solve();


        if (nonOrth == nPcNonOrthCorr_)
        { 
            sgPc_-=pcEqn.flux();
        }
		pc_.correctBoundaryConditions();

    }


	if(!std::isfinite(average(pc_).value()))
	{// code used to crash
		Info<<"!std::isfinite(average(pc_)), sgPc_:"<<min(sgPc_)<<max(sgPc_)<<endl;
		Info<<"!std::isfinite(average(pc_)), sgPcErr_:"<<min(sgPcErr_)<<max(sgPcErr_)<<endl;
		Info<<"!std::isfinite(average(pc_)), stf:"<<min(stf)<<max(stf)<<endl;
		Info<<"pcRefValue:"<<pcRefValue<<endl;
		Info<<"pcRefCell:"<<pcRefCell<<endl;
		sleep(1);
		std::terminate();
	}
	//ensure(std::isfinite(average(pc_)),"",2);
    
	surfaceScalarField stfTypicalLow=mag((sigma_*0.1)*deltaS_*deltaS_*mesh.magSf()); // Note: this filtering is only applied to the interface region, but not away from interface
	surfaceScalarField stfThreshold=((0.5*fcdFilter_)*( mag(stf)+ stfTypicalLow)); // Note: this filtering is only applied to the interface region, but not away from interface
	//sgPcErrn_ = (min (max(sgPc_, -stfThreshold),stfThreshold));
	//sgPc_-=sgPcErrn_;  
    
        
    gPc_= fvc::reconstruct((sgPc_));
    gPc_.correctBoundaryConditions();



    volVectorField vgPcllInterface=(gPc_)-((gPc_)&(nS))*(nS);    //surfaceVectorField gPcf=fvc::interpolate(gPc_,"smoothSchemegPc");

	surfaceVectorField gPcllInterface=fvc::interpolate(vgPcllInterface,"smoothSchemegPc");    //surfaceVectorField gPcf=fvc::interpolate(gPc_,"smoothSchemegPc");
	gPcllInterface=(gPcllInterface)-((gPcllInterface)&(nHatfv))*(nHatfv);    //surfaceVectorField gPcf=fvc::interpolate(gPc_,"smoothSchemegPc");
	gPcllInterface*=(1.0-mag(nHatf_)/mesh.magSf());
	sgPcErr_=dInterf*(min(2.0*(1.4-mag(nHatf_)/mesh.magSf()), 1.0)*(sgPcErr_+fcCorrectTangent_*(gPcllInterface& mesh.Sf())) );//  *( 1.0-(mag(nHatf_)/mesh.magSf()) )   ); /*aSolid_cPcCorrRelax_* */ 


	long long sgPcErr_checkErrs=0;
	forAll(sgPcErr_,i)
	{
		if (sgPc_[i]*sgPcErr_[i] < -1.0e-6*stfTypicalLow[i])  sgPcErr_[i]  *= 0.99;
		if( mag(sgPc_[i])<0.5*mag(sgPcErr_[i]) ) {   sgPcErr_[i] = 0.5*sgPc_[i]; } //++sgPcErr_checkErrs;    Info<<"bEC:"<<sgPcErr_[i]<<" "<<sgPc_[i]<<endl; 
	}
	forAll(boundary, patchI)
	{
		if (!boundary[patchI].coupled())
		{
			sgPcErr_.boundaryField()[patchI] = 0.0;
		}
		else
		{
			Field<scalar> & psgPcErr_ = sgPcErr_.boundaryField()[patchI];
			Field<scalar> & psgPc_ = sgPc_.boundaryField()[patchI];
			forAll(psgPcErr_,i)
			{
				if(psgPc_[i]*psgPcErr_[i] < -1.0e-6*stfTypicalLow[i])  psgPcErr_[i] *= 0.99;
				if( mag(psgPc_[i])<0.5*mag(psgPcErr_[i]) ) {  psgPcErr_[i] = 0.5*psgPc_[i];} // ++sgPcErr_checkErrs;   Info<<"bEC:"<<psgPcErr_[i]<<" "<<psgPc_[i]<<endl;
			}
		}
	}


	if(sgPcErr_checkErrs)  Info<<"\n\n sgPcErr_checkErrs: "<<sgPcErr_checkErrs<<endl<<endl;

	//stfThreshold=0.5*(stfThreshold+(1.*fcdFilter_)*mag(stf));
	stfThreshold=((1.0*fcdFilter_)*( mag(stf)+ mag((sigma_*0.1)*deltaS_*deltaS_*mesh.magSf()) ));
	sgPcErrn_ = (min (max(sgPc_, -stfThreshold),stfThreshold) );
	sgPc_-=sgPcErrn_;  

	gPc_= fvc::reconstruct((sgPc_));
	gPc_.correctBoundaryConditions();


	volVectorField gPcS=fvc::average(fvc::interpolate(gPc_)*Internalfaces1_)/fvc::average(Internalfaces1_);
	stfThreshold = 2.0*mag((fvc::interpolate(gPcS)) & mesh.Sf());
	sgPc_ = (min (max(sgPc_, -stfThreshold),stfThreshold));

	gPc_= fvc::reconstruct((sgPc_));
	gPc_.correctBoundaryConditions();



//this part can be uncommented but I haven't seen any noticable difference

//   gPc_= fvc::reconstruct((sgPc_));
        //{
            //sgPc_=(stf -sgPcErr_-2.*sgPcErrn_);
            //for(int nonOrth=0; nonOrth<=nPcNonOrthCorr_; nonOrth++)
            //{
                //fvScalarMatrix pcEqn
                //(
                    //fvm::laplacian( pc_) == fvc::div(sgPc_)
                //);
                //pcEqn.setReference(pcRefCell_, pcRefValue_);
                //pcEqn.solve();
                //if (nonOrth == nPcNonOrthCorr_)
                //{ 
                    //sgPc_-=pcEqn.flux();
                //}
            //}
        //}

//   gPc_= fvc::reconstruct((sgPc_));

	Info<< "~~~Pc computation time = " << alpha1_.time( ).elapsedCpuTime()-stime << " s"	<< endl;

}

