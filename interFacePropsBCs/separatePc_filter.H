

// TODO clean up & convert to a function, sync_with:findRefCell.C
#define SET_REF_(rfc, rfv, sfx, _alpha1SAtanh_) 			if (Pstream::master()) { \
	constexpr scalar aSp01=std::atanh(1.8*0.01-0.9); /* SYN123433*/ \
	for(rfc = rfc##sfx; rfc<_alpha1SAtanh_.size()-5; ++rfc) \
		if ((_alpha1SAtanh_[rfc]<aSp01 || _alpha1SAtanh_[rfc]>-aSp01) && (magGradAlpha[rfc]<1000.*deltaN_.value()) && IsRefCandid_[rfc])  break; \
	if (rfc==_alpha1SAtanh_.size()) Info<<"Error cannot find suitable refCell for pc solver"<<endl;\
	rfv=0.9999*pc_[rfc]+0.0001*rfv##sfx; \
	} else rfc=-1;


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
	surfaceVectorField nHatfv = fvc::interpolate(gradAlpha)/fvc::interpolate(magGradAlpha); nHatfv/=mag(nHatfv)+1e-15;
	//correctContactAngle(nHatfv.boundaryFieldRef(),gradAlpha.boundaryFieldRef(),nS.boundaryFieldRef(),alpha1S_.boundaryFieldRef(), a1a2.boundaryFieldRef());
	nHatf_ = nHatfv & Sf;

	surfaceScalarField dInterf=(mag(deltaS_)/(mag(deltaS_)+1e-4*deltaN_));
	//surfaceScalarField dInterf=/*0.001*(mag(deltaS_)/(mag(deltaS_)+deltaN_))+0.999**/(mag(snGradaS)/(mag(snGradaS)+deltaN_));
	sgPcErr_*=fcCorrectTangentRelax_*dInterf;
	sgPc_=(stf -sgPcErr_);


	pc_ = 0.9*pc_+0.1*fvc::average(fvc::interpolate(pc_));// hoping to avoid solver crash for bad quality meshes

	SET_REF_(pcRefCell,pcRefValue,Orig_,alpha1S_);

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
		Info<<"!std::isfinite(average(pc_)), sgPc_   :"<<min(sgPc_)   <<"    "<<max(sgPc_)<<endl;
		Info<<"!std::isfinite(average(pc_)), sgPcErr_:"<<min(sgPcErr_)<<"   "<<max(sgPcErr_)<<endl;
		Info<<"!std::isfinite(average(pc_)), stf:    "<<min(stf)      <<"   "<<max(stf)<<endl;
		Info<<"!std::isfinite(average(pc_)), alpha1_:"<<min(alpha1_)  <<"   "<<max(alpha1_)<<endl;
		Info<<"pcRefValue:"<<pcRefValue<<endl;
		Info<<"pcRefCell:"<<pcRefCell<<endl;
		Info<<"pc_[pcRefCell]:"<<pc_[pcRefCell]<<endl;
		Info<<"alpha1S_[pcRefCell]:"<<alpha1S_[pcRefCell]<<endl;
		Info<<"alpha1_[pcRefCell]:"<<alpha1_[pcRefCell]<<endl;
		Info<<endl;
		sleep(1);
		std::terminate();
	}
	//ensure(std::isfinite(average(pc_)),"",2);
    
	surfaceScalarField stfTypicalLow=mag((sigma_*0.1)*deltaS_*deltaS_*msh.magSf()); // Note: this filtering is only applied to the interface region, but not away from interface
	surfaceScalarField stfThreshold=((0.5*fcdFilter_)*( mag(stf)+ stfTypicalLow)); // Note: this filtering is only applied to the interface region, but not away from interface
	//sgPcErrn_ = (min (max(sgPc_, -stfThreshold),stfThreshold));
	//sgPc_-=sgPcErrn_;  
    
        
    gPc_= fvc::reconstruct((sgPc_));
    gPc_.correctBoundaryConditions();



    volVectorField vgPcllInterface=(gPc_)-((gPc_)&(nS))*(nS);    //surfaceVectorField gPcf=fvc::interpolate(gPc_,"smoothSchemegPc");

	surfaceVectorField gPcllInterface=fvc::interpolate(vgPcllInterface,"smoothSchemegPc");    //surfaceVectorField gPcf=fvc::interpolate(gPc_,"smoothSchemegPc");
	gPcllInterface=(gPcllInterface)-((gPcllInterface)&(nHatfv))*(nHatfv);    //surfaceVectorField gPcf=fvc::interpolate(gPc_,"smoothSchemegPc");
	gPcllInterface*=(1.-mag(nHatf_)/msh.magSf());
	sgPcErr_=dInterf*(min(2.*(1.4-mag(nHatf_)/msh.magSf()), 1.)*(sgPcErr_+fcCorrectTangent_*(gPcllInterface& msh.Sf())) );//  *( 1.-(mag(nHatf_)/msh.magSf()) )   ); /*aSolid_cPcCorrRelax_* */ 


	//label sgPcErr_checkErrs=0;
	forAll(sgPcErr_,i)
	{
		if (sgPc_[i]*sgPcErr_[i] < -1e-6*stfTypicalLow[i])  sgPcErr_[i]  *= 0.99;
		if( mag(sgPc_[i])<0.5*mag(sgPcErr_[i]) ) {   sgPcErr_[i] = 0.5*sgPc_[i]; } //++sgPcErr_checkErrs;    Info<<"bEC:"<<sgPcErr_[i]<<" "<<sgPc_[i]<<endl; 
	}
	forAll(boundary, bI)
	{
		if (!boundary[bI].coupled())
		{
			sgPcErr_.boundaryFieldRef()[bI] = 0.;
		}
		else
		{
			Field<scalar> & psgPcErr_ = sgPcErr_.boundaryFieldRef()[bI];
			Field<scalar> & psgPc_ = sgPc_.boundaryFieldRef()[bI];
			forAll(psgPcErr_,i)
			{
				if(psgPc_[i]*psgPcErr_[i] < -1e-6*stfTypicalLow[i])  psgPcErr_[i] *= 0.99;
				if( mag(psgPc_[i])<0.5*mag(psgPcErr_[i]) ) {  psgPcErr_[i] = 0.5*psgPc_[i];} // ++sgPcErr_checkErrs;   Info<<"bEC:"<<psgPcErr_[i]<<" "<<psgPc_[i]<<endl;
			}
		}
	}


	//if(sgPcErr_checkErrs)  Info<<"\n\n sgPcErr_checkErrs: "<<sgPcErr_checkErrs<<endl<<endl;

	//stfThreshold=0.5*(stfThreshold+(1.*fcdFilter_)*mag(stf));
	stfThreshold=((1.*fcdFilter_)*( mag(stf)+ mag((sigma_*0.1)*deltaS_*deltaS_*msh.magSf()) ));
	sgPcErrn_ = (min (max(sgPc_, -stfThreshold),stfThreshold) );
	sgPc_-=sgPcErrn_;  

	gPc_= fvc::reconstruct((sgPc_));
	gPc_.correctBoundaryConditions();


	volVectorField gPcS=fvc::average(fvc::interpolate(gPc_)*Internalfaces1_)/fvc::average(Internalfaces1_);
	stfThreshold = 2.*mag((fvc::interpolate(gPcS)) & msh.Sf());
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

