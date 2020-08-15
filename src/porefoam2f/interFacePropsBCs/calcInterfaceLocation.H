
///////////Interface reconstruction based on cellToPoint leastSquare interlolation////

			//distPointInterface_ =	dimensionedVector("distPointInterface_",dimLength,vector::zero);
	//forAll(interfPoints,pi) if(!interfPoints[pi]) nSpCurl_[pi]*=0.0; ;;


	//pointVectorField Up = pointInterpolate(U_,pMesh_);
	//nSpCurl_=0.1*nHatSp^pointCurl(nHatSp,Up,fvc::interpolate(U_),alpha1_.mesh().time().deltaT().value());
	//nSpCurl_/=0.1*mag(nSpCurl_)+0.99*mag(nHatSp)+0.01;
	//nHatSp += nSpCurl_;
	//nHatSp /= mag(nHatSp) +1.0e-15;

//nSpCurl_=nHatSp;



	// calculate the distance between interface and point
	forAll(interfPoints,pi)
	{ if(interfPoints[pi])
		{
			//if(interfPointsOld_[pi])
			//{
				//vector disp = Up[pi]*mesh.time().deltaT().value() ;
				//disp-=1.2*(disp-(disp&nHatSp[pi])*nHatSp[pi]);
				//distPointInterface_[pi] = (distPointInterface_[pi]+disp)*0.9 + 0.1*((scalar(0.5) - alpha1P[pi])/(magGradAlphaP[pi]))*nHatSp[pi];
			//}
			//else
			{
				scalar ap=alpha1P[pi];
				//distPointInterface_[pi] = min(1.0,max(5.0*ap*(1.0-ap)-0.1,0.0))*((scalar(0.5) - ap)/(magGradAlphaP[pi]))*nHatSp[pi];
				distPointInterface_[pi] = ((0.0- ap)/(magGradAlphaP[pi]))*nHatSp[pi];
			}
		}
		else distPointInterface_[pi]*=0.0;
			
	}

	// calculate the distance between interface and point
	  //distPointInterface_ = ((scalar(0.5) - alpha1P)/(magGradAlphaP))*nHatSp;    	//distPointInterface_ = (((1.0/0.9/0.9)*(alpha1P+0.4)*(1.4-alpha1P)))*((0.5 - alpha1P)/(magGradAlphaP))*nHatSp;






		///. delete gap between interface and wall
		const polyBoundaryMesh& boundaryMesh = mesh.boundaryMesh();
		//const fvBoundaryMesh& patches = mesh.boundary();
		forAll(boundaryMesh, patchi)
		{
			const labelList& patchPoints = boundaryMesh[patchi].meshPoints();
			if (patches[patchi].coupled()) {}
			else	if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[patchi]))
			{
				primitivePatchInterpolation pinterpolator(boundaryMesh[patchi]);
				pointField nWp= pinterpolator.faceToPointInterpolate(nw_.boundaryField()[patchi]);
				nWp/=mag(nWp)+1.0e-15;
				forAll(patchPoints, ppI)
				{
					label pI = patchPoints[ppI];
					vector ns = nHatSp[pI];
					vector ss =  ns -  (ns & nWp[ppI])*nWp[ppI];  ss/=mag(ss)+1.0e-6;
					const scalar cosNsSs = ns&ss;
					distPointInterface_[pI] = 1.1*((distPointInterface_[pI] & ss)/(cosNsSs+0.1)) * ss;
				}
			}
			else forAll(patchPoints, ppI) distPointInterface_[patchPoints[ppI]]=vector::zero;
			
		}


