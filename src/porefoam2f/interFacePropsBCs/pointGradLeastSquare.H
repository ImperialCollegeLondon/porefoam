namespace Foam
{
	Foam::pointVectorField pointGardLeastSquare
	(
		const surfaceScalarField& alphaf,
		const pointScalarField alpha1P,
		const pointField points,
		const labelList& interfPoints
	)
	{
		const fvMesh& mesh = alphaf.mesh();
		//const labelListList& pointCells = mesh.pointCells();
		const pointMesh &  pMesh = alpha1P.mesh();

		pointVectorField gradAlphaP
		(
		IOobject  (  "gradAlphaP",  mesh.time().timeName(),  mesh  ),
			pMesh,
			dimensionedVector("gradAlphaP",dimless/dimLength,vector::zero)
		);


		//surfaceScalarField alphaf = fvc::interpolate(alpha1);
		const faceList &  faces=mesh.faces();
		const surfaceVectorField & Cf=mesh.Cf();
		//const surfaceVectorField & deltaCC_=mesh.delta();
		//const vectorField& Cc=mesh.C();
		const surfaceScalarField & magSf=mesh.magSf();
		//const surfaceScalarField & snGAlpha=fvc::snGrad(alpha1);
		//const labelList& owners = mesh.faceOwner();
		//const labelList& neibors = mesh.faceNeighbour();
		//  const pointField &  points=mesh.points();
		//  vectorField gradAlphaP(points.size(),vector::zero);

		// normal surface-gradient of the sharpened alpha
		//surfaceScalarField deltaS_(fvc::snGrad(alpha1Sharp));
		//surfaceScalarField DelSdX = (deltaS_)/mesh.deltaCoeffs();
		//surfaceScalarField magdS = mag(DelSdX);



		symmTensorField GG(points.size(),symmTensor(1.0e-36,0.0,0.0,1.0e-36,0.0,1.0e-36));



		forAll(Cf, fI) // if(magdS[fI]>1.0e-12)
		{
			vector Cpi = Cf[fI];
			scalar weight=magSf[fI];  weight*=weight;// *magdS[fI]
			const face& ff = faces[fI];
			forAll(ff, ip)
			{	label pI = ff[ip];
				vector pp = points[pI];
				scalar weightpd2 = weight/magSqr(pp-Cpi); //weightpd2*=weightpd2;
				GG[pI] += weightpd2*sqr((Cpi - pp));
			}
		}


		const fvBoundaryMesh& patches = mesh.boundary();
		forAll(patches, patchI)
		{
			label pStart=patches[patchI].patch().start();
			const scalarField & pmagSf = magSf.boundaryField()[patchI];
			const vectorField & Cfp = Cf.boundaryField()[patchI];
			if (patches[patchI].coupled())
			{
				forAll(Cfp, pfI) // if(pmagdS[pfI]>1.0e-12)
				{
					vector Cpi = Cfp[pfI];
					scalar weight=pmagSf[pfI];   weight*=weight*0.5; //*pmagdS[pfI]
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{	vector pp = points[ff[ip]];
						scalar weightpd2 = weight/magSqr(pp-Cpi); //weightpd2*=weightpd2;
						GG[ff[ip]] += weightpd2*sqr((Cpi - pp));
					}
				}
			}
			else
			{
				forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
				{
					vector Cpi = Cfp[pfI];
					scalar weight=pmagSf[pfI];  weight*=weight;  //*pmagdS[pfI]
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{	vector pp = points[ff[ip]];
						scalar weightpd2 = weight/magSqr(pp-Cpi); //weightpd2*=weightpd2;
						GG[ff[ip]] += weightpd2*sqr((Cpi - pp));
					}
				}
			}
		}

		syncTools::syncPointList(mesh, GG, plusEqOp<symmTensor>(), pTraits<symmTensor>::zero, false);





		  symmTensorField invGG(points.size(),symmTensor::zero);
		  forAll(points,pI)
		  {
				 if(interfPoints[pI])
					  invGG[pI] = inv(GG[pI]);
		  }


		forAll(alphaf, fI) // if(magdS[fI]>1.0e-12)
		{
			scalar weight=magSf[fI];  weight*=weight;// *magdS[fI]
			vector Cfi = Cf[fI];
			const face& ff = faces[fI];
			forAll(ff, ip)
			{	label pI = ff[ip];
				scalar weightpd2 = weight/magSqr(points[pI]-Cfi); //weightpd2*=weightpd2;
				gradAlphaP[pI] += weightpd2*(invGG[pI]&(Cfi - points[pI]))*(alphaf[fI] - alpha1P[pI]);
			}


		}


		forAll(patches, patchI)
		{
			label pStart=patches[patchI].patch().start();
			const scalarField & pmagSf = magSf.boundaryField()[patchI];
			const scalarField & alphafp = alphaf.boundaryField()[patchI];
			const vectorField & Cfp = Cf.boundaryField()[patchI];
			if (patches[patchI].coupled())
			 {
				forAll(Cfp, pfI) // if(pmagdS[pfI]>1.0e-12)
				{
					scalar weight=pmagSf[pfI];weight*=weight*0.5;//*pmagdS[pfI]
					scalar alphafi=alphafp[pfI] ;
					vector Cfi = Cfp[pfI];
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{	label pI = ff[ip];
						scalar weightpd2 = weight/magSqr(points[pI]-Cfi); //weightpd2*=weightpd2;
						gradAlphaP[pI] += weightpd2*(invGG[pI]&(Cfi - points[pI]))*(alphafi - alpha1P[pI]);
					}
				}
			 }
			 else
			 {
				forAll(Cfp, pfI)  //if (pmagdS[pfI]>1.0e-8)
				{
					scalar weight=pmagSf[pfI]; weight*=weight;//*pmagdS[pfI]
					scalar alphafi=alphafp[pfI] ;
					vector Cfi = Cfp[pfI];
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{	label pI = ff[ip];
						scalar weightpd2 = weight/magSqr(points[pI]-Cfi); //weightpd2*=weightpd2;
						gradAlphaP[pI] += weightpd2*(invGG[pI]&(Cfi - points[pI]))*(alphafi - alpha1P[pI]);
					}
				}
			 }
		  }

		syncTools::syncPointList(mesh, gradAlphaP, plusEqOp<vector>(), pTraits<vector>::zero, false);





		  return gradAlphaP;
	 }





}
