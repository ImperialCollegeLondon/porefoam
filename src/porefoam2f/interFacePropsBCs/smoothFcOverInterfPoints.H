{
	const dictionary pimple =alpha1_.mesh().solutionDict().subDict("PIMPLE");

    smoothingKernel_=   readLabel( pimple.lookup("smoothingKernel") )%10;
    smoothingRelaxFactor_=   readScalar( pimple.lookup("smoothingRelaxFactor") );

	pointVectorField pointForce
	(	IOobject( "pointForce", alpha1_.time().timeName(), alpha1_.mesh() ),
		 pMesh_,  dimensionedVector("pointForce", curvatureForcef.dimensions(), vector::zero)
	);
	pointScalarField interfPointArea
	(	IOobject( "interfPointArea",  alpha1_.time().timeName(),  alpha1_.mesh() ),
		pMesh_,  dimensionedScalar("interfPointArea", magSf.dimensions(), 0.0)
	);

//Info<<endl<< sum( mag(curvatureForcef*magDelS) ).value()<<" sfpfs ";;


	for(int i=0;i<smoothingKernel_;++i)
	{

		pointForce*=0.0;
		interfPointArea*=0.0;


		forAll(magDelS,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			const face& ff = faces[fI];
			forAll(ff, ip)
			{
				label pI = ff[ip];
				pointForce[pI]+=(1.0/ff.size())*magDelS[fI]*curvatureForcef[fI];
				interfPointArea[pI]+=magDelS[fI]*magInterfaceSf[fI];
			}
		}
		forAll(patches, patchI)
		{
			label pStart=patches[patchI].patch().start();
			const scalarField & pInterfaceSf = magInterfaceSf.boundaryField()[patchI];
			const scalarField & pMagDelS = magDelS.boundaryField()[patchI];
			const vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[patchI];
			if (patches[patchI].coupled())
			 {
				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						pointForce[pI]+=(0.5/ff.size())*pMagDelS[pfI]*pCurvatureForcef[pfI];		
						interfPointArea[pI]+=0.5*pMagDelS[pfI]*pInterfaceSf[pfI];
					}
				}
			 }
			 else		 {			 }
		}

		syncTools::syncPointList(mesh, pointForce, plusEqOp<vector>(), pTraits<vector>::zero, false);
		syncTools::syncPointList(mesh, interfPointArea, plusEqOp<scalar>(), pTraits<scalar>::zero, false);






		forAll(magDelS,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			curvatureForcef[fI] *= 1.0-smoothingRelaxFactor_;
			const face& ff = faces[fI];
			forAll(ff, ip)
			{
				label pI = ff[ip];
				curvatureForcef[fI] += smoothingRelaxFactor_*pointForce[pI]*magInterfaceSf[fI]/(mag(interfPointArea[pI])+1.0e-36);
			}
		}

		forAll(patches, patchI)
		{
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[patchI];
			const scalarField & pMagDelS = magDelS.boundaryField()[patchI];
			if (patches[patchI].coupled())
			{
				label pStart=patches[patchI].patch().start();
				vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[patchI];
				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					pCurvatureForcef[pfI] *= 1.0-smoothingRelaxFactor_;
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						pCurvatureForcef[pfI] += smoothingRelaxFactor_*pointForce[pI]*pMagInterfaceSf[pfI]/(mag(interfPointArea[pI])+1.0e-36);
						//pCurvatureForcef[pfI] += 1.0*pCellNybrForce[pfI]*pMagInterfaceSf[pfI]/(pCellNybrArea[pfI]);
					}
				}

			 }
			 else		 {			 }
		}


		//syncTools::syncFaceList(mesh, curvatureForcef, plusEqOp<vector>(), pTraits<vector>::zero, false);

	}


//Info<< sum( mag(curvatureForcef*magDelS) ).value()<<endl;;


}

