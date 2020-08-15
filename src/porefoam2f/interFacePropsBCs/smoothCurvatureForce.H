{
	volVectorField cellForce
	(
		 IOobject
		 (
			  "cellForce",
			  alpha1_.time().timeName(),
			  alpha1_.mesh()
		 ),
		 alpha1_.mesh(),
		 dimensionedVector("cellForce", curvatureForcef.dimensions(), vector::zero)
	);

	const labelList& owners = mesh.owner();
	const labelList& nybors = mesh.neighbour();

	volScalarField cellWeight=0.0*alpha1_+1.0;
	forAll(nybors,fI) if (magDelS[fI] > CONTRAST_FACTOR)
	{
		cellWeight[owners[fI]]+=100.0*min(1.0, cellWeight[owners[fI]]);
		cellWeight[nybors[fI]]+=100.0*min(1.0, cellWeight[nybors[fI]]);
	}
	
	for(int i=0;i<2;++i)
	{

		//scalar SumForces = 0.0;
		cellForce*=0.0;
		volScalarField interfCellArea=0.0*fvc::average(magSf);



		forAll(nybors,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			label ownr = owners[fI];
			label nybr = nybors[fI];
			cellForce[ownr]+=0.5*cellWeight[ownr]*magDelS[fI]*curvatureForcef[fI];
			cellForce[nybr]+=0.5*cellWeight[nybr]*magDelS[fI]*curvatureForcef[fI];
			interfCellArea[ownr]+=cellWeight[ownr]*magDelS[fI]*magInterfaceSf[fI];
			interfCellArea[nybr]+=cellWeight[nybr]*magDelS[fI]*magInterfaceSf[fI];

		}
		forAll(patches, patchI)
		{
			label pStart=patches[patchI].patch().start();
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[patchI];
			const scalarField & pMagDelS = magDelS.boundaryField()[patchI];
			const vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[patchI];
			if (patches[patchI].coupled())
			 {
				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					label ownr = owners[pStart + pfI];
					cellForce[ownr]+=0.5*cellWeight[ownr]*pMagDelS[pfI]*pCurvatureForcef[pfI];		
					interfCellArea[ownr]+=cellWeight[ownr]*pMagDelS[pfI]*pMagInterfaceSf[pfI];
				}
			 }
			 else		 {			 }
		}

		interfCellArea.correctBoundaryConditions();
		cellForce.correctBoundaryConditions();






		forAll(nybors,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			label ownr = owners[fI];
			label nybr = nybors[fI];
			curvatureForcef[fI] *= 0.01;
			curvatureForcef[fI] += 0.99*cellForce[ownr]*magInterfaceSf[fI]/(interfCellArea[ownr]);
			curvatureForcef[fI] += 0.99*cellForce[nybr]*magInterfaceSf[fI]/(interfCellArea[nybr]);
		}
		//for(int fI=nybors.size(); fI<owners.size();++fI)
		//{
				//curvatureForcef[fI] += 0.5*0.998*cellForce[owners[fI]]*magInterfaceSf[fI]/(interfCellArea[owners[fI]]);
		//}
		forAll(patches, patchI)
		{
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[patchI];
			const scalarField & pMagDelS = magDelS.boundaryField()[patchI];
			if (patches[patchI].coupled())
			{
				vectorField pCellOwnrForce = cellForce.boundaryField()[patchI].patchInternalField();
				vectorField pCellNybrForce = cellForce.boundaryField()[patchI].patchNeighbourField();
				scalarField pCellOwnrArea = interfCellArea.boundaryField()[patchI].patchInternalField();
				scalarField pCellNybrArea = interfCellArea.boundaryField()[patchI].patchNeighbourField();
				vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[patchI];

				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					pCurvatureForcef[pfI] *= 0.01;
					pCurvatureForcef[pfI] += 0.99*pCellOwnrForce[pfI]*pMagInterfaceSf[pfI]/(pCellOwnrArea[pfI]);
					pCurvatureForcef[pfI] += 0.99*pCellNybrForce[pfI]*pMagInterfaceSf[pfI]/(pCellNybrArea[pfI]);
				}

			 }
			 else		 {			 }
		}


		//syncTools::syncFaceList(mesh, curvatureForcef, plusEqOp<vector>(), pTraits<vector>::zero, false);

	}


}



