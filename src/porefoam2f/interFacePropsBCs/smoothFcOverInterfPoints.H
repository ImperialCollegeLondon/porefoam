/*-------------------------------------------------------------------------*\

 Curvature force smoothing over interface points

 Copyright (C) 2014-2020  Mosayeb Shams
 Copyright (C) 2017-2020  Ali Qaseminejad Raeini 

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/

{
	const dictionary pimple =msh.solutionDict().subDict("PIMPLE");

    smoothingKernel_=   readLabel( pimple.lookup("smoothingKernel") )%10;
    smoothingRelaxFactor_=   readScalar( pimple.lookup("smoothingRelaxFactor") );

	pointVectorField pointForce
	(	IOobject( "pointForce", timeName(), msh ),
		 pMesh_,  dimensionedVector("pointForce", curvatureForcef.dimensions(), vector::zero)
	);
	pointScalarField interfPointArea
	(	IOobject( "interfPointArea",  timeName(),  msh ),
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
		forAll(patches, bI)
		{
			label pStart=patches[bI].patch().start();
			const scalarField & pInterfaceSf = magInterfaceSf.boundaryField()[bI];
			const scalarField & pMagDelS = magDelS.boundaryField()[bI];
			const vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[bI];
			if (patches[bI].coupled())
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

		syncTools::syncPointList(msh, pointForce, plusEqOp<vector>(), pTraits<vector>::zero, false);
		syncTools::syncPointList(msh, interfPointArea, plusEqOp<scalar>(), pTraits<scalar>::zero, false);






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

		forAll(patches, bI)
		{
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[bI];
			const scalarField & pMagDelS = magDelS.boundaryField()[bI];
			if (patches[bI].coupled())
			{
				label pStart=patches[bI].patch().start();
				vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[bI];
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


		//syncTools::syncFaceList(msh, curvatureForcef, plusEqOp<vector>(), pTraits<vector>::zero, false);

	}


//Info<< sum( mag(curvatureForcef*magDelS) ).value()<<endl;;


}

