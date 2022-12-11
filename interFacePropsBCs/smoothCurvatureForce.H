/*-------------------------------------------------------------------------*\

 Curvature smoothing over interface faces

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
	volVectorField cellForce
	(
		 IOobject("cellForce", timeName(), msh ),
		 msh, dimensionedVector("cellForce", curvatureForcef.dimensions(), vector::zero)
	);

	const labelList& owners = msh.owner();
	const labelList& nybors = msh.neighbour();

	volScalarField cellWeight=0.*alpha1_+1.;
	forAll(nybors,fI) if (magDelS[fI] > CONTRAST_FACTOR)
	{
		cellWeight[owners[fI]]+=100.*min(1., cellWeight[owners[fI]]);
		cellWeight[nybors[fI]]+=100.*min(1., cellWeight[nybors[fI]]);
	}

	for(int i=0;i<2;++i)
	{

		//scalar SumForces = 0.;
		cellForce*=dimensionedScalar(0.);
		volScalarField interfCellArea=0.*fvc::average(magSf);



		forAll(nybors,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			label ownr = owners[fI];
			label nybr = nybors[fI];
			cellForce[ownr]+=0.5*cellWeight[ownr]*magDelS[fI]*curvatureForcef[fI];
			cellForce[nybr]+=0.5*cellWeight[nybr]*magDelS[fI]*curvatureForcef[fI];
			interfCellArea[ownr]+=cellWeight[ownr]*magDelS[fI]*magInterfaceSf[fI];
			interfCellArea[nybr]+=cellWeight[nybr]*magDelS[fI]*magInterfaceSf[fI];

		}
		forAll(patches, bI)
		{
			label pStart=patches[bI].patch().start();
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[bI];
			const scalarField & pMagDelS = magDelS.boundaryField()[bI];
			const vectorField & pCurvatureForcef = curvatureForcef.boundaryField()[bI];
			if (patches[bI].coupled())
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
		forAll(patches, bI)
		{
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[bI];
			const scalarField & pMagDelS = magDelS.boundaryField()[bI];
			if (patches[bI].coupled())
			{
				vectorField pCellOwnrForce = cellForce.boundaryField()[bI].patchInternalField();
				vectorField pCellNybrForce = cellForce.boundaryField()[bI].patchNeighbourField();
				scalarField pCellOwnrArea = interfCellArea.boundaryField()[bI].patchInternalField();
				scalarField pCellNybrArea = interfCellArea.boundaryField()[bI].patchNeighbourField();
				vectorField & pCurvatureForcef = curvatureForcef.boundaryFieldRef()[bI];

				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					pCurvatureForcef[pfI] *= 0.01;
					pCurvatureForcef[pfI] += 0.99*pCellOwnrForce[pfI]*pMagInterfaceSf[pfI]/(pCellOwnrArea[pfI]);
					pCurvatureForcef[pfI] += 0.99*pCellNybrForce[pfI]*pMagInterfaceSf[pfI]/(pCellNybrArea[pfI]);
				}

			 }
			 else		 {			 }
		}


		//syncTools::syncFaceList(msh, curvatureForcef, plusEqOp<vector>(), pTraits<vector>::zero, false);

	}


}
