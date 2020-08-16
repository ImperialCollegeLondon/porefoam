/*-------------------------------------------------------------------------*\

 Smoothing normal vectors over interface points

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

void Foam::interfaceProperties::smoothOverInterfPointsNotUsed
(
	surfaceVectorField& curvatureForcefDiff, 
	const surfaceScalarField& magDelS, 
	const surfaceScalarField& magInterfaceSf, 
	label smoothingKernel, scalar smoothingRelaxFactor
)
{
	
	const scalar CONTRAST_FACTOR = 1.0e-14*deltaN_.value();
	const fvMesh& msh = alpha1_.mesh();
		const fvBoundaryMesh& patches = msh.boundary(); 

	const faceList& faces = msh.faces();

	const surfaceScalarField&  magSf = msh.magSf();




			//Info<<" * "<< sum( mag(curvatureForcefDiff) ).value()<<" * ";

	
	const dictionary pimple = msh.solutionDict().subDict("PIMPLE");

	pointVectorField pointForce
	(	IOobject( "pointForce", timeName(), msh ),
		 pMesh_,  dimensionedVector("pointForce", curvatureForcefDiff.dimensions(), vector::zero)
	);
	pointScalarField interfPointArea
	(	IOobject( "interfPointArea",  timeName(),  msh ),
		pMesh_,  dimensionedScalar("interfPointArea", magSf.dimensions(), 0.0)
	);



	for(int i=0;i<smoothingKernel;++i)
	{

		pointForce*=0.0;
		interfPointArea*=0.0;


		forAll(magDelS,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			const face& ff = faces[fI];
			forAll(ff, ip)
			{
				label pI = ff[ip];
				pointForce[pI]+=(1.0/ff.size())*curvatureForcefDiff[fI];
				interfPointArea[pI]+=magInterfaceSf[fI];
			}
		}
		forAll(patches, bI)
		{
			label pStart=patches[bI].patch().start();
			const scalarField & pInterfaceSf = magInterfaceSf.boundaryField()[bI];
			const scalarField & pMagDelS = magDelS.boundaryField()[bI];
			const vectorField & pCurvatureForcef = curvatureForcefDiff.boundaryField()[bI];
			if (patches[bI].coupled())
			 {
				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						pointForce[pI]+=(0.5/ff.size())*pCurvatureForcef[pfI];		
						interfPointArea[pI]+=0.5*pInterfaceSf[pfI];
					}
				}
			 }
			 else		 {			 }
		}

		syncTools::syncPointList(msh, pointForce, plusEqOp<vector>(), pTraits<vector>::zero, false);
		syncTools::syncPointList(msh, interfPointArea, plusEqOp<scalar>(), pTraits<scalar>::zero, false);



			//Info<<" *"<< sum( mag(pointForce) ).value()<<"* ";



		forAll(magDelS,fI) if (magDelS[fI] > CONTRAST_FACTOR)
		{
			curvatureForcefDiff[fI] *= 1.0-smoothingRelaxFactor;
			const face& ff = faces[fI];
			forAll(ff, ip)
			{
				label pI = ff[ip];
				curvatureForcefDiff[fI] += smoothingRelaxFactor*pointForce[pI]*magInterfaceSf[fI]/(mag(interfPointArea[pI])+1.0e-36);
			}
		}

		forAll(patches, bI)
		{
			const scalarField & pMagInterfaceSf = magInterfaceSf.boundaryField()[bI];
			const scalarField & pMagDelS = magDelS.boundaryField()[bI];
			if (patches[bI].coupled())
			{
				label pStart=patches[bI].patch().start();
				vectorField & pCurvatureForcef = curvatureForcefDiff.boundaryField()[bI];
				forAll(pMagDelS, pfI) if (pMagDelS[pfI] > CONTRAST_FACTOR)
				{
					pCurvatureForcef[pfI] *= 1.0-smoothingRelaxFactor;
					const face& ff = faces[pStart+pfI];
					forAll(ff, ip)
					{
						label pI = ff[ip];
						pCurvatureForcef[pfI] += smoothingRelaxFactor*pointForce[pI]*pMagInterfaceSf[pfI]/(mag(interfPointArea[pI])+1.0e-36);
						//pCurvatureForcef[pfI] += 1.0*pCellNybrForce[pfI]*pMagInterfaceSf[pfI]/(pCellNybrArea[pfI]);
					}
				}

			 }
			 else		 {			 }
		}

			//Info<<" * "<< sum( mag(curvatureForcefDiff) ).value()<<" * ";

		//syncTools::syncFaceList(msh, curvatureForcefDiff, plusEqOp<vector>(), pTraits<vector>::zero, false);

	}




}



void Foam::interfaceProperties::smoothNSOverInterfPoints
(
	pointVectorField& nHatSp, 
	const surfaceScalarField& magDelS, 
	const surfaceScalarField& magInterfaceSf, 
	label smoothingKernel, scalar smoothingRelaxFactor,
	const labelList& interfPoints
)
{
	
	const scalar CONTRAST_FACTOR = 1.0e-4*deltaN_.value();
	const fvMesh& msh = mesh();
	//const fvBoundaryMesh& patches = msh.boundary(); 
	const fvBoundaryMesh& boundary = msh.boundary();

	const faceList& faces = msh.faces();
	const pointField& points = msh.points();

	//const surfaceScalarField&  magSf = msh.magSf();
	const surfaceVectorField& Cf = msh.Cf();


	const labelListList & faceEdges=msh.faceEdges();


			//Info<<" * "<< sum( mag(nHatSf) ).value()<<" * ";

	
	//const dictionary pimple =msh.solutionDict().subDict("PIMPLE");


	pointVectorField interfFaceNSpS
	(	IOobject( "interfFaceNSpS",  timeName(),  msh ),
		pMesh_,  dimensionedVector("interfFaceArea", dimless, vector::zero),"calculated"
	);
	//pointScalarField interfFaceArea
	//(	IOobject( "interfFaceArea",  timeName(),  msh ),
		//pMesh_,  dimensionedScalar("interfFaceArea", magSf.dimensions(), 0.0),"calculated"
	//);



	for(int i=0;i<smoothingKernel;++i)
	{

    forAll(magDelS,fI)
    {
		  double weightDelS = magDelS[fI];
        if (weightDelS > CONTRAST_FACTOR)
        {
            const face& ff = faces[fI];

 				const labelList& fes = faceEdges[fI];
          forAll(ff,pI)
            {
					label ip1(ff[pI]);
					label ip2(ff.nextLabel(pI));

					point point1 = points[ip1];// + distPointInterface_[ip1];

					point point2 = points[ip2];// + distPointInterface_[ip2];
					
					label bp1(edgemarks_[fes[(-1+pI+ff.size())%ff.size()]]);
					label bp2(edgemarks_[fes[(pI+1)%ff.size()]]);

					vector nHatpEdge = 0.5 * (nHatSp[ip2]*(1.0+2.0*bp2) + nHatSp[ip1]*(1.0+2.0*bp1));

					nHatpEdge /= mag(nHatpEdge) + 1.0e-35;
					
					scalar area =  weightDelS/mag( point2 - point1 );//*mag(( point2 - point1 ) ^ (Cf[fI]-point1)) ;
					//interfFaceArea[ip1] +=  area ;
					//interfFaceArea[ip2] +=  area ;
					interfFaceNSpS[ip1] +=  area*nHatpEdge*mag((Cf[fI]-point2))  ;
					interfFaceNSpS[ip2] +=  area*nHatpEdge*mag((Cf[fI]-point1))  ;
            }
        }
    }
    forAll(boundary, bI)
    {

		if (boundary[bI].coupled())
		{       
			const scalarField& pMagDelS = magDelS.boundaryField()[bI];
        //vectorField& pCurvatureNSf = curvatureNSf.boundaryField()[bI];
        const vectorField& pCf = Cf.boundaryField()[bI];


        forAll(pMagDelS,pfI)
        {
			  double weightDelS = pMagDelS[pfI];
            if ( weightDelS  > CONTRAST_FACTOR)
            {
					label fI = boundary[bI].patch().start() + pfI;
					const face& ff = faces[fI];
					const labelList& fes = faceEdges[fI];
					forAll(ff,pI)
					{
						label ip1(ff[pI]);
						label ip2(ff.nextLabel(pI));

						point point1 = points[ip1];// + distPointInterface_[ip1];

						point point2 = points[ip2];// + distPointInterface_[ip2];

						label bp1(edgemarks_[fes[(-1+pI+ff.size())%ff.size()]]);
						label bp2(edgemarks_[fes[(pI+1)%ff.size()]]);

						vector nHatpEdge = 0.5 * (nHatSp[ip2]*(1.0+2.0*bp2) + nHatSp[ip1]*(1.0+2.0*bp1));

						nHatpEdge /= mag(nHatpEdge)+1.0e-35;
					
						scalar area =  weightDelS*0.5/mag( point2 - point1 );//*mag(( point2 - point1 ) ^ (pCf[pfI]-point1)) ;
						//interfFaceArea[ip1] +=  area ;
						//interfFaceArea[ip2] +=  area ;
						interfFaceNSpS[ip1] +=  area*nHatpEdge*mag((pCf[pfI]-point2)) ;
						interfFaceNSpS[ip2] +=  area*nHatpEdge*mag((pCf[pfI]-point1)) ;
					}

            }
        }
		}
    }

		//syncTools::syncPointList(msh, interfFaceArea, plusEqOp<scalar>(), pTraits<scalar>::zero, false);
		syncTools::syncPointList(msh, interfFaceNSpS, plusEqOp<vector>(), pTraits<vector>::zero, false);





	nHatSp = (1.0-smoothingRelaxFactor)*nHatSp + smoothingRelaxFactor*interfFaceNSpS/ (mag(interfFaceNSpS) + 1.0e-18);
	
	nHatSp/=mag(nHatSp)+1.0e-14;
	correctForContactAngle(alpha1_,U_,nw_,nHatSp,interfPoints);


	}




}

