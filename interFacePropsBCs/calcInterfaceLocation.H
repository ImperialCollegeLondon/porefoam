/*-------------------------------------------------------------------------*\

 Interface reconstruction using gradient of alpha

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

//! This file is part of calcCurvatureForce.C



	//distPointInterface_ =	dimensionedVector("distPointInterface_",dimLength,vector::zero);
	//forAll(interfPoints,pi) if(!interfPoints[pi]) nSpCurl_[pi]*=0.; ;;


	//pointVectorField Up = pointInterpolate(U_,pMesh_);
	//nSpCurl_=0.1*nHatSp^pointCurl(nHatSp,Up,fvc::interpolate(U_),msh.time().deltaT().value());
	//nSpCurl_/=0.1*mag(nSpCurl_)+0.99*mag(nHatSp)+0.01;
	//nHatSp += nSpCurl_;
	//nHatSp /= mag(nHatSp) +1e-15;

	//nSpCurl_=nHatSp;



	// calculate the distance between interface and point
	forAll(interfPoints,pi)
	{ if(interfPoints[pi])
		{
			//if(interfPointsOld_[pi])
			//{
				//vector disp = Up[pi]*msh.time().deltaT().value() ;
				//disp-=1.2*(disp-(disp&nHatSp[pi])*nHatSp[pi]);
				//distPointInterface_[pi] = (distPointInterface_[pi]+disp)*0.9 + 0.1*((scalar(0.5) - alpha1P[pi])/(magGradAlphaP[pi]))*nHatSp[pi];
			//}
			//else
			{
				scalar ap=alpha1P[pi];
				//distPointInterface_[pi] = min(1.,max(5.*ap*(1.-ap)-0.1,0.))*((scalar(0.5) - ap)/(magGradAlphaP[pi]))*nHatSp[pi];
				distPointInterface_[pi] = ((0.- ap)/(magGradAlphaP[pi]))*nHatSp[pi];
			}
		}
		else distPointInterface_[pi]*=0.;

	}

	// calculate the distance between interface and point
	  //distPointInterface_ = ((scalar(0.5) - alpha1P)/(magGradAlphaP))*nHatSp;    	//distPointInterface_ = (((1./0.9/0.9)*(alpha1P+0.4)*(1.4-alpha1P)))*((0.5 - alpha1P)/(magGradAlphaP))*nHatSp;






		///. delete gap between interface and wall
		const polyBoundaryMesh& boundaryMesh = msh.boundaryMesh();
		//const fvBoundaryMesh& patches = msh.boundary();
		forAll(boundaryMesh, bi)
		{
			const labelList& patchPoints = boundaryMesh[bi].meshPoints();
			if (patches[bi].coupled()) {}
			else	if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[bi]))
			{
				primitivePatchInterpolation pinterpolator(boundaryMesh[bi]);
				pointField nWp= pinterpolator.faceToPointInterpolate(nw_.boundaryField()[bi]);
				nWp/=mag(nWp)+1e-15;
				forAll(patchPoints, ppI)
				{
					label pI = patchPoints[ppI];
					vector ns = nHatSp[pI];
					vector ss =  ns -  (ns & nWp[ppI])*nWp[ppI];  ss/=mag(ss)+1e-6;
					const scalar cosNsSs = ns&ss;
					distPointInterface_[pI] = 1.1*((distPointInterface_[pI] & ss)/(cosNsSs+0.1)) * ss;
				}
			}
			else forAll(patchPoints, ppI) distPointInterface_[patchPoints[ppI]]=vector::zero;

		}
