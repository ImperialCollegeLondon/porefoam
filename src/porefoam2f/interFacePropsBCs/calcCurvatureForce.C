/*-------------------------------------------------------------------------*\

 Interface force computation

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


#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "primitivePatchInterpolation.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvmLaplacian.H"
#include "fvCFD.H"
#include "OFstream.H"
#include "syncTools.H"
#include "correctForContactAngle.H"
#include "pointGradLeastSquare.H"
//#include "pointSnGrad_notUsed.H"
#include "pointInterpolate.H"
//#include "edgeGrad.H"
//#include "pointCurl.H"

#include "smoothOverInterfPoints.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define curtailBADOFSET(a,b) (min (max(a,b),(1.0-(b)))	)





Foam::surfaceScalarField Foam::interfaceProperties::calcCurvatureFConservative
(
		//surfaceVectorField& nSHatfv,
		//const surfaceScalarField& deltaS_
)
{
	const scalar CONTRAST_FACTOR = 1e-4*deltaN_.value();
	const fvMesh& msh = mesh();
	const fvBoundaryMesh& boundary = msh.boundary();
	const fvBoundaryMesh& patches = msh.boundary(); 

	const faceList& faces = msh.faces();
	const pointField& points = msh.points();

	const surfaceVectorField& Sf = msh.Sf();
	//const surfaceVectorField& Cf = msh.Cf();
	const surfaceScalarField&  magSf = msh.magSf();
	const surfaceVectorField SfHat(Sf/magSf);


	// normal surface-gradient of the sharpened alpha
	surfaceScalarField DelSdX = (deltaS_)/msh.deltaCoeffs();
	surfaceScalarField magDelS = mag(deltaS_);











	// volume to point interpolation of alpha
	pointScalarField alpha1P = pointInterpolate(fvc::interpolate(alpha1S_),pMesh_);


	labelList interfPoints(points.size(),0);
	forAll(magDelS,fI)
	{
		if (magDelS[fI] > CONTRAST_FACTOR)
		{
			const face& fac = faces[fI];
			forAll(fac,pI)
			{
					interfPoints[fac[pI]]=1;
			}
		}
	}
	forAll(patches, bI)
	{
		const scalarField & deltaS_p = deltaS_.boundaryField()[bI];
		if (patches[bI].coupled())
		{
			label pStart = patches[bI].patch().start();
			forAll(deltaS_p, pfI) if (mag(deltaS_p[pfI]) > CONTRAST_FACTOR)
			{
			  const face& fac = faces[pStart+pfI];
			  forAll(fac,pI)
			  {
					interfPoints[fac[pI]]=1;
			  }
			}
		}
	}

	syncTools::syncPointList(msh, interfPoints, maxEqOp<label>(), 0, false);

	pointVectorField gradAlphaP = pointGardLeastSquare (fvc::interpolate(alpha1S_), alpha1P, msh.points(),interfPoints);
	//pointVectorField gradAlphaP = pointGardLeastSquare (fvc::interpolate(alpha1S_), alpha1P, msh.points(),interfPoints);
	//pointVectorField gradAlphaP = pointInterpolate(gradAlpha,pMesh_); 
	////////pointVectorField gradAlphaP = pointSnGrad(alpha1S_,U_, interfPoints, pMesh_);
	//vectorField eGradAlpha = edgeGrad(gradAlphaP,nSHatfv,alpha1S_,alpha1P,magDelS, CONTRAST_FACTOR);
	correctForContactAngle(alpha1_,U_,nw_,gradAlphaP,interfPoints);
	//eGradAlpha /= mag(eGradAlpha)+deltaN_.value();
	pointScalarField magGradAlphaP(mag(gradAlphaP));
	magGradAlphaP += 0.1*deltaN_;


	pointVectorField nHatSp = gradAlphaP/magGradAlphaP;
	smoothNSOverInterfPoints(nHatSp, magDelS, msh.magSf(), smoothingKernel_, smoothingRelaxFactor_,interfPoints);
	nHatSp /= mag(nHatSp) + 1e-12;

	#include "calcInterfaceLocation.H"

	#include "calcInterfaceSurfaceArea.H"



	surfaceVectorField curvatureForcef
	(
		IOobject("curvatureForcef", timeName(), msh),
		msh, dimensionedVector("curvatureForcef", dimLength, vector::zero)
	);
	surfaceVectorField nSHatfv // normal to the interface
	(
		IOobject("nSHatfv", timeName(), msh),
		msh, dimensionedVector("nSHatfv", dimless, vector::zero)
	);



{
	// Calculate curvature force at faces
	forAll(magDelS,fI)
	{
		if (magDelS[fI] > CONTRAST_FACTOR)
		{
			const face& fac = faces[fI];

			forAll(fac,pI)
			{
					label ip1(fac[pI]),  ip2(fac.nextLabel(pI));
					vector Le = (points[ip2] - points[ip1]) + (distPointInterface_[ip2] - distPointInterface_[ip1]); // edge vector

					vector Ne = 0.5*(nHatSp[ip2]+nHatSp[ip1]);    Ne /= mag(Ne)+1e-15;


					curvatureForcef[fI] +=  (2.0*(deltaS_[fI]>=0)-1.0) * Le^Ne;
					nSHatfv[fI] +=  gradAlphaP[ip1];
					//nSHatfv[fI] +=  mag(Le) * (Ne);
			}
		}
	}

	forAll(boundary, bI)
	{
		vectorField&  patchCurveIntegralf = curvatureForcef.boundaryField()[bI];
		vectorField&  nSHatfvp = nSHatfv.boundaryField()[bI];
		const scalarField& pMagDelS = magDelS.boundaryField()[bI];
		const scalarField& deltaS_p = deltaS_.boundaryField()[bI];

		forAll(pMagDelS,bfI)
		{
			if (pMagDelS[bfI] > CONTRAST_FACTOR)
			{
				label fI = boundary[bI].patch().start() + bfI;
				const face& fac = faces[fI];

				forAll(fac,pI)
				{
					label ip1(fac[pI]),  ip2(fac.nextLabel(pI));
					vector Le = (points[ip2]-points[ip1])  + (distPointInterface_[ip2]-distPointInterface_[ip1]);

					vector Ne = 0.5*(nHatSp[ip2]+nHatSp[ip1]);   Ne /= mag(Ne)+1e-15;

					patchCurveIntegralf[bfI] += (2.0*(deltaS_p[bfI]>=0)-1.0) * Le^Ne;
					nSHatfvp[bfI] += gradAlphaP[ip1];
					//nHatfp[bfI] += mag(Le) * Ne;
				}
			}
		}
	}

	nSHatfv /= mag(nSHatfv)+1e-15;
	//nHatf = nSHatfv;
	//nSHatfv=nHatf;

}

//Info<<endl<< sum( curvatureForcef*magDelS )<<endl;;


	//#include "smoothCurvatureForce.H"

	surfaceScalarField magInterfaceSf(mag(interfaceSf & nSHatfv)+ 1e-14*magSf);
	#include "smoothFcOverInterfPoints.H"

	#include "smoothCurvatureForce.H"



	//interfPointsOld_=interfPoints;

	return (curvatureForcef & nSHatfv)  /  (magInterfaceSf) *deltaS_*msh.magSf() ;


///.  Kfbsf from balanced total curvature algorithm:(fSf&msh.Sf())/ (mag(msh.Sf()&nSf));
	//return  ((curvatureForcef&msh.Sf()) + (msh.magSf() - mag(nSHatfv&msh.Sf()))*(curvatureForcef&nSHatfv))*deltaS_/ msh.magSf()  ;

	//surfaceVectorField faceCrvForce = ((curvatureForcef & nHatf)  /  (magInterfaceSf))*deltaS_*msh.Sf();
	//#include "preserveTotalForce.H"
	//#include "preserveTotalForce.H"
	//return (faceCrvForce & msh.Sf())  /  (msh.magSf());

}
