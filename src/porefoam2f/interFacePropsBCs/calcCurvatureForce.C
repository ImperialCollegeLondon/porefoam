/*-------------------------------------------------------------------------*\
This code is part of poreFOAM, a suite of codes written using OpenFOAM
for direct simulation of flow at the pore scale. 
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ali Qaseminejad Raeini as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:	 a.q.raeini@imperial.ac.uk
Mosayeb Shams:	m.shams14@imperial.ac.uk
Branko Bijeljic:  b.bijeljic@imperial.ac.uk
Martin J Blunt:   m.blunt@imperial.ac.uk

Description:
	interface force and capillary pressure computation

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
	const scalar CONTRAST_FACTOR = 1.0e-4*deltaN_.value();
	const fvMesh& mesh = alpha1_.mesh();
	const fvBoundaryMesh& boundary = mesh.boundary();
	const fvBoundaryMesh& patches = mesh.boundary(); 

	const faceList& faces = mesh.faces();
	const pointField& points = mesh.points();

	const surfaceVectorField& Sf = mesh.Sf();
	//const surfaceVectorField& Cf = mesh.Cf();
	const surfaceScalarField&  magSf = mesh.magSf();
	const surfaceVectorField SfHat(Sf/magSf);


	// normal surface-gradient of the sharpened alpha
	surfaceScalarField DelSdX = (deltaS_)/mesh.deltaCoeffs();
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

	syncTools::syncPointList(mesh, interfPoints, maxEqOp<label>(), 0, false);

	pointVectorField gradAlphaP = pointGardLeastSquare (fvc::interpolate(alpha1S_), alpha1P, mesh.points(),interfPoints);
	//pointVectorField gradAlphaP = pointGardLeastSquare (fvc::interpolate(alpha1S_), alpha1P, mesh.points(),interfPoints);
	//pointVectorField gradAlphaP = pointInterpolate(gradAlpha,pMesh_); 
	////////pointVectorField gradAlphaP = pointSnGrad(alpha1S_,U_, interfPoints, pMesh_);
	//vectorField eGradAlpha = edgeGrad(gradAlphaP,nSHatfv,alpha1S_,alpha1P,magDelS, CONTRAST_FACTOR);
	correctForContactAngle(alpha1_,U_,nw_,gradAlphaP,interfPoints);
	//eGradAlpha /= mag(eGradAlpha)+deltaN_.value();
	pointScalarField magGradAlphaP(mag(gradAlphaP));
	magGradAlphaP += 0.1*deltaN_;


	pointVectorField nHatSp = gradAlphaP/magGradAlphaP;
	smoothNSOverInterfPoints(nHatSp, magDelS, mesh.magSf(), smoothingKernel_, smoothingRelaxFactor_,interfPoints);
	nHatSp /= mag(nHatSp) + 1.0e-12;

	#include "calcInterfaceLocation.H"

	#include "calcInterfaceSurfaceArea.H"



	surfaceVectorField curvatureForcef
	(
		IOobject("curvatureForcef", alpha1_.time().timeName(), alpha1_.mesh()),
		alpha1_.mesh(), dimensionedVector("curvatureForcef", dimLength, vector::zero)
	);
	surfaceVectorField nSHatfv // normal to the interface
	(
		IOobject("nSHatfv", alpha1_.time().timeName(), alpha1_.mesh()),
		alpha1_.mesh(), dimensionedVector("nSHatfv", dimless, vector::zero)
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

					vector Ne = 0.5*(nHatSp[ip2]+nHatSp[ip1]);    Ne /= mag(Ne)+1.0e-15;


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

					vector Ne = 0.5*(nHatSp[ip2]+nHatSp[ip1]);   Ne /= mag(Ne)+1.0e-15;

					patchCurveIntegralf[bfI] += (2.0*(deltaS_p[bfI]>=0)-1.0) * Le^Ne;
					nSHatfvp[bfI] += gradAlphaP[ip1];
					//nHatfp[bfI] += mag(Le) * Ne;
				}
			}
		}
	}

	nSHatfv /= mag(nSHatfv)+1.0e-15;
	//nHatf = nSHatfv;
	//nSHatfv=nHatf;

}

//Info<<endl<< sum( curvatureForcef*magDelS )<<endl;;


	//#include "smoothCurvatureForce.H"

	surfaceScalarField magInterfaceSf(mag(interfaceSf & nSHatfv)+ 1.0e-14*magSf);
	#include "smoothFcOverInterfPoints.H"

	#include "smoothCurvatureForce.H"



	//interfPointsOld_=interfPoints;

	return (curvatureForcef & nSHatfv)  /  (magInterfaceSf) *deltaS_*mesh.magSf() ;


///.  Kfbsf from balanced total curvature algorithm:(fSf&mesh.Sf())/ (mag(mesh.Sf()&nSf));
	//return  ((curvatureForcef&mesh.Sf()) + (mesh.magSf() - mag(nSHatfv&mesh.Sf()))*(curvatureForcef&nSHatfv))*deltaS_/ mesh.magSf()  ;

	//surfaceVectorField faceCrvForce = ((curvatureForcef & nHatf)  /  (magInterfaceSf))*deltaS_*mesh.Sf();
	//#include "preserveTotalForce.H"
	//#include "preserveTotalForce.H"
	//return (faceCrvForce & mesh.Sf())  /  (mesh.magSf());

}
