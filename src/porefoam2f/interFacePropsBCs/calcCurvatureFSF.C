/*-------------------------------------------------------------------------*\
 Compute surface force using CSF/SSF/FCF algorithms

 Copyright (C) 2010-2020  Ali Qaseminejad Raeini 

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

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
    surfaceVectorField::Boundary& nsb,
    volVectorField::Boundary& gradAlphab,
    volVectorField::Boundary& nSb,
    volScalarField::Boundary& alphaSb
) const
{
	const fvMesh& msh = mesh();
	const volScalarField::Boundary& abf = alpha1_.boundaryField();

	const fvBoundaryMesh& boundary = msh.boundary();

	forAll(boundary, bI)
	{
		if (isA<alphaContactAngleFvPatchScalarField>(abf[bI]))
		{
			alphaContactAngleFvPatchScalarField& acap =
				 const_cast<alphaContactAngleFvPatchScalarField&>
				 (
					  refCast<const alphaContactAngleFvPatchScalarField>
					  (abf[bI])
				 );

			

			gradAlphab[bI] = gradAlphab[bI].patchInternalField();
			gradAlphab[bI] == gradAlphab[bI].patchInternalField();

			vectorField gAlphaS = gradAlphab[bI];
			primitivePatchInterpolation pinterpolator(msh.boundaryMesh()[bI]);
			for (int i=0; i<1;i++)
			{	
				gAlphaS = 0.1*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(gAlphaS) )  //+0.5*alpha2sqrtp*gradAlphab[bI];
							 + 0.9*gradAlphab[bI];
			}
			vectorField nsp = gAlphaS/(mag(gAlphaS) + deltaN_.value());


			
			//vectorField nsp0 = nsb[bI];
			vectorField nf = nw_.boundaryField()[bI];
			scalarField theta =
			convertToRad*acap.theta(U_.boundaryField()[bI], nsp, nf);

			// vectorField nf = boundary[bI].nf();
			// Reset nsp to correspond to the contact angle




			vectorField ss=nsp-(nsp&nf)*nf;
			ss/=mag(ss)+1e-8;
			nsp=sin(theta)*ss+cos(theta)*nf;


			nSb[bI] == nsp;
			nsb[bI] == nsp;

			acap.gradient() == 0.5*(1.-acap*acap)*((boundary[bI].nf() & nsp)*mag(gradAlphab[bI]));

			alphaContactAngleFvPatchScalarField& alphaSbCap =
				 const_cast<alphaContactAngleFvPatchScalarField&>
				 ( refCast<const alphaContactAngleFvPatchScalarField>(alphaSb[bI]) );            
			alphaSbCap.gradient() = (boundary[bI].nf() & nsp)*mag(gAlphaS);
			alphaSbCap.gradient() == (boundary[bI].nf() & nsp)*mag(gAlphaS);
			gradAlphab[bI] == nsp*mag(gradAlphab[bI]);
			alphaSbCap.evaluate();

		}
	}
}




void Foam::interfaceProperties::calcCurvatureFSF
(
		surfaceScalarField&  stf,
		const surfaceScalarField& delS,
		const volScalarField& a1a2
)
{

    const fvMesh& msh = stf.mesh();
    const surfaceVectorField& Sf = msh.Sf();
    
	const dictionary pimple = msh.solutionDict().subDict("PIMPLE");

	smoothingKernel_ = readLabel( pimple.lookup("smoothingKernel") )%10;
    smoothingRelaxFactor_ = readScalar( pimple.lookup("smoothingRelaxFactor") );


	volVectorField gradAlpha = fvc::reconstruct(fvc::snGrad(alpha1S_)*msh.magSf());//fvc::grad(alpha1S_,"smoothScheme");
	//volVectorField gradAlpha = fvc::grad(alpha1S_);//fvc::grad(alpha1S_,"smoothScheme");
	gradAlpha.correctBoundaryConditions();
	volScalarField magGradAlpha = mag(gradAlpha) + deltaN_; 

	volVectorField nS_ = a1a2*gradAlpha/magGradAlpha; nS_ = nS_/(mag(nS_) + 1e-8);
	surfaceVectorField nHatfv = fvc::interpolate(nS_,"smoothScheme"); nHatfv/=mag(nHatfv)+1e-12;
	correctContactAngle(nHatfv.boundaryField(),gradAlpha.boundaryField(),nS_.boundaryField(), alpha1S_.boundaryField());    




	volScalarField a1a2Relaxed = smoothingRelaxFactor_*a1a2;
	volScalarField a1xa2 = a1a2*(1.-0.000001)+0.000001;




	for (int i=0; i<smoothingKernel_;i++) 
	{	
		nS_=(1.-a1a2Relaxed)*nS_+a1a2Relaxed*fvc::average(fvc::interpolate(a1a2Relaxed*nS_,"smoothScheme"));	
		nS_.correctBoundaryConditions();
		nS_=(1.-a1a2Relaxed)*nS_+a1a2Relaxed*fvc::average(fvc::interpolate(nS_));	
		nS_.correctBoundaryConditions();

		nHatfv = fvc::interpolate(nS_,"smoothScheme");
		nHatfv = nHatfv/(mag(nHatfv) + 1e-12);
		correctContactAngle(nHatfv.boundaryField(),gradAlpha.boundaryField(),nS_.boundaryField(),alpha1S_.boundaryField()	);

	}
	nS_ = nS_/(mag(nS_) + 1e-12);




	nHatf_ = nHatfv & Sf;

	surfaceScalarField snGradAlpha = fvc::snGrad(alpha1S_);
	forAll(nHatf_,i)  if (nHatf_[i]*snGradAlpha[i] < 0) nHatf_[i]*=0.9;


	volScalarField K_ = -fvc::div(nHatf_);



	nHatfv = fvc::interpolate(nS_);  // the scheme is different from above
	nHatf_ = nHatfv & Sf;
	forAll(nHatf_,i)   if (nHatf_[i]*snGradAlpha[i] < 0) nHatf_[i]*=-0.1; 




	if (smoothingKernel_) ///. smoothing Kc
	{
		// this correction helps stablizing also improves the accuracy for capillary pressure,
		// it has some theoretical basis but the coefficients here are chosen emperically
		volScalarField mgK=mag(K_);
		K_/=mag( 1. + 0.1*((alpha1S_-0.5)+0.1*K_/(magGradAlpha+5.*mgK)) * K_/(magGradAlpha+2.*mgK) )+1e-12;

		{
			 K_.correctBoundaryConditions();
			 surfaceScalarField WKf=0.02+0.08*fvc::interpolate(a1xa2)+(mag(delS)/(mag(delS)+deltaN_));
			 K_ = fvc::average(fvc::interpolate(K_)*WKf )/fvc::average(WKf );
			 K_ = fvc::average(fvc::interpolate(K_,"smoothSchemeK")*WKf )/fvc::average(WKf );

			 K_.correctBoundaryConditions();
		}

		stf=(sigma_)*((fvc::interpolate(K_*a1xa2,"smoothSchemeK"))/fvc::interpolate(a1xa2,"smoothSchemeK"))  *delS*msh.magSf(); 

	}
	else
	{	
		stf=sigma_*(fvc::interpolate(K_))*delS*msh.magSf(); 
		Info<<"SK:"<<0<<endl;	
	}



}

