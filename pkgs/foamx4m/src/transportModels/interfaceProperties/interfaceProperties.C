/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.1
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
	This file is part of foam-extend.

	foam-extend is free software: you can redistribute it and/or modify it
	under the terms of the GNU General Public License as published by the
	Free Software Foundation, either version 3 of the License, or (at your
	option) any later version.

	foam-extend is distributed in the hope that it will be useful, but
	WITHOUT ANY WARRANTY; without even the implied warranty of
	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
	General Public License for more details.

	You should have received a copy of the GNU General Public License
	along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"

// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
	Foam::mathematicalConstant::pi/180.0;


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

// Correction for the boundary condition on the unit normal nHat on
// walls to produce the correct contact angle.

// The dynamic contact angle is calculated from the component of the
// velocity on the direction of the interface, parallel to the wall.

void Foam::interfaceProperties::correctContactAngle
(
	surfaceVectorField::Boundary& nHatb
) const
{
	const fvMesh& mesh = alpha1_.mesh();
	const volScalarField::Boundary& gbf = alpha1_.boundaryField();

	const fvBoundaryMesh& boundary = mesh.boundary();

	forAll(boundary, patchi)
	{
		if (isA<alphaContactAngleFvPatchScalarField>(gbf[patchi]))
		{
			const alphaContactAngleFvPatchScalarField& gcap =
				refCast<const alphaContactAngleFvPatchScalarField>
				(gbf[patchi]);

			fvsPatchVectorField& nHatp = nHatb[patchi];
			scalarField theta =
				convertToRad*gcap.theta(U_.boundaryField()[patchi], nHatp);

			vectorField nf = boundary[patchi].nf();

			// Reset nHatp to correspond to the contact angle

			scalarField a12 = nHatp & nf;

			scalarField b1 = cos(theta);

			scalarField b2(nHatp.size());

			forAll(b2, facei)
			{
				b2[facei] = cos(acos(a12[facei]) - theta[facei]);
			}

			scalarField det = 1.0 - a12*a12;

			scalarField a = (b1 - a12*b2)/det;
			scalarField b = (b2 - a12*b1)/det;

			nHatp = a*nf + b*nHatp;

			nHatp /= (mag(nHatp) + deltaN_.value());
		}
	}
}


void Foam::interfaceProperties::calculateK()
{
	const fvMesh& mesh = alpha1_.mesh();
	const surfaceVectorField& Sf = mesh.Sf();

	// Cell gradient of alpha
	const tmp<volVectorField> gradAlpha = fvc::grad(alpha1_);

	// Interpolated face-gradient of alpha
	surfaceVectorField gradAlphaf = fvc::interpolate(gradAlpha);
	//gradAlphaf -=
	//    (mesh.Sf()/mesh.magSf())
	//   *(fvc::snGrad(alpha1_) - (mesh.Sf() & gradAlphaf)/mesh.magSf());

	// Face unit interface normal
	surfaceVectorField nHatfv = gradAlphaf/(mag(gradAlphaf) + deltaN_);
	correctContactAngle(nHatfv.boundaryField());

	// Face unit interface normal flux
	nHatf_ = nHatfv & Sf;

	// Simple expression for curvature
	K_ = -fvc::div(nHatf_);

	// Complex expression for curvature.
	// Correction is formally zero but numerically non-zero.
	/*
	volVectorField nHat = gradAlpha/(mag(gradAlpha) + deltaN_);
	forAll(nHat.boundaryField(), patchi)
	{
		nHat.boundaryField()[patchi] = nHatfv.boundaryField()[patchi];
	}

	K_ = -fvc::div(nHatf_) + (nHat & fvc::grad(nHatfv) & nHat);
	*/
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
	const volScalarField& alpha1,
	const volVectorField& U,
	const IOdictionary& dict
)
:
	transportPropertiesDict_(dict),
	cAlpha_
	(
		readScalar
		(
			alpha1.mesh().solutionDict().subDict("PISO").lookup("cAlpha")
		)
	),
	sigma_(dict.lookup("sigma")),

	deltaN_
	(
		"deltaN",
		1e-8/cbrt(average(alpha1.mesh().V()))
	),

	alpha1_(alpha1),
	U_(U),

	nHatf_
	(
		IOobject
		(
			"nHatf",
			alpha1_.time().timeName(),
			alpha1_.mesh()
		),
		alpha1_.mesh(),
		dimensionedScalar("nHatf", dimArea, 0.0)
	),

	K_
	(
		IOobject
		(
			"K",
			alpha1_.time().timeName(),
			alpha1_.mesh()
		),
		alpha1_.mesh(),
		dimensionedScalar("K", dimless/dimLength, 0.0)
	)
{
	calculateK();
}


// ************************************************************************* //
