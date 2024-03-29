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

#include "Leonov.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(Leonov, 0);
	addToRunTimeSelectionTable(viscoelasticLaw, Leonov, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Leonov::Leonov
(
	const word& name,
	const volVectorField& U,
	const surfaceScalarField& phi,
	const dictionary& dict
)
:
	viscoelasticLaw(name, U, phi),
	sigma_
	(
		IOobject
		(
			"sigma" + name,
			U.time().timeName(),
			U.mesh(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		U.mesh()
	),
	tau_
	(
		IOobject
		(
			"tau" + name,
			U.time().timeName(),
			U.mesh(),
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		U.mesh(),
		dimensionedSymmTensor
		(
			"zero",
			dimensionSet(1, -1, -2, 0, 0, 0, 0),
			symmTensor::zero
		)
	),
	I_
	(
		dimensionedSymmTensor
		(
			"I",
			dimensionSet(0, 0, 0, 0, 0, 0, 0),
			symmTensor
			(
				1, 0, 0,
				   1, 0,
					  1
			)
		)
	),
	rho_(dict.lookup("rho")),
	etaS_(dict.lookup("etaS")),
	etaP_(dict.lookup("etaP")),
	lambda_(dict.lookup("lambda"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::Leonov::divTau(volVectorField& U) const
{
	dimensionedScalar etaPEff = etaP_;

	return
	(
		fvc::div(tau_/rho_, "div(tau)")
	  - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
	  + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
	);
}


void Foam::Leonov::correct()
{
	// Velocity gradient tensor
	const tmp<volTensorField> tL = fvc::grad(U());
	const volTensorField& L = tL();

	// Convected derivate term
	volTensorField C = sigma_ & L;

	// Twice the rate of deformation tensor
	volSymmTensorField twoD = twoSymm(L);

	 // Stress transport equation
	fvSymmTensorMatrix sigmaEqn
	(
		fvm::ddt(sigma_)
	  + fvm::div(phi(), sigma_)
	 ==
		twoSymm(C)
	  - 1/etaP_/2*(symm(sigma_ & sigma_) - Foam::sqr(etaP_/lambda_)*I_)
	  + fvm::Sp
		(
			1/etaP_/6*
			(
				tr(sigma_)
			  - Foam::sqr(etaP_/lambda_) * tr(inv(sigma_))
			),
			sigma_
		)
	);


	sigmaEqn.relax();
	sigmaEqn.solve();

	// Viscoelastic stress
	tau_ = sigma_ - etaP_/lambda_*I_;
}


// ************************************************************************* //
