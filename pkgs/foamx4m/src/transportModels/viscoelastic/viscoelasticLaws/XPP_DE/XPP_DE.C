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

#include "XPP_DE.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(XPP_DE, 0);
	addToRunTimeSelectionTable(viscoelasticLaw, XPP_DE, dictionary);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::XPP_DE::XPP_DE
(
	const word& name,
	const volVectorField& U,
	const surfaceScalarField& phi,
	const dictionary& dict
)
:
	viscoelasticLaw(name, U, phi),
	S_
	(
		IOobject
		(
			"S" + name,
			U.time().timeName(),
			U.mesh(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		U.mesh()
	),
	Lambda_
	(
		IOobject
		(
			"Lambda" + name,
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
	alpha_(dict.lookup("alpha")),
	lambdaOb_(dict.lookup("lambdaOb")),
	lambdaOs_(dict.lookup("lambdaOs")),
	q_(dict.lookup("q"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::XPP_DE::divTau(volVectorField& U) const
{
	 dimensionedScalar etaPEff = etaP_;

	return
	(
		fvc::div(tau_/rho_, "div(tau)")
	  - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
	  + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
	);

}


void Foam::XPP_DE::correct()
{
	// Velocity gradient tensor
	const tmp<volTensorField> tL = fvc::grad(U());
	const volTensorField& L = tL();

	// Convected derivate term
	volTensorField C = S_ & L;

	// Twice the rate of deformation tensor
	volSymmTensorField twoD = twoSymm(L);

	// Evolution of orientation
	tmp<fvSymmTensorMatrix> SEqn
	(
		fvm::ddt(S_)
	  + fvm::div(phi(), S_)
	 ==
		twoSymm(C)
	  - fvm::Sp((twoD && S_) , S_)
	  - fvm::Sp
		(
			1/lambdaOb_/Foam::sqr(Lambda_)*
			(1 - alpha_ - 3*alpha_*Foam::pow(Lambda_, 4)*tr(S_ & S_)),
			S_
		)
		- 1/lambdaOb_/Foam::sqr(Lambda_)*
		(3*alpha_*Foam::pow(Lambda_, 4)*symm(S_ & S_) - (1 - alpha_)/3*I_)
	);

	SEqn().relax();
	solve(SEqn);


	 // Evolution of the backbone stretch
	fvScalarMatrix LambdaEqn
	(
		fvm::ddt(Lambda_)
	  + fvm::div(phi(), Lambda_)
	 ==
		fvm::Sp((twoD && S_)/2 , Lambda_)
	  - fvm::Sp(Foam::exp( 2/q_*(Lambda_ - 1))/lambdaOs_ , Lambda_)
	  + Foam::exp(2/q_*(Lambda_ - 1))/lambdaOs_
	);

	LambdaEqn.relax();
	LambdaEqn.solve();

	// Viscoelastic stress
	tau_ = etaP_/lambdaOb_*(3*Foam::sqr(Lambda_)*S_ - I_);
}


// ************************************************************************* //
