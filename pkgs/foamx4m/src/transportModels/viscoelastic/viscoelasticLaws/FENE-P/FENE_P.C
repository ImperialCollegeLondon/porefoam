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

#include "FENE_P.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(FENE_P, 0);
	addToRunTimeSelectionTable(viscoelasticLaw, FENE_P, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::FENE_P::FENE_P
(
	const word& name,
	const volVectorField& U,
	const surfaceScalarField& phi,
	const dictionary& dict
)
:
	viscoelasticLaw(name, U, phi),
	tau_
	(
		IOobject
		(
			"tau" + name,
			U.time().timeName(),
			U.mesh(),
			IOobject::MUST_READ,
			IOobject::AUTO_WRITE
		),
		U.mesh()
	),
	rho_(dict.lookup("rho")),
	etaS_(dict.lookup("etaS")),
	etaP_(dict.lookup("etaP")),
	L2_(dict.lookup("L2")),
	lambda_(dict.lookup("lambda"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::FENE_P::divTau(volVectorField& U) const
{
	dimensionedScalar etaPEff = etaP_;

	return
	(
		fvc::div(tau_/rho_, "div(tau)")
	  - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
	  + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
	);

}


void Foam::FENE_P::correct()
{
	// Velocity gradient tensor
	const tmp<volTensorField> tL = fvc::grad(U());
	const volTensorField& L = tL();

	// Convected derivate term
	volTensorField C = tau_ & L;

	// Twice the rate of deformation tensor
	volSymmTensorField twoD = twoSymm(L);

	// Stress transport equation
	fvSymmTensorMatrix tauEqn
	(
		fvm::ddt(tau_)
	  + fvm::div(phi(), tau_)
	 ==
		(1/lambda_/(1 - 3/L2_))*etaP_*twoD
	  + twoSymm(C)
	  - fvm::Sp
		(
			1/lambda_ + (3/lambda_/(1 - 3/L2_) + tr(tau_)/etaP_)/(L2_),
			tau_
		)
	);

	tauEqn.relax();
	tauEqn.solve();
}


// ************************************************************************* //
