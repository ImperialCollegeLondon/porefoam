

#include "Maxwell.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
	defineTypeNameAndDebug(Maxwell, 0);
	addToRunTimeSelectionTable(viscoelasticLaw, Maxwell, dictionary);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::Maxwell::Maxwell
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
	lambda_(dict.lookup("lambda"))
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::fvVectorMatrix> Foam::Maxwell::divTau(volVectorField& U) const
{
	dimensionedScalar etaPEff = etaP_;

	return
	(
		fvc::div(tau_/rho_, "div(tau)")
	  - fvc::laplacian(etaPEff/rho_, U, "laplacian(etaPEff,U)")
	  + fvm::laplacian( (etaPEff + etaS_)/rho_, U, "laplacian(etaPEff+etaS,U)")
	);

}


void Foam::Maxwell::correct()
{
	// Velocity gradient tensor
	volTensorField L = fvc::grad(U());

	// Twice the rate of deformation tensor
	volSymmTensorField twoD = twoSymm(L);

	 // Stress transport equation
	fvSymmTensorMatrix tauEqn
	(
		fvm::ddt(tau_)
	 ==
		etaP_/lambda_*twoD
	  - fvm::Sp( 1/lambda_, tau_ )
	);

	tauEqn.relax();
	tauEqn.solve();
}


// ************************************************************************* //
