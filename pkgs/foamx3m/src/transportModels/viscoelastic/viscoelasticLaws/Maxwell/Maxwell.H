

#ifndef Maxwell_H
#define Maxwell_H

#include "viscoelasticLaw.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


class Maxwell
:
	public viscoelasticLaw
{
	// Private data

		//- Transported viscoelastic stress
		volSymmTensorField tau_;


		// Model constants

			//- Density
			dimensionedScalar rho_;

			//- Solvent viscosity
			dimensionedScalar etaS_;

			//- Zero shear rate polymer viscosity
			dimensionedScalar etaP_;

			//- Relaxation time
			dimensionedScalar lambda_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		Maxwell(const Maxwell&);

		//- Disallow default bitwise assignment
		void operator=(const Maxwell&);


public:

	//- Runtime type information
	TypeName("Maxwell");

	// Constructors

		//- Construct from components
		Maxwell
		(
			const word& name,
			const volVectorField& U,
			const surfaceScalarField& phi,
			const dictionary& dict
		);


	// Destructor

		virtual ~Maxwell()
		{}


	// Member Functions

		//- Return the viscoelastic stress tensor
		virtual tmp<volSymmTensorField> tau() const
		{
			return tau_;
		}

		//- Return the coupling term for the momentum equation
		virtual tmp<fvVectorMatrix> divTau(volVectorField& U) const;

		//- Correct the viscoelastic stress
		virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
