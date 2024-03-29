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

Class
	Foam::fv::SLTSDdtScheme

Description
	Stabilised local time-step first-order Euler implicit/explicit ddt.
	The time-step is adjusted locally so that an advective equations remains
	diagonally dominant.

	This scheme should only be used for steady-state computations
	using transient codes where local time-stepping is preferably to
	under-relaxation for transport consistency reasons.

	See also CoEulerDdtScheme.

SourceFiles
	SLTSDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef SLTSDdtScheme_H
#define SLTSDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class SLTSDdtScheme
:
	public fv::ddtScheme<Type>
{
	// Private Data

		//- Name of the flux field used to calculate the local time-step
		word phiName_;

		//- Name of the density field used to obtain the volumetric flux
		//  from the mass flux if required
		word rhoName_;

		//- Under-relaxation factor
		scalar alpha_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		SLTSDdtScheme(const SLTSDdtScheme&);

		//- Disallow default bitwise assignment
		void operator=(const SLTSDdtScheme&);

		//- Calculate a relaxed diagonal from the given flux field
		void relaxedDiag(scalarField& rD, const surfaceScalarField& phi) const;

		//- Return the reciprocal of the stabilised local time-step
		tmp<volScalarField> SLrDeltaT() const;


public:

	//- Runtime type information
	TypeName("SLTS");


	// Constructors

		//- Construct from mesh and Istream
		SLTSDdtScheme(const fvMesh& mesh, Istream& is)
		:
			ddtScheme<Type>(mesh, is),
			phiName_(is),
			rhoName_(is),
			alpha_(readScalar(is))
		{}


	// Member Functions

		//- Return mesh reference
		const fvMesh& mesh() const
		{
			return fv::ddtScheme<Type>::mesh();
		}

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const dimensioned<Type>&
		);

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const dimensionedScalar&,
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const volScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<fvMatrix<Type> > fvmDdt
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<fvMatrix<Type> > fvmDdt
		(
			const dimensionedScalar&,
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		tmp<fvMatrix<Type> > fvmDdt
		(
			const volScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		);

		typedef typename ddtScheme<Type>::fluxFieldType fluxFieldType;

		tmp<fluxFieldType> fvcDdtPhiCorr
		(
			const volScalarField& rA,
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const fluxFieldType& phi
		);

		tmp<fluxFieldType> fvcDdtPhiCorr
		(
			const volScalarField& rA,
			const volScalarField& rho,
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const fluxFieldType& phi
		);


		// Member functions for the new time consistent formulation
		tmp<fluxFieldType> fvcDdtConsistentPhiCorr
		(
			const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const surfaceScalarField& rAUf
		);


		tmp<surfaceScalarField> meshPhi
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		);
};


template<>
tmp<surfaceScalarField> SLTSDdtScheme<scalar>::fvcDdtPhiCorr
(
	const volScalarField& rA,
	const volScalarField& U,
	const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> SLTSDdtScheme<scalar>::fvcDdtPhiCorr
(
	const volScalarField& rA,
	const volScalarField& rho,
	const volScalarField& U,
	const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> SLTSDdtScheme<scalar>::fvcDdtConsistentPhiCorr
(
	const surfaceScalarField& faceU,
	const volScalarField& U,
	const surfaceScalarField& rAUf
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "SLTSDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
