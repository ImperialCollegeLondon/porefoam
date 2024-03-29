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
	Foam::fv::CrankNicolsonDdtScheme

Description
	Second-oder CrankNicolson implicit ddt using the current and
	previous time-step fields as well as the previous time-step ddt.

SourceFiles
	CrankNicolsonDdtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef CrankNicolsonDdtScheme_H
#define CrankNicolsonDdtScheme_H

#include "ddtScheme.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class CrankNicolsonDdtScheme
:
	public fv::ddtScheme<Type>
{
	// Private Data

		//- Class to store the ddt0 fields on the objectRegistry for use in the
		//  next time-step.  The start-time index of the CN scheme is also
		//  stored to help handle the transition from Euler to CN
		template<class GeoField>
		class DDt0Field
		:
			public GeoField
		{
			label startTimeIndex_;

		public:

			//- Constructor from file for restart.
			DDt0Field
			(
				const IOobject& io,
				const fvMesh& mesh
			);

			//- Constructor from components, initisalised to zero with given
			//  dimensions.
			DDt0Field
			(
				const IOobject& io,
				const fvMesh& mesh,
				const dimensioned<typename GeoField::value_type>& dimType
			);

			//- Return the start-time index
			label startTimeIndex() const;

			//- Cast to the underlying GeoField
			GeoField& operator()();

			//- Assignment to a GeoField
			void operator=(const GeoField& gf);
		};


		//- Off-centering coefficient, 1 -> CN, less than one blends with EI
		scalar ocCoeff_;


	// Private Member Functions

		//- Disallow default bitwise copy construct
		CrankNicolsonDdtScheme(const CrankNicolsonDdtScheme&);

		//- Disallow default bitwise assignment
		void operator=(const CrankNicolsonDdtScheme&);

		template<class GeoField>
		DDt0Field<GeoField>& ddt0_
		(
			const word& name,
			const dimensionSet& dims
		);

		//- Check if the ddt0 needs to be evaluated for this time-step
		template<class GeoField>
		bool evaluate(const DDt0Field<GeoField>& ddt0) const;

		//- Return the coefficient for Euler scheme for the first time-step
		//  for and CN thereafter
		template<class GeoField>
		scalar coef_(const DDt0Field<GeoField>&) const;

		//- Return the old time-step coefficient for Euler scheme for the
		//  second time-step and for CN thereafter
		template<class GeoField>
		scalar coef0_(const DDt0Field<GeoField>&) const;

		//- Return the reciprocal time-step coefficient for Euler for the
		//  first time-step and CN thereafter
		template<class GeoField>
		dimensionedScalar rDtCoef_(const DDt0Field<GeoField>&) const;

		//- Return the reciprocal old time-step coefficient for Euler for the
		//  second time-step and CN thereafter
		template<class GeoField>
		dimensionedScalar rDtCoef0_(const DDt0Field<GeoField>&) const;

		//- Return ddt0 multiplied by the off-centreing coefficient
		template<class GeoField>
		tmp<GeoField> offCentre_(const GeoField& ddt0) const;


public:

	//- Runtime type information
	TypeName("CrankNicolson");


	// Constructors

		//- Construct from mesh
		CrankNicolsonDdtScheme(const fvMesh& mesh)
		:
			ddtScheme<Type>(mesh),
			ocCoeff_(1.0)
		{}

		//- Construct from mesh and Istream
		CrankNicolsonDdtScheme(const fvMesh& mesh, Istream& is)
		:
			ddtScheme<Type>(mesh, is),
			ocCoeff_(readScalar(is))
		{
			if (ocCoeff_ < 0 || ocCoeff_ > 1)
			{
				FatalIOErrorIn
				(
					"CrankNicolsonDdtScheme(const fvMesh& mesh, Istream& is)",
					is
				)   << "coefficient = " << ocCoeff_
					<< " should be >= 0 and <= 1"
					<< exit(FatalIOError);
			}
		}


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
		virtual tmp<fluxFieldType> fvcDdtConsistentPhiCorr
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
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtPhiCorr
(
	const volScalarField& rA,
	const volScalarField& U,
	const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtPhiCorr
(
	const volScalarField& rA,
	const volScalarField& rho,
	const volScalarField& U,
	const surfaceScalarField& phi
);


template<>
tmp<surfaceScalarField> CrankNicolsonDdtScheme<scalar>::fvcDdtConsistentPhiCorr
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
#	include "CrankNicolsonDdtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
