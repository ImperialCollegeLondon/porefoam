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
	Foam::fv::ddtScheme

Description
	Abstract base class for ddt schemes.

SourceFiles
	ddtScheme.C

\*---------------------------------------------------------------------------*/

#ifndef ddtScheme_H
#define ddtScheme_H

#include "tmp.H"
#include "dimensionedType.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

template<class Type>
class fvMatrix;

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class ddtScheme
:
	public refCount
{

protected:

	// Protected data

		const fvMesh& mesh_;


	// Private Member Functions

		//- Disallow copy construct
		ddtScheme(const ddtScheme&);

		//- Disallow default bitwise assignment
		void operator=(const ddtScheme&);


public:

	//- Runtime type information
	virtual const word& type() const = 0;


	// Declare run-time constructor selection tables

		declareRunTimeSelectionTable
		(
			tmp,
			ddtScheme,
			Istream,
			(const fvMesh& mesh, Istream& schemeData),
			(mesh, schemeData)
		);


	// Constructors

		//- Construct from mesh
		ddtScheme(const fvMesh& mesh)
		:
			mesh_(mesh)
		{}

		//- Construct from mesh and Istream
		ddtScheme(const fvMesh& mesh, Istream&)
		:
			mesh_(mesh)
		{}


	// Selectors

		//- Return a pointer to a new ddtScheme created on freestore
		static tmp<ddtScheme<Type> > New
		(
			const fvMesh& mesh,
			Istream& schemeData
		);


	// Destructor

		virtual ~ddtScheme();


	// Member Functions

		//- Return mesh reference
		const fvMesh& mesh() const
		{
			return mesh_;
		}

		virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const dimensioned<Type>&
		) = 0;

		virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;

		virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const dimensionedScalar&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;

		virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDdt
		(
			const volScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;

		virtual tmp<fvMatrix<Type> > fvmDdt
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;

		virtual tmp<fvMatrix<Type> > fvmDdt
		(
			const dimensionedScalar&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;

		virtual tmp<fvMatrix<Type> > fvmDdt
		(
			const volScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;


		typedef GeometricField
		<
			typename flux<Type>::type,
			fvsPatchField,
			surfaceMesh
		> fluxFieldType;

		tmp<surfaceScalarField> fvcDdtPhiCoeff
		(
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const fluxFieldType& phi,
			const fluxFieldType& phiCorr
		);

		tmp<surfaceScalarField> fvcDdtPhiCoeff
		(
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const fluxFieldType& phi
		);

		virtual tmp<fluxFieldType> fvcDdtPhiCorr
		(
			const volScalarField& rA,
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const fluxFieldType& phi
		) = 0;

		tmp<surfaceScalarField> fvcDdtPhiCoeff
		(
			const volScalarField& rho,
			const GeometricField<Type, fvPatchField, volMesh>& rhoU,
			const fluxFieldType& phi
		);

		virtual tmp<fluxFieldType> fvcDdtPhiCorr
		(
			const volScalarField& rA,
			const volScalarField& rho,
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const fluxFieldType& phi
		) = 0;


		// Member functions for the new time consistent formulation
		virtual tmp<fluxFieldType> fvcDdtConsistentPhiCorr
		(
			const GeometricField<Type, fvsPatchField, surfaceMesh>& faceU,
			const GeometricField<Type, fvPatchField, volMesh>& U,
			const surfaceScalarField& rAUf
		) = 0;


		virtual tmp<surfaceScalarField> meshPhi
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvDdtTypeScheme(SS, Type)                                          \
					                                                           \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                              \
					                                                           \
ddtScheme<Type>::addIstreamConstructorToTable<SS<Type> >                       \
	add##SS##Type##IstreamConstructorToTable_;


#define makeFvDdtScheme(SS)                                                    \
					                                                           \
makeFvDdtTypeScheme(SS, scalar)                                                \
makeFvDdtTypeScheme(SS, vector)                                                \
makeFvDdtTypeScheme(SS, sphericalTensor)                                       \
makeFvDdtTypeScheme(SS, symmTensor)                                            \
makeFvDdtTypeScheme(SS, symmTensor4thOrder)                                    \
makeFvDdtTypeScheme(SS, diagTensor)                                            \
makeFvDdtTypeScheme(SS, tensor)                                                \
					                                                           \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtPhiCorr                              \
(                                                                              \
	const volScalarField& rA,                                                  \
	const volScalarField& U,                                                   \
	const surfaceScalarField& phi                                              \
)                                                                              \
{                                                                              \
	notImplemented(#SS"<scalar>::fvcDdtPhiCorr");                              \
	return surfaceScalarField::null();                                         \
}                                                                              \
					                                                           \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtPhiCorr                              \
(                                                                              \
	const volScalarField& rA,                                                  \
	const volScalarField& rho,                                                 \
	const volScalarField& U,                                                   \
	const surfaceScalarField& phi                                              \
)                                                                              \
{                                                                              \
	notImplemented(#SS"<scalar>::fvcDdtPhiCorr");                              \
	return surfaceScalarField::null();                                         \
}                                                                              \
					                                                           \
template<>                                                                     \
tmp<surfaceScalarField> SS<scalar>::fvcDdtConsistentPhiCorr                    \
(                                                                              \
	const surfaceScalarField& faceU,                                           \
	const volScalarField& U,                                                   \
	const surfaceScalarField& rAUf                                             \
)                                                                              \
{                                                                              \
	notImplemented(#SS"<scalar>::fvcDdtConsistentPhiCorr");                    \
	return surfaceScalarField::null();                                         \
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "ddtScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
