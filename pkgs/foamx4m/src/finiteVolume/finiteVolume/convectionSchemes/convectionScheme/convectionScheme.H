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
	Foam::fv::convectionScheme

Description
	Abstract base class for convection schemes.

SourceFiles
	convectionScheme.C

\*---------------------------------------------------------------------------*/

#ifndef convectionScheme_H
#define convectionScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"
#include "multivariateSurfaceInterpolationScheme.H"

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
class convectionScheme
:
	public refCount
{
	// Private data

		const fvMesh& mesh_;


public:

	//- Runtime type information
	virtual const word& type() const = 0;


	// Declare run-time constructor selection tables

		declareRunTimeSelectionTable
		(
			tmp,
			convectionScheme,
			Istream,
			(
				const fvMesh& mesh,
				const surfaceScalarField& faceFlux,
				Istream& schemeData
			),
			(mesh, faceFlux, schemeData)
		);

		declareRunTimeSelectionTable
		(
			tmp,
			convectionScheme,
			Multivariate,
			(
				const fvMesh& mesh,
				const typename multivariateSurfaceInterpolationScheme<Type>::
					fieldTable& fields,
				const surfaceScalarField& faceFlux,
				Istream& schemeData
			),
			(mesh, fields, faceFlux, schemeData)
		);


	// Constructors

		//- Copy construct
		convectionScheme(const convectionScheme&);

		//- Construct from mesh, flux and Istream
		convectionScheme
		(
			const fvMesh& mesh,
			const surfaceScalarField&
		)
		:
			mesh_(mesh)
		{}


	// Selectors

		//- Return a pointer to a new convectionScheme created on freestore
		static tmp<convectionScheme<Type> > New
		(
			const fvMesh& mesh,
			const surfaceScalarField& faceFlux,
			Istream& schemeData
		);


		//- Return a pointer to a new multivariate convectionScheme
		//  created on freestore
		static tmp<convectionScheme<Type> > New
		(
			const fvMesh& mesh,
			const typename multivariateSurfaceInterpolationScheme<Type>::
				fieldTable& fields,
			const surfaceScalarField& faceFlux,
			Istream& schemeData
		);


	// Destructor

		virtual ~convectionScheme();


	// Member Functions

		//- Return mesh reference
		const fvMesh& mesh() const
		{
			return mesh_;
		}

		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		interpolate
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const = 0;

		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > flux
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const = 0;

		virtual tmp<fvMatrix<Type> > fvmDiv
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const = 0;

		virtual tmp<GeometricField<Type, fvPatchField, volMesh> > fvcDiv
		(
			const surfaceScalarField&,
			const GeometricField<Type, fvPatchField, volMesh>&
		) const = 0;


	// Member operators

		void operator=(const convectionScheme<Type>&);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeFvConvectionTypeScheme(SS, Type)                                  \
					                                                          \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                             \
					                                                          \
convectionScheme<Type>::addIstreamConstructorToTable<SS<Type> >               \
	add##SS##Type##IstreamConstructorToTable_;

#define makeFvConvectionScheme(SS)                                            \
					                                                          \
makeFvConvectionTypeScheme(SS, scalar)                                        \
makeFvConvectionTypeScheme(SS, vector)                                        \
makeFvConvectionTypeScheme(SS, sphericalTensor)                               \
makeFvConvectionTypeScheme(SS, symmTensor)                                    \
makeFvConvectionTypeScheme(SS, tensor)


#define makeMultivariateFvConvectionTypeScheme(SS, Type)                      \
					                                                          \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                             \
					                                                          \
convectionScheme<Type>::addMultivariateConstructorToTable<SS<Type> >          \
	add##SS##Type##MultivariateConstructorToTable_;


#define makeMultivariateFvConvectionScheme(SS)                                \
					                                                          \
makeMultivariateFvConvectionTypeScheme(SS, scalar)                            \
makeMultivariateFvConvectionTypeScheme(SS, vector)                            \
makeMultivariateFvConvectionTypeScheme(SS, sphericalTensor)                   \
makeMultivariateFvConvectionTypeScheme(SS, symmTensor)                        \
makeMultivariateFvConvectionTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "convectionScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
