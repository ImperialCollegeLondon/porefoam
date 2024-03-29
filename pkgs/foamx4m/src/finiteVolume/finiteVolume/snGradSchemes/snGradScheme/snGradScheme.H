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
	Foam::fv::snGradScheme

Description
	Abstract base class for snGrad schemes.

SourceFiles
	snGradScheme.C

\*---------------------------------------------------------------------------*/

#ifndef snGradScheme_H
#define snGradScheme_H

#include "tmp.H"
#include "volFieldsFwd.H"
#include "surfaceFieldsFwd.H"
#include "typeInfo.H"
#include "runTimeSelectionTables.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class fvMesh;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace fv
{


template<class Type>
class snGradScheme
:
	public refCount
{
	// Private data

		//- Hold reference to mesh
		const fvMesh& mesh_;


	// Private Member Functions

		//- Disallow default bitwise assignment
		void operator=(const snGradScheme&);


public:

	//- Runtime type information
	virtual const word& type() const = 0;


	// Declare run-time constructor selection tables

		declareRunTimeSelectionTable
		(
			tmp,
			snGradScheme,
			Mesh,
			(const fvMesh& mesh, Istream& schemeData),
			(mesh, schemeData)
		);


	// Constructors

		//- Construct from mesh
		snGradScheme(const fvMesh& mesh)
		:
			mesh_(mesh)
		{}


	// Selectors

		//- Return new tmp interpolation scheme
		static tmp<snGradScheme<Type> > New
		(
			const fvMesh& mesh,
			Istream& schemeData
		);


	// Destructor

		virtual ~snGradScheme();


	// Member Functions

		//- Return mesh reference
		const fvMesh& mesh() const
		{
			return mesh_;
		}


		//- Return the snGrad of the given cell field with the given deltaCoeffs
		static tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		snGrad
		(
			const GeometricField<Type, fvPatchField, volMesh>&,
			const tmp<surfaceScalarField>&,
			const word& snGradName = "snGrad"
		);

		//- Return the interpolation weighting factors for the given field
		virtual tmp<surfaceScalarField> deltaCoeffs
		(
			const GeometricField<Type, fvPatchField, volMesh>&
		) const = 0;

		//- Return true if this scheme uses an explicit correction
		virtual bool corrected() const
		{
			return false;
		}

		//- Return the explicit correction to the snGrad
		//  for the given field
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		correction(const GeometricField<Type, fvPatchField, volMesh>&) const
		{
			return tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >(nullptr);
		}

		//- Return the snGrad of the given cell field
		//  with explicit correction
		virtual tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		snGrad(const GeometricField<Type, fvPatchField, volMesh>&) const;

		//- Return the snGrad of the given tmp cell field
		//  with explicit correction
		tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
		snGrad
		(
			const tmp<GeometricField<Type, fvPatchField, volMesh> >&
		) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Add the patch constructor functions to the hash tables

#define makeSnGradTypeScheme(SS, Type)                                        \
					                                                          \
defineNamedTemplateTypeNameAndDebug(SS<Type>, 0);                             \
					                                                          \
snGradScheme<Type>::addMeshConstructorToTable<SS<Type> >                      \
	add##SS##Type##MeshConstructorToTable_;

#define makeSnGradScheme(SS)                                                  \
					                                                          \
makeSnGradTypeScheme(SS, scalar)                                              \
makeSnGradTypeScheme(SS, vector)                                              \
makeSnGradTypeScheme(SS, sphericalTensor)                                     \
makeSnGradTypeScheme(SS, symmTensor)                                          \
makeSnGradTypeScheme(SS, symmTensor4thOrder)                                  \
makeSnGradTypeScheme(SS, diagTensor)                                          \
makeSnGradTypeScheme(SS, tensor)


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "snGradScheme.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
