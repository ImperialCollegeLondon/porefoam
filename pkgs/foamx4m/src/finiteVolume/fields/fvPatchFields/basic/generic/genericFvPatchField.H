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
	Foam::genericFvPatchField

Description
	Foam::genericFvPatchField

SourceFiles
	genericFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef genericFvPatchField_H
#define genericFvPatchField_H

#include "calculatedFvPatchField.H"
#include "HashPtrTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class genericFvPatchField
:
	public calculatedFvPatchField<Type>
{
	// Private data

		word actualTypeName_;
		dictionary dict_;

		HashPtrTable<scalarField> scalarFields_;
		HashPtrTable<vectorField> vectorFields_;
		HashPtrTable<sphericalTensorField> sphericalTensorFields_;
		HashPtrTable<symmTensorField> symmTensorFields_;
		HashPtrTable<symmTensor4thOrderField> symmTensor4thOrderFields_;
		HashPtrTable<diagTensorField> diagTensorFields_;
		HashPtrTable<tensorField> tensorFields_;


public:

	//- Runtime type information
	TypeName("generic");


	// Constructors

		//- Construct from patch and internal field
		genericFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		genericFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given patchField<Type> onto a new patch
		genericFvPatchField
		(
			const genericFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		genericFvPatchField
		(
			const genericFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new genericFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		genericFvPatchField
		(
			const genericFvPatchField<Type>&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchField<Type> > clone
		(
			const DimensionedField<Type, volMesh>& iF
		) const
		{
			return tmp<fvPatchField<Type> >
			(
				new genericFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Mapping functions

			//- Map (and resize as needed) from self given a mapping object
			virtual void autoMap
			(
				const fvPatchFieldMapper&
			);

			//- Reverse map the given fvPatchField onto this fvPatchField
			virtual void rmap
			(
				const fvPatchField<Type>&,
				const labelList&
			);


		// Evaluation functions

			//- Evaluate the patch field, sets Updated to false
			virtual void evaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			);


			//- Return the matrix diagonal coefficients corresponding to the
			//  evaluation of the value of this patchField with given weights
			virtual tmp<Field<Type> > valueInternalCoeffs
			(
				const tmp<scalarField>&
			) const;

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the value of this patchField with given weights
			virtual tmp<Field<Type> > valueBoundaryCoeffs
			(
				const tmp<scalarField>&
			) const;

			//- Return the matrix diagonal coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			tmp<Field<Type> > gradientInternalCoeffs() const;

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			tmp<Field<Type> > gradientBoundaryCoeffs() const;


		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "genericFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
