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
	Foam::basicSymmetryFvPatchField

Description
	A symmetry patch

SourceFiles
	basicSymmetryFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef basicSymmetryFvPatchField_H
#define basicSymmetryFvPatchField_H

#include "transformFvPatchFields.H"
#include "symmetryFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class basicSymmetryFvPatchField
:
	public transformFvPatchField<Type>
{
	// Private Data

		//- Skew correction
		Switch skewCorrected_;

		//- Second order accuracy correction
		Switch secondOrder_;


public:

	//- Runtime type information

		//- Construct from patch and internal field
		basicSymmetryFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		basicSymmetryFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given basicSymmetryFvPatchField onto a
		//  new patch
		basicSymmetryFvPatchField
		(
			const basicSymmetryFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		basicSymmetryFvPatchField
		(
			const basicSymmetryFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new basicSymmetryFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		basicSymmetryFvPatchField
		(
			const basicSymmetryFvPatchField<Type>&,
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
				new basicSymmetryFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Evaluation functions

			//- Return gradient at boundary
			virtual tmp<Field<Type> > snGrad() const;

			//- Evaluate the patch field
			//  Default argument needed to allow call in constructors
			//  HJ, 30/Jun/2009
			virtual void evaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			);

			//- Return face-gradient transform diagonal
			virtual tmp<Field<Type> > snGradTransformDiag() const;


		// I-O

			//- Write
			virtual void write(Ostream& os) const;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField> basicSymmetryFvPatchField<scalar>::snGrad() const;

template<>
void basicSymmetryFvPatchField<scalar>::evaluate
(
	const Pstream::commsTypes commsType
);

template<>
tmp<vectorField> basicSymmetryFvPatchField<vector>::snGrad() const;

template<>
void basicSymmetryFvPatchField<vector>::evaluate
(
	const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "basicSymmetryFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
