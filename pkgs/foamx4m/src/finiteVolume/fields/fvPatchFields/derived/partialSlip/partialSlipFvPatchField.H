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
	Foam::partialSlipFvPatchField

Description
	Foam::partialSlipFvPatchField

SourceFiles
	partialSlipFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef partialSlipFvPatchField_H
#define partialSlipFvPatchField_H

#include "transformFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class partialSlipFvPatchField
:
	public transformFvPatchField<Type>
{
	// Private data

		//- Fraction (0-1) of value used for boundary condition
		scalarField valueFraction_;


public:

	//- Runtime type information
	TypeName("partialSlip");


	// Constructors

		//- Construct from patch and internal field
		partialSlipFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		partialSlipFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given partialSlipFvPatchField onto a new patch
		partialSlipFvPatchField
		(
			const partialSlipFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		partialSlipFvPatchField
		(
			const partialSlipFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new partialSlipFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		partialSlipFvPatchField
		(
			const partialSlipFvPatchField<Type>&,
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
				new partialSlipFvPatchField<Type>(*this, iF)
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


		// Return defining fields

			virtual scalarField& valueFraction()
			{
				return valueFraction_;
			}

			virtual const scalarField& valueFraction() const
			{
				return valueFraction_;
			}


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


		//- Write
		virtual void write(Ostream&) const;


	// Member operators

		virtual void operator=(const UList<Type>&) {}

		virtual void operator=(const fvPatchField<Type>&) {}
		virtual void operator+=(const fvPatchField<Type>&) {}
		virtual void operator-=(const fvPatchField<Type>&) {}
		virtual void operator*=(const fvPatchField<scalar>&) {}
		virtual void operator/=(const fvPatchField<scalar>&) {}

		virtual void operator+=(const Field<Type>&) {}
		virtual void operator-=(const Field<Type>&) {}

		virtual void operator*=(const Field<scalar>&) {}
		virtual void operator/=(const Field<scalar>&) {}

		virtual void operator=(const Type&) {}
		virtual void operator+=(const Type&) {}
		virtual void operator-=(const Type&) {}
		virtual void operator*=(const scalar) {}
		virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "partialSlipFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
