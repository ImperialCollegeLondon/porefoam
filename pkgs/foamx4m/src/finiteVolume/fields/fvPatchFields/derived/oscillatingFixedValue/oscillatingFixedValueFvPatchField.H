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
	Foam::oscillatingFixedValueFvPatchField

Description
	Foam::oscillatingFixedValueFvPatchField

SourceFiles
	oscillatingFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingFixedValueFvPatchField_H
#define oscillatingFixedValueFvPatchField_H

#include "Random.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class oscillatingFixedValueFvPatchField
:
	public fixedValueFvPatchField<Type>
{
	// Private data

		//- Reference value
		Field<Type> refValue_;

		//- Amplitude
		scalar amplitude_;

		//- Frequency
		scalar frequency_;

		//- Current time index
		label curTimeIndex_;


	// Private member functions

		//- Return current scale
		scalar currentScale() const;


public:

	//- Runtime type information
	TypeName("oscillatingFixedValue");


	// Constructors

		//- Construct from patch and internal field
		oscillatingFixedValueFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		oscillatingFixedValueFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given oscillatingFixedValueFvPatchField
		//  onto a new patch
		oscillatingFixedValueFvPatchField
		(
			const oscillatingFixedValueFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		oscillatingFixedValueFvPatchField
		(
			const oscillatingFixedValueFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new oscillatingFixedValueFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		oscillatingFixedValueFvPatchField
		(
			const oscillatingFixedValueFvPatchField<Type>&,
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
				new oscillatingFixedValueFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return the ref value
			const Field<Type>& refValue() const
			{
				return refValue_;
			}

			//- Return reference to the ref value to allow adjustment
			Field<Type>& refValue()
			{
				return refValue_;
			}

			//- Return amplitude
			scalar amplitude() const
			{
				return amplitude_;
			}

			scalar& amplitude()
			{
				return amplitude_;
			}

			//- Return frequency
			scalar frequency() const
			{
				return frequency_;
			}

			scalar& frequency()
			{
				return frequency_;
			}


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

			//- Update the coefficients associated with the patch field
			virtual void updateCoeffs();


		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "oscillatingFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
