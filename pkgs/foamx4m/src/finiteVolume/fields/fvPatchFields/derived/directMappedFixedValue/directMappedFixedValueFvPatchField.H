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
	Foam::directMappedFixedValueFvPatchField

Description
	Recycles the value at a set of internal faces back to *this.

SourceFiles
	directMappedFixedValueFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef directMappedFixedValueFvPatchField_H
#define directMappedFixedValueFvPatchField_H

#include "fixedValueFvPatchFields.H"
#include "directMappedFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class directMappedFixedValueFvPatchField
:
	public fixedValueFvPatchField<Type>
{
	// Private data

		//- If true adjust the mapped field to maintain average value average_
		bool setAverage_;

		//- Average value the mapped field is adjusted to maintain if
		//  setAverage_ is set true
		Type average_;

		//- Current time index
		label curTimeIndex_;

		//- New values from interpolation
		Field<Type> newValues_;

	// Private member functions

		//- Map field data
		void mapField();


public:

	//- Runtime type information
	TypeName("directMapped");


	// Constructors

		//- Construct from patch and internal field
		directMappedFixedValueFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		directMappedFixedValueFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given directMappedFixedValueFvPatchField
		//  onto a new patch
		directMappedFixedValueFvPatchField
		(
			const directMappedFixedValueFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		directMappedFixedValueFvPatchField
		(
			const directMappedFixedValueFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new directMappedFixedValueFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		directMappedFixedValueFvPatchField
		(
			const directMappedFixedValueFvPatchField<Type>&,
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
				new directMappedFixedValueFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Evaluation functions

			//- Initialise the evaluation of the patch field
			virtual void initEvaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			);


			//- Evaluate the patch field, sets Updated to false
			virtual void evaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			);


		//- Write
		virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "directMappedFixedValueFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
