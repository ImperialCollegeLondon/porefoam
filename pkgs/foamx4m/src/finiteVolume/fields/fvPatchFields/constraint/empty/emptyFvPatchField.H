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
	Foam::emptyFvPatchField

Description
	Foam::emptyFvPatchField

SourceFiles
	emptyFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef emptyFvPatchField_H
#define emptyFvPatchField_H

#include "fvPatchField.H"
#include "emptyFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class emptyFvPatchField
:
	public fvPatchField<Type>
{

public:

	//- Runtime type information
	TypeName(emptyFvPatch::typeName_());


	// Constructors

		//- Construct from patch and internal field
		emptyFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		emptyFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given emptyFvPatchField onto a new patch
		emptyFvPatchField
		(
			const emptyFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		emptyFvPatchField
		(
			const emptyFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new emptyFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		emptyFvPatchField
		(
			const emptyFvPatchField<Type>&,
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
				new emptyFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Mapping functions

			//- Map (and resize as needed) from self given a mapping object
			virtual void autoMap
			(
				const fvPatchFieldMapper&
			)
			{}

			//- Reverse map the given fvPatchField onto this fvPatchField
			virtual void rmap
			(
				const fvPatchField<Type>&,
				const labelList&
			)
			{}


		// Evaluation functions

			//- Update the coefficients associated with the patch field
			//  This only checks to see the case is actually 1D or 2D
			//  for which this boundary condition is valid
			void updateCoeffs();


			//- Return the matrix diagonal coefficients corresponding to the
			//  evaluation of the value of this patchField with given weights
			virtual tmp<Field<Type> > valueInternalCoeffs
			(
				const tmp<scalarField>&
			) const
			{
				return tmp<Field<Type> >(new Field<Type>(0));
			}

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the value of this patchField with given weights
			virtual tmp<Field<Type> > valueBoundaryCoeffs
			(
				const tmp<scalarField>&
			) const
			{
				return tmp<Field<Type> >(new Field<Type>(0));
			}

			//- Return the matrix diagonal coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			tmp<Field<Type> > gradientInternalCoeffs() const
			{
				return tmp<Field<Type> >(new Field<Type>(0));
			}

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			tmp<Field<Type> > gradientBoundaryCoeffs() const
			{
				return tmp<Field<Type> >(new Field<Type>(0));
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "emptyFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
