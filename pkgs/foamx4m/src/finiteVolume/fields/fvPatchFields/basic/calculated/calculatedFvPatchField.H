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
	Foam::calculatedFvPatchField

Description
	Foam::calculatedFvPatchField

SourceFiles
	calculatedFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef calculatedFvPatchField_H
#define calculatedFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class calculatedFvPatchField
:
	public fvPatchField<Type>
{

public:

	//- Runtime type information
	TypeName("calculated");


	// Constructors

		//- Construct from patch and internal field
		calculatedFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		calculatedFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&,
			const bool valueRequired=false
		);

		//- Construct by mapping given patchField<Type> onto a new patch
		calculatedFvPatchField
		(
			const calculatedFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		calculatedFvPatchField
		(
			const calculatedFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new calculatedFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		calculatedFvPatchField
		(
			const calculatedFvPatchField<Type>&,
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
				new calculatedFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Access

			//- Return true if this patch field fixes a value.
			//  Needed to check if a level has to be specified while solving
			//  Poissons equations.
			virtual bool fixesValue() const
			{
				return true;
			}


		// Evaluation functions

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
#	include "calculatedFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
