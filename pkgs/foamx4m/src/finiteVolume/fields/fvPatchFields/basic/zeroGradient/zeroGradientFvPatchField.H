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
	Foam::zeroGradientFvPatchField

Description
	Foam::zeroGradientFvPatchField

SourceFiles
	zeroGradientFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroGradientFvPatchField_H
#define zeroGradientFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class zeroGradientFvPatchField
:
	public fvPatchField<Type>
{

public:

	//- Runtime type information
	TypeName("zeroGradient");


	// Constructors

		//- Construct from patch and internal field
		zeroGradientFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		zeroGradientFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping the given zeroGradientFvPatchField<Type>
		//  onto a new patch
		zeroGradientFvPatchField
		(
			const zeroGradientFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		zeroGradientFvPatchField
		(
			const zeroGradientFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new zeroGradientFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		zeroGradientFvPatchField
		(
			const zeroGradientFvPatchField<Type>&,
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
				new zeroGradientFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Evaluation functions

			//- Return gradient at boundary
			virtual tmp<Field<Type> > snGrad() const
			{
				return tmp<Field<Type> >
				(
					new Field<Type>(this->size(), pTraits<Type>::zero)
				);
			}

			//- Evaluate the patch field
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
			virtual tmp<Field<Type> > gradientInternalCoeffs() const;

			//- Return the matrix source coefficients corresponding to the
			//  evaluation of the gradient of this patchField
			virtual tmp<Field<Type> > gradientBoundaryCoeffs() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "zeroGradientFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
