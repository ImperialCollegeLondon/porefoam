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
	Foam::transformFvPatchField

Description
	Foam::transformFvPatchField

SourceFiles
	transformFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef transformFvPatchField_H
#define transformFvPatchField_H

#include "fvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class transformFvPatchField
:
	public fvPatchField<Type>
{

public:

	//- Runtime type information
	TypeName("transform");


	// Constructors

		//- Construct from patch and internal field
		transformFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		transformFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping the given transformFvPatchField<Type>
		//  onto a new patch
		transformFvPatchField
		(
			const transformFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		transformFvPatchField
		(
			const transformFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const = 0;

		//- Construct as copy setting internal field reference
		transformFvPatchField
		(
			const transformFvPatchField<Type>&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchField<Type> > clone
		(
			const DimensionedField<Type, volMesh>&
		) const = 0;


	// Member functions

		// Evaluation functions

			//- Return gradient at boundary
			virtual tmp<Field<Type> > snGrad() const = 0;

			//- Return face-gradient transform diagonal
			virtual tmp<Field<Type> > snGradTransformDiag() const = 0;

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


	// Member operators

		virtual void operator=(const fvPatchField<Type>&);
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField>
transformFvPatchField<scalar>::gradientInternalCoeffs() const;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "transformFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
