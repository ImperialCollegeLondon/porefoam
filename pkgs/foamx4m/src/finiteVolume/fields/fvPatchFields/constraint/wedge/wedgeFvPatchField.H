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
	Foam::wedgeFvPatchField

Description
	Foam::wedgeFvPatchField

SourceFiles
	wedgeFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef wedgeFvPatchField_H
#define wedgeFvPatchField_H

#include "transformFvPatchFields.H"
#include "wedgeFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class wedgeFvPatchField
:
	public transformFvPatchField<Type>
{

public:

	//- Runtime type information
	TypeName(wedgeFvPatch::typeName_());


	// Constructors

		//- Construct from patch and internal field
		wedgeFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		wedgeFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given wedgeFvPatchField onto a new patch
		wedgeFvPatchField
		(
			const wedgeFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy
		wedgeFvPatchField
		(
			const wedgeFvPatchField<Type>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new wedgeFvPatchField<Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		wedgeFvPatchField
		(
			const wedgeFvPatchField<Type>&,
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
				new wedgeFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		// Evaluation functions

			//- Return gradient at boundary
			virtual tmp<Field<Type> > snGrad() const;

			//- Evaluate the patch field
			virtual void evaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			);

			//- Return face-gradient transform diagonal
			virtual tmp<Field<Type> > snGradTransformDiag() const;
};


// * * * * * * * * * * * Template Specialisations  * * * * * * * * * * * * * //

template<>
tmp<scalarField > wedgeFvPatchField<scalar>::snGrad() const;

template<>
void wedgeFvPatchField<scalar>::evaluate
(
	const Pstream::commsTypes commsType
);


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "wedgeFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
