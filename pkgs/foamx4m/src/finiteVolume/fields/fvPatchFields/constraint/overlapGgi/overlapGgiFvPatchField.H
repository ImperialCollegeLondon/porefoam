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
	overlapGgiFvPatchField

Description
	Partial overlap generalized grid interface patch field.

Author
	Hrvoje Jasak, Wikki Ltd.  All rights reserved.

SourceFiles
	overlapGgiFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef overlapGgiFvPatchField_H
#define overlapGgiFvPatchField_H

#include "coupledFvPatchField.H"
#include "overlapGGILduInterfaceField.H"
#include "overlapGgiFvPatch.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template<class Type>
class overlapGgiFvPatchField
:
	public overlapGGILduInterfaceField,
	public coupledFvPatchField<Type>
{
	// Private data

		//- Local reference cast into the overlapGgi patch
		const overlapGgiFvPatch& overlapGgiPatch_;


public:

	//- Runtime type information
	TypeName(overlapGgiFvPatch::typeName_());


	// Constructors

		//- Construct from patch and internal field
		overlapGgiFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct from patch, internal field and dictionary
		overlapGgiFvPatchField
		(
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const dictionary&
		);

		//- Construct by mapping given overlapGgiFvPatchField onto a new patch
		overlapGgiFvPatchField
		(
			const overlapGgiFvPatchField<Type>&,
			const fvPatch&,
			const DimensionedField<Type, volMesh>&,
			const fvPatchFieldMapper&
		);

		//- Construct as copy setting internal field reference
		overlapGgiFvPatchField
		(
			const overlapGgiFvPatchField<Type>&,
			const DimensionedField<Type, volMesh>&
		);

		//- Construct and return a clone
		virtual tmp<fvPatchField<Type> > clone() const
		{
			return tmp<fvPatchField<Type> >
			(
				new overlapGgiFvPatchField<Type>(*this)
			);
		}

		//- Construct and return a clone setting internal field reference
		virtual tmp<fvPatchField<Type> > clone
		(
			const DimensionedField<Type, volMesh>& iF
		) const
		{
			return tmp<fvPatchField<Type> >
			(
				new overlapGgiFvPatchField<Type>(*this, iF)
			);
		}


	// Member functions

		//- Access

			// Return overlapGgiFvPatch
			const overlapGgiFvPatch& overlapGgiPatch() const
			{
				return overlapGgiPatch_;
			}


		// Evaluation functions

			//- Return neighbour field given internal cell data
			//  This is a dummy function: returning patch values
			//  due to scheduling issues in parallel communication.
			//  See processofFvPatchField for details HJ, 11/May/2009
			virtual tmp<Field<Type> > patchNeighbourField() const;

			//- Return neighbour field of internal field, untransformed
			virtual tmp<scalarField>
			untransformedInterpolate(const direction cmpt) const;

			//- Initialise the evaluation of the patch field
			virtual void initEvaluate(const Pstream::commsTypes commsType);

			//- Evaluate the patch field
			virtual void evaluate
			(
				const Pstream::commsTypes commsType
			);

			//- Transform given patch component field
			virtual void transformCoupleField
			(
				scalarField& f,
				const direction cmpt
			) const
			{
				overlapGGILduInterfaceField::transformCoupleField(f, cmpt);
			}

			//- Initialise neighbour matrix update
			virtual void initInterfaceMatrixUpdate
			(
				const scalarField& psiInternal,
				scalarField& result,
				const lduMatrix& m,
				const scalarField& coeffs,
				const direction cmpt,
				const Pstream::commsTypes commsType,
				const bool switchToLhs
			) const;

			//- Update result field based on interface functionality
			virtual void updateInterfaceMatrix
			(
				const scalarField& psiInternal,
				scalarField& result,
				const lduMatrix&,
				const scalarField& coeffs,
				const direction cmpt,
				const Pstream::commsTypes commsType,
				const bool switchToLhs
			) const;


		//- GGI coupled interface functions

			//- Does the patch field perform the transfromation
			virtual bool doTransform() const
			{
				return false;
			}

			//- Return face transformation tensor
			virtual const tensorField& forwardT() const
			{
				return overlapGgiPatch_.forwardT();
			}

			//- Return neighbour-cell transformation tensor
			virtual const tensorField& reverseT() const
			{
				return overlapGgiPatch_.reverseT();
			}

			//- Return rank of component for transform
			virtual int rank() const
			{
				return pTraits<Type>::rank;
			}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	include "overlapGgiFvPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
