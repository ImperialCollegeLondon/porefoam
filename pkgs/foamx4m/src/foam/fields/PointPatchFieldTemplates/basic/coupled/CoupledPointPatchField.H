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
	CoupledPointPatchField

Description

SourceFiles
	CoupledPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef CoupledPointPatchField_H
#define CoupledPointPatchField_H

#include "lduInterfaceField.H"
#include "PointPatchFieldTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	class CoupledPointPatch,
	template<class> class MatrixType,
	class Type
>
class CoupledPointPatchField
:
	public lduInterfaceField,
	public PatchField<Type>
{

public:

	//- Runtime type information
	TypeName(CoupledPointPatch::typeName_());


	// Constructors

		//- Construct from patch and internal field
		CoupledPointPatchField
		(
			const PointPatch&,
			const DimensionedField<Type, Mesh>&
		);

		//- Construct by mapping given patchField<Type> onto a new patch
		CoupledPointPatchField
		(
			const CoupledPointPatchField
			<
				PatchField,
				Mesh,
				PointPatch,
				CoupledPointPatch,
				MatrixType,
				Type
			>&,
			const PointPatch&,
			const DimensionedField<Type, Mesh>&,
			const PointPatchFieldMapper&
		);

		//- Construct as copy
		CoupledPointPatchField
		(
			const CoupledPointPatchField
			<
				PatchField,
				Mesh,
				PointPatch,
				CoupledPointPatch,
				MatrixType,
				Type
			>&
		);

		//- Construct and return a clone
		virtual autoPtr<PatchField<Type> > clone() const = 0;

		//- Construct as copy setting internal field reference
		CoupledPointPatchField
		(
			const CoupledPointPatchField
			<
				PatchField,
				Mesh,
				PointPatch,
				CoupledPointPatch,
				MatrixType,
				Type
			>&,
			const DimensionedField<Type, Mesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual autoPtr<PatchField<Type> > clone
		(
			const DimensionedField<Type, Mesh>& iF
		) const = 0;


	// Member functions

		// Access

			//- Return true if this patch field is derived from
			//  CoupledPointPatchField
			virtual bool coupled() const
			{
				return true;
			}

			//- Evaluate the patch field
			virtual void evaluate
			(
				const Pstream::commsTypes commsType = Pstream::blocking
			) = 0;

			//- Add field to internal field
			virtual void addField(Field<Type>&) const = 0;


		// Matrix construction (completion) functionality

			//- Add diagonal/source contributions
			virtual void addDiag(scalarField&) const = 0;
			virtual void addSource(scalarField&) const = 0;

			//- Add upper/lower contributions
			virtual void addUpperLower(scalarField&) const = 0;

			virtual void eliminateUpperLower(scalarField&) const = 0;

			//- Get the cut edge coefficients in Amul order
			virtual tmp<scalarField> cutBouCoeffs(const lduMatrix&) const = 0;

			//- Get the cut edge coefficients in Tmul order
			virtual tmp<scalarField> cutIntCoeffs(const lduMatrix&) const = 0;


		// Coupled interface functionality

			//- Transform neighbour field
			virtual void transformCoupleField
			(
				scalarField& pnf,
				const direction cmpt
			) const
			{}

			//- Initialise neighbour matrix update
			virtual void initInterfaceMatrixUpdate
			(
				const scalarField&,
				scalarField&,
				const lduMatrix&,
				const scalarField&,
				const direction,
				const Pstream::commsTypes commsType,
				const bool switchToLhs
			) const
			{}

			//- Update result field based on interface functionality
			virtual void updateInterfaceMatrix
			(
				const scalarField&,
				scalarField&,
				const lduMatrix&,
				const scalarField&,
				const direction,
				const Pstream::commsTypes commsType,
				const bool switchToLhs
			) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	 include "CoupledPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
