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
	ZeroGradientPointPatchField

Description

SourceFiles
	ZeroGradientPointPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef ZeroGradientPointPatchField_H
#define ZeroGradientPointPatchField_H

#include "PointPatchFieldTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{


template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	template<class> class MatrixType,
	class Type
>
class ZeroGradientPointPatchField
:
	public PatchField<Type>
{

public:

	//- Runtime type information
	TypeName("zeroGradient");


	// Constructors

		//- Construct from patch and internal field
		ZeroGradientPointPatchField
		(
			const PointPatch&,
			const DimensionedField<Type, Mesh>&
		);

		//- Construct from patch, internal field and dictionary
		ZeroGradientPointPatchField
		(
			const PointPatch&,
			const DimensionedField<Type, Mesh>&,
			const dictionary&
		);

		//- Construct by mapping given patchField<Type> onto a new patch
		ZeroGradientPointPatchField
		(
			const ZeroGradientPointPatchField
				<PatchField, Mesh, PointPatch, MatrixType, Type>&,
			const PointPatch&,
			const DimensionedField<Type, Mesh>&,
			const PointPatchFieldMapper&
		);

		//- Construct as copy
		ZeroGradientPointPatchField
		(
			const ZeroGradientPointPatchField
			<PatchField, Mesh, PointPatch, MatrixType, Type>&
		);

		//- Construct and return a clone
		virtual autoPtr<PatchField<Type> > clone() const
		{
			return autoPtr<PatchField<Type> >
			(
				new ZeroGradientPointPatchField
				<PatchField, Mesh, PointPatch, MatrixType, Type>(*this)
			);
		}

		//- Construct as copy setting internal field reference
		ZeroGradientPointPatchField
		(
			const ZeroGradientPointPatchField
			<PatchField, Mesh, PointPatch, MatrixType, Type>&,
			const DimensionedField<Type, Mesh>&
		);

		//- Construct and return a clone setting internal field reference
		virtual autoPtr<PatchField<Type> > clone
		(
			const DimensionedField<Type, Mesh>& iF
		) const
		{
			return autoPtr<PatchField<Type> >
			(
				new ZeroGradientPointPatchField
				<PatchField, Mesh, PointPatch, MatrixType, Type>(*this, iF)
			);
		}


	// Member functions

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#	 include "ZeroGradientPointPatchField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
