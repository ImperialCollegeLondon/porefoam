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

\*---------------------------------------------------------------------------*/

#include "SlipPointPatchField.H"
#include "transformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	template<class> class MatrixType,
	class Type
>
SlipPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
SlipPointPatchField
(
	const PointPatch& p,
	const DimensionedField<Type, Mesh>& iF
)
:
	BasicSymmetryPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{}


template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	template<class> class MatrixType,
	class Type
>
SlipPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
SlipPointPatchField
(
	const PointPatch& p,
	const DimensionedField<Type, Mesh>& iF,
	const dictionary&
)
:
	BasicSymmetryPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{}


template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	template<class> class MatrixType,
	class Type
>
SlipPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
SlipPointPatchField
(
	const SlipPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>&,
	const PointPatch& p,
	const DimensionedField<Type, Mesh>& iF,
	const PointPatchFieldMapper&
)
:
	BasicSymmetryPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>(p, iF)
{}


template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	template<class> class MatrixType,
	class Type
>
SlipPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
SlipPointPatchField
(
	const SlipPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>& ptf
)
:
	BasicSymmetryPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>(ptf)
{}


template
<
	template<class> class PatchField,
	class Mesh,
	class PointPatch,
	template<class> class MatrixType,
	class Type
>
SlipPointPatchField<PatchField, Mesh, PointPatch, MatrixType, Type>::
SlipPointPatchField
(
	const SlipPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>& ptf,
	const DimensionedField<Type, Mesh>& iF
)
:
	BasicSymmetryPointPatchField
	<PatchField, Mesh, PointPatch, MatrixType, Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
