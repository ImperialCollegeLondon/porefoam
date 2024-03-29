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

#include "uniformFixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(p, iF),
	uniformValue_(pTraits<Type>::zero)
{}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
	const uniformFixedValueFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper&
)
:
	fixedValueFvPatchField<Type>(p, iF),
	uniformValue_(ptf.uniformValue_)
{
	fvPatchField<Type>::operator==(uniformValue_);
}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchField<Type>(p, iF),
	uniformValue_(pTraits<Type>(dict.lookup("uniformValue")))
{
	fvPatchField<Type>::operator==(uniformValue_);
}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
	const uniformFixedValueFvPatchField<Type>& ptf
)
:
	fixedValueFvPatchField<Type>(ptf),
	uniformValue_(ptf.uniformValue_)
{
	fvPatchField<Type>::operator==(uniformValue_);
}


template<class Type>
uniformFixedValueFvPatchField<Type>::uniformFixedValueFvPatchField
(
	const uniformFixedValueFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(ptf, iF),
	uniformValue_(ptf.uniformValue_)
{
	fvPatchField<Type>::operator==(uniformValue_);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void uniformFixedValueFvPatchField<Type>::autoMap
(
	const fvPatchFieldMapper& m
)
{
	this->setSize(m.size());
	fvPatchField<Type>::operator==(uniformValue_);
}


template<class Type>
void uniformFixedValueFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	os.writeKeyword("uniformValue")
		<< uniformValue_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
