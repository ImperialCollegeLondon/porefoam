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

#include "directionMixedFvPatchField.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	transformFvPatchField<Type>(p, iF),
	refValue_(p.size()),
	refGrad_(p.size()),
	valueFraction_(p.size())
{}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	transformFvPatchField<Type>(p, iF, dict),
	refValue_("refValue", dict, p.size()),
	refGrad_("refGradient", dict, p.size()),
	valueFraction_("valueFraction", dict, p.size())
{
	evaluate();
}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
	const directionMixedFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	transformFvPatchField<Type>(ptf, p, iF, mapper),
	refValue_(ptf.refValue_, mapper),
	refGrad_(ptf.refGrad_, mapper),
	valueFraction_(ptf.valueFraction_, mapper)
{}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
	const directionMixedFvPatchField<Type>& ptf
)
:
	transformFvPatchField<Type>(ptf),
	refValue_(ptf.refValue_),
	refGrad_(ptf.refGrad_),
	valueFraction_(ptf.valueFraction_)
{}


template<class Type>
directionMixedFvPatchField<Type>::directionMixedFvPatchField
(
	const directionMixedFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	transformFvPatchField<Type>(ptf, iF),
	refValue_(ptf.refValue_),
	refGrad_(ptf.refGrad_),
	valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void directionMixedFvPatchField<Type>::autoMap
(
	const fvPatchFieldMapper& m
)
{
	transformFvPatchField<Type>::autoMap(m);
	refValue_.autoMap(m);
	refGrad_.autoMap(m);
	valueFraction_.autoMap(m);
}


template<class Type>
void directionMixedFvPatchField<Type>::rmap
(
	const fvPatchField<Type>& ptf,
	const labelList& addr
)
{
	transformFvPatchField<Type>::rmap(ptf, addr);

	const directionMixedFvPatchField<Type>& dmptf =
		refCast<const directionMixedFvPatchField<Type> >(ptf);

	refValue_.rmap(dmptf.refValue_, addr);
	refGrad_.rmap(dmptf.refGrad_, addr);
	valueFraction_.rmap(dmptf.valueFraction_, addr);
}


template<class Type>
tmp<Field<Type> > directionMixedFvPatchField<Type>::snGrad() const
{
	Field<Type> pif = this->patchInternalField();

	Field<Type> normalValue = transform(valueFraction_, refValue_);

	Field<Type> gradValue = pif + refGrad_/this->patch().deltaCoeffs();

	Field<Type> transformGradValue =
		transform(I - valueFraction_, gradValue);

	return
		(normalValue + transformGradValue - pif)*
		this->patch().deltaCoeffs();
}


template<class Type>
void directionMixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
	if (!this->updated())
	{
		this->updateCoeffs();
	}

	Field<Type> normalValue = transform(valueFraction_, refValue_);

	Field<Type> gradValue =
		this->patchInternalField() + refGrad_/this->patch().deltaCoeffs();

	Field<Type> transformGradValue =
		transform(I - valueFraction_, gradValue);

	Field<Type>::operator=(normalValue + transformGradValue);

	transformFvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > directionMixedFvPatchField<Type>::snGradTransformDiag() const
{
	vectorField diag(valueFraction_.size());

	diag.replace
	(
		vector::X,
		sqrt(mag(valueFraction_.component(symmTensor::XX)))
	);
	diag.replace
	(
		vector::Y,
		sqrt(mag(valueFraction_.component(symmTensor::YY)))
	);
	diag.replace
	(
		vector::Z,
		sqrt(mag(valueFraction_.component(symmTensor::ZZ)))
	);

	return transformFieldMask<Type>(pow<vector, pTraits<Type>::rank>(diag));
}


template<class Type>
void directionMixedFvPatchField<Type>::write(Ostream& os) const
{
	transformFvPatchField<Type>::write(os);
	refValue_.writeEntry("refValue", os);
	refGrad_.writeEntry("refGradient", os);
	valueFraction_.writeEntry("valueFraction", os);
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
