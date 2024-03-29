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

#include "mixedFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
mixedFvPatchField<Type>::mixedFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fvPatchField<Type>(p, iF),
	refValue_(p.size()),
	refGrad_(p.size()),
	valueFraction_(p.size())
{}


template<class Type>
mixedFvPatchField<Type>::mixedFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fvPatchField<Type>(p, iF, dict),
	refValue_("refValue", dict, p.size()),
	refGrad_("refGradient", dict, p.size()),
	valueFraction_("valueFraction", dict, p.size())
{
	// Call evaluate only if the value is not found. Used to avoid evaluating
	// when we have incomplete meshes during Parallel Load Balancing. When
	// shipping the field over to another processor, we first call write, making
	// sure that the value is written and read it on the other side (see
	// write member function). If this proves to be problematic, we can always
	// initialize with patch internal field for the start-up. VV, 12/Apr/2019.
	if (!dict.found("value"))
	{
		evaluate();
	}
}


template<class Type>
mixedFvPatchField<Type>::mixedFvPatchField
(
	const mixedFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fvPatchField<Type>(ptf, p, iF, mapper),
	refValue_(ptf.refValue_, mapper),
	refGrad_(ptf.refGrad_, mapper),
	valueFraction_(ptf.valueFraction_, mapper)
{}


template<class Type>
mixedFvPatchField<Type>::mixedFvPatchField
(
	const mixedFvPatchField<Type>& ptf
)
:
	fvPatchField<Type>(ptf),
	refValue_(ptf.refValue_),
	refGrad_(ptf.refGrad_),
	valueFraction_(ptf.valueFraction_)
{}


template<class Type>
mixedFvPatchField<Type>::mixedFvPatchField
(
	const mixedFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fvPatchField<Type>(ptf, iF),
	refValue_(ptf.refValue_),
	refGrad_(ptf.refGrad_),
	valueFraction_(ptf.valueFraction_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void mixedFvPatchField<Type>::autoMap
(
	const fvPatchFieldMapper& m
)
{
	fvPatchField<Type>::autoMap(m);
	refValue_.autoMap(m);
	refGrad_.autoMap(m);
	valueFraction_.autoMap(m);
}


template<class Type>
void mixedFvPatchField<Type>::rmap
(
	const fvPatchField<Type>& ptf,
	const labelList& addr
)
{
	fvPatchField<Type>::rmap(ptf, addr);

	const mixedFvPatchField<Type>& mptf =
		refCast<const mixedFvPatchField<Type> >(ptf);

	refValue_.rmap(mptf.refValue_, addr);
	refGrad_.rmap(mptf.refGrad_, addr);
	valueFraction_.rmap(mptf.valueFraction_, addr);
}


template<class Type>
void mixedFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
	if (!this->updated())
	{
		this->updateCoeffs();
	}

	Field<Type>::operator=
	(
		valueFraction_*refValue_
	  + (1.0 - valueFraction_)*
		(
			this->patchInternalField()
		  + refGrad_/this->patch().deltaCoeffs()
		)
	);

	fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > mixedFvPatchField<Type>::snGrad() const
{
	return
		valueFraction_
	   *(refValue_ - this->patchInternalField())
	   *this->patch().deltaCoeffs()
	  + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
tmp<Field<Type> > mixedFvPatchField<Type>::valueInternalCoeffs
(
	const tmp<scalarField>&
) const
{
	return pTraits<Type>::one*(1.0 - valueFraction_);

}


template<class Type>
tmp<Field<Type> > mixedFvPatchField<Type>::valueBoundaryCoeffs
(
	const tmp<scalarField>&
) const
{
	return
		 valueFraction_*refValue_
	   + (1.0 - valueFraction_)*refGrad_/this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > mixedFvPatchField<Type>::gradientInternalCoeffs() const
{
	return -pTraits<Type>::one*valueFraction_*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > mixedFvPatchField<Type>::gradientBoundaryCoeffs() const
{
	return
		valueFraction_*this->patch().deltaCoeffs()*refValue_
	  + (1.0 - valueFraction_)*refGrad_;
}


template<class Type>
void mixedFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	refValue_.writeEntry("refValue", os);
	refGrad_.writeEntry("refGradient", os);
	valueFraction_.writeEntry("valueFraction", os);
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
