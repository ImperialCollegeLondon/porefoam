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

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fvPatchField<Type>(p, iF)
{}


template<class Type>
fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fvPatchField<Type>(p, iF, dict, true)
{}


template<class Type>
fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
	const fixedValueFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
	const fixedValueFvPatchField<Type>& ptf
)
:
	fvPatchField<Type>(ptf)
{}


template<class Type>
fixedValueFvPatchField<Type>::fixedValueFvPatchField
(
	const fixedValueFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::valueInternalCoeffs
(
	const tmp<scalarField>&
) const
{
	return tmp<Field<Type> >
	(
		new Field<Type>(this->size(), pTraits<Type>::zero)
	);
}


template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::valueBoundaryCoeffs
(
	const tmp<scalarField>&
) const
{
	return *this;
}


template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::gradientInternalCoeffs() const
{
	return -pTraits<Type>::one*this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > fixedValueFvPatchField<Type>::gradientBoundaryCoeffs() const
{
	return this->patch().deltaCoeffs()*(*this);
}


template<class Type>
void fixedValueFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
