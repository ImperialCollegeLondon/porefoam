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

#include "fixedGradientFvPatchField.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fvPatchField<Type>(p, iF),
	gradient_(p.size(), pTraits<Type>::zero)
{}


template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
	const fixedGradientFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fvPatchField<Type>(ptf, p, iF, mapper),
	gradient_(ptf.gradient_, mapper)
{}


template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fvPatchField<Type>(p, iF, dict),
	gradient_("gradient", dict, p.size())
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
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
	const fixedGradientFvPatchField<Type>& ptf
)
:
	fvPatchField<Type>(ptf),
	gradient_(ptf.gradient_)
{}


template<class Type>
fixedGradientFvPatchField<Type>::fixedGradientFvPatchField
(
	const fixedGradientFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fvPatchField<Type>(ptf, iF),
	gradient_(ptf.gradient_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void fixedGradientFvPatchField<Type>::autoMap
(
	const fvPatchFieldMapper& m
)
{
	fvPatchField<Type>::autoMap(m);
	gradient_.autoMap(m);
}


template<class Type>
void fixedGradientFvPatchField<Type>::rmap
(
	const fvPatchField<Type>& ptf,
	const labelList& addr
)
{
	fvPatchField<Type>::rmap(ptf, addr);

	const fixedGradientFvPatchField<Type>& fgptf =
		refCast<const fixedGradientFvPatchField<Type> >(ptf);

	gradient_.rmap(fgptf.gradient_, addr);
}


template<class Type>
void fixedGradientFvPatchField<Type>::evaluate(const Pstream::commsTypes)
{
	if (!this->updated())
	{
		this->updateCoeffs();
	}

	Field<Type>::operator=
	(
		this->patchInternalField() + gradient_/this->patch().deltaCoeffs()
	);

	fvPatchField<Type>::evaluate();
}


template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::valueInternalCoeffs
(
	const tmp<scalarField>&
) const
{
	return tmp<Field<Type> >
	(
		new Field<Type>(this->size(), pTraits<Type>::one)
	);
}


template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::valueBoundaryCoeffs
(
	const tmp<scalarField>&
) const
{
	return gradient()/this->patch().deltaCoeffs();
}


template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::
gradientInternalCoeffs() const
{
	return tmp<Field<Type> >
	(
		new Field<Type>(this->size(), pTraits<Type>::zero)
	);
}


template<class Type>
tmp<Field<Type> > fixedGradientFvPatchField<Type>::
gradientBoundaryCoeffs() const
{
	return gradient();
}


template<class Type>
void fixedGradientFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	gradient_.writeEntry("gradient", os);

	// Write value along with the gradient in order to avoid calling evaluate
	// on construction (from dictionary constructor) for Parallel Load Balancing
	// runs. See dictionary constructor for details. VV, 12/Apr/2019.
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
