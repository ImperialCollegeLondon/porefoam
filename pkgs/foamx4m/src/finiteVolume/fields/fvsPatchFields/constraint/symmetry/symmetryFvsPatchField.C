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

#include "symmetryFvsPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
symmetryFvsPatchField<Type>::symmetryFvsPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, surfaceMesh>& iF
)
:
	fvsPatchField<Type>(p, iF)
{}


template<class Type>
symmetryFvsPatchField<Type>::symmetryFvsPatchField
(
	const symmetryFvsPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, surfaceMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fvsPatchField<Type>(ptf, p, iF, mapper)
{
	if (!isType<symmetryFvPatch>(this->patch()))
	{
		FatalErrorIn
		(
			"symmetryFvsPatchField<Type>::symmetryFvsPatchField\n"
			"(\n"
			"    const symmetryFvsPatchField<Type>& ptf,\n"
			"    const fvPatch& p,\n"
			"    const DimensionedField<Type, surfaceMesh>& iF,\n"
			"    const fvPatchFieldMapper& mapper\n"
			")\n"
		)   << "Field type does not correspond to patch type for patch "
			<< this->patch().index() << "." << endl
			<< "Field type: " << typeName << endl
			<< "Patch type: " << this->patch().type()
			<< exit(FatalError);
	}
}


template<class Type>
symmetryFvsPatchField<Type>::symmetryFvsPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, surfaceMesh>& iF,
	const dictionary& dict
)
:
	fvsPatchField<Type>(p, iF, dict)
{
	if (!isType<symmetryFvPatch>(p))
	{
		FatalIOErrorIn
		(
			"symmetryFvsPatchField<Type>::symmetryFvsPatchField\n"
			"(\n"
			"    const fvPatch& p,\n"
			"    const Field<Type>& field,\n"
			"    const dictionary& dict\n"
			")\n",
			dict
		)   << "patch " << this->patch().index() << " not symmetry type. "
			<< "Patch type = " << p.type()
			<< exit(FatalIOError);
	}
}


template<class Type>
symmetryFvsPatchField<Type>::symmetryFvsPatchField
(
	const symmetryFvsPatchField<Type>& ptf
)
:
	fvsPatchField<Type>(ptf)
{}


template<class Type>
symmetryFvsPatchField<Type>::symmetryFvsPatchField
(
	const symmetryFvsPatchField<Type>& ptf,
	const DimensionedField<Type, surfaceMesh>& iF
)
:
	fvsPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Write
template<class Type>
void symmetryFvsPatchField<Type>::write(Ostream& os) const
{
	fvsPatchField<Type>::write(os);
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
