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

#include "cellMotionFvPatchField.H"
#include "fvMesh.H"
#include "volMesh.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(p, iF)
{}


template<class Type>
cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
	const cellMotionFvPatchField<Type>& ptf,
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchField<Type>(ptf, p, iF, mapper)
{}


template<class Type>
cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
	const fvPatch& p,
	const DimensionedField<Type, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchField<Type>(p, iF)
{
	fvPatchField<Type>::operator=(Field<Type>("value", dict, p.size()));
}


template<class Type>
cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
	const cellMotionFvPatchField<Type>& ptf
)
:
	fixedValueFvPatchField<Type>(ptf)
{}


template<class Type>
cellMotionFvPatchField<Type>::cellMotionFvPatchField
(
	const cellMotionFvPatchField<Type>& ptf,
	const DimensionedField<Type, volMesh>& iF
)
:
	fixedValueFvPatchField<Type>(ptf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void cellMotionFvPatchField<Type>::updateCoeffs()
{
	if (this->updated())
	{
		return;
	}

	const fvPatch& p = this->patch();
	const polyPatch& pp = p.patch();
	const fvMesh& mesh = this->dimensionedInternalField().mesh();
	const pointField& points = mesh.points();

	word pfName = this->dimensionedInternalField().name();
	pfName.replace("cell", "point");

	const GeometricField<Type, pointPatchField, pointMesh>& pointMotion =
		this->db().objectRegistry::template
		lookupObject<GeometricField<Type, pointPatchField, pointMesh> >
		(
			pfName
		);

	forAll(p, i)
	{
		this->operator[](i) = pp[i].average(points, pointMotion);
	}

	fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void cellMotionFvPatchField<Type>::write(Ostream& os) const
{
	fvPatchField<Type>::write(os);
	this->writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
