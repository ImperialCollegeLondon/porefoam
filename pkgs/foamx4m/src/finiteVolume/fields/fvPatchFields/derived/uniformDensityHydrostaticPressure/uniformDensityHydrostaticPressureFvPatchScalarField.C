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

#include "uniformDensityHydrostaticPressureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "uniformDimensionedFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF
)
:
	fixedValueFvPatchScalarField(p, iF),
	rho_(0.0),
	pRefValue_(0.0),
	pRefPoint_(vector::zero)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchScalarField(p, iF),
	rho_(readScalar(dict.lookup("rho"))),
	pRefValue_(readScalar(dict.lookup("pRefValue"))),
	pRefPoint_(dict.lookup("pRefPoint"))
{
	if (dict.found("value"))
	{
		fvPatchField<scalar>::operator=
		(
			scalarField("value", dict, p.size())
		);
	}
	else
	{
		// Evaluation not allowed here since g may not be available
		// Bug fix, HJ, 19/May/2011
		fvPatchField<scalar>::operator=(patchInternalField());
//		 evaluate();
	}
}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
	const uniformDensityHydrostaticPressureFvPatchScalarField& ptf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchScalarField(ptf, p, iF, mapper),
	rho_(ptf.rho_),
	pRefValue_(ptf.pRefValue_),
	pRefPoint_(ptf.pRefPoint_)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
	const uniformDensityHydrostaticPressureFvPatchScalarField& ptf
)
:
	fixedValueFvPatchScalarField(ptf),
	rho_(ptf.rho_),
	pRefValue_(ptf.pRefValue_),
	pRefPoint_(ptf.pRefPoint_)
{}


Foam::uniformDensityHydrostaticPressureFvPatchScalarField::
uniformDensityHydrostaticPressureFvPatchScalarField
(
	const uniformDensityHydrostaticPressureFvPatchScalarField& ptf,
	const DimensionedField<scalar, volMesh>& iF
)
:
	fixedValueFvPatchScalarField(ptf, iF),
	rho_(ptf.rho_),
	pRefValue_(ptf.pRefValue_),
	pRefPoint_(ptf.pRefPoint_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::uniformDensityHydrostaticPressureFvPatchScalarField::updateCoeffs()
{
	if (updated())
	{
		return;
	}

	const uniformDimensionedVectorField& g =
		db().lookupObject<uniformDimensionedVectorField>("g");

	operator==
	(
		pRefValue_
	  + rho_*((g.value() & patch().Cf()) - (g.value() & pRefPoint_))
	);

	fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::uniformDensityHydrostaticPressureFvPatchScalarField::write
(
	Ostream& os
) const
{
	fvPatchScalarField::write(os);
	os.writeKeyword("rho") << rho_ << token::END_STATEMENT << nl;
	os.writeKeyword("pRefValue") << pRefValue_ << token::END_STATEMENT << nl;
	os.writeKeyword("pRefPoint") << pRefPoint_ << token::END_STATEMENT << nl;
	writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		uniformDensityHydrostaticPressureFvPatchScalarField
	);
}

// ************************************************************************* //
