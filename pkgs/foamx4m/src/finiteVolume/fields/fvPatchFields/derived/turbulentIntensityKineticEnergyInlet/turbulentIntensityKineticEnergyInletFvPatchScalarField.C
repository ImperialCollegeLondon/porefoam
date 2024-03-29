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

#include "turbulentIntensityKineticEnergyInletFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF
)
:
	inletOutletFvPatchScalarField(p, iF),
	intensity_(0.0),
	UName_("undefined-U"),
	phiName_("undefined-phi")
{
	this->refValue() = 0.0;
	this->refGrad() = 0.0;
	this->valueFraction() = 0.0;
}


Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
	const turbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	inletOutletFvPatchScalarField(ptf, p, iF, mapper),
	intensity_(ptf.intensity_),
	UName_(ptf.UName_),
	phiName_(ptf.phiName_)
{}


Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const dictionary& dict
)
:
	inletOutletFvPatchScalarField(p, iF),
	intensity_(readScalar(dict.lookup("intensity"))),
	UName_(dict.lookupOrDefault<word>("U", "U")),
	phiName_(dict.lookupOrDefault<word>("phi", "phi"))
{
	if (intensity_ < 0 || intensity_ > 1)
	{
		FatalErrorIn
		(
			"turbulentIntensityKineticEnergyInletFvPatchScalarField::"
			"turbulentIntensityKineticEnergyInletFvPatchScalarField"
			"(const fvPatch& p, const DimensionedField<scalar, volMesh>& iF, "
			"const dictionary& dict)"
		)   << "Turbulence intensity should be specified as a fraction 0-1 "
			   "of the mean velocity\n"
			   "    value given is " << intensity_
			<< "\n    on patch " << this->patch().name()
			<< " of field " << this->dimensionedInternalField().name()
			<< " in file " << this->dimensionedInternalField().objectPath()
			<< exit(FatalError);
	}

	fvPatchScalarField::operator=(scalarField("value", dict, p.size()));

	this->refValue() = 0.0;
	this->refGrad() = 0.0;
	this->valueFraction() = 0.0;
}


Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
	const turbulentIntensityKineticEnergyInletFvPatchScalarField& ptf
)
:
	inletOutletFvPatchScalarField(ptf),
	intensity_(ptf.intensity_),
	UName_(ptf.UName_),
	phiName_(ptf.phiName_)
{}


Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
turbulentIntensityKineticEnergyInletFvPatchScalarField
(
	const turbulentIntensityKineticEnergyInletFvPatchScalarField& ptf,
	const DimensionedField<scalar, volMesh>& iF
)
:
	inletOutletFvPatchScalarField(ptf, iF),
	intensity_(ptf.intensity_),
	UName_(ptf.UName_),
	phiName_(ptf.phiName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::
updateCoeffs()
{
	if (updated())
	{
		return;
	}

	const fvPatchVectorField& Up =
		lookupPatchField<volVectorField, vector>(UName_);

	const fvsPatchScalarField& phip =
		lookupPatchField<surfaceScalarField, scalar>(phiName_);

	this->refValue() = SMALL + 1.5*sqr(intensity_)*magSqr(Up);
	this->valueFraction() = neg(phip);

	inletOutletFvPatchScalarField::updateCoeffs();
}


void Foam::turbulentIntensityKineticEnergyInletFvPatchScalarField::write
(
	Ostream& os
) const
{
	fvPatchScalarField::write(os);
	os.writeKeyword("intensity") << intensity_ << token::END_STATEMENT << nl;
	os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
	os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
	writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		turbulentIntensityKineticEnergyInletFvPatchScalarField
	);
}

// ************************************************************************* //
