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

#include "inletOutletTotalTemperatureFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF
)
:
	mixedFvPatchScalarField(p, iF),
	UName_("U"),
	phiName_("phi"),
	psiName_("psi"),
	gamma_(0.0),
	T0_(p.size(), 0.0)
{
	this->refValue() = pTraits<scalar>::zero;
	this->refGrad() = pTraits<scalar>::zero;
	this->valueFraction() = 0.0;
}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const dictionary& dict
)
:
	mixedFvPatchScalarField(p, iF),
	UName_(dict.lookupOrDefault<word>("U", "U")),
	phiName_(dict.lookupOrDefault<word>("phi", "phi")),
	psiName_(dict.lookupOrDefault<word>("psi", "psi")),
	gamma_(readScalar(dict.lookup("gamma"))),
	T0_("T0", dict, p.size())
{
	this->refValue() = pTraits<scalar>::zero;
	if (dict.found("value"))
	{
		fvPatchField<scalar>::operator=
		(
			scalarField("value", dict, p.size())
		);
	}
	else
	{
		fvPatchField<scalar>::operator=(T0_);
	}

	this->refGrad() = pTraits<scalar>::zero;
	this->valueFraction() = 0.0;
}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
	const inletOutletTotalTemperatureFvPatchScalarField& ptf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	mixedFvPatchScalarField(ptf, p, iF, mapper),
	UName_(ptf.UName_),
	phiName_(ptf.phiName_),
	psiName_(ptf.psiName_),
	gamma_(ptf.gamma_),
	T0_(ptf.T0_, mapper)
{}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
	const inletOutletTotalTemperatureFvPatchScalarField& tppsf
)
:
	mixedFvPatchScalarField(tppsf),
	UName_(tppsf.UName_),
	phiName_(tppsf.phiName_),
	psiName_(tppsf.psiName_),
	gamma_(tppsf.gamma_),
	T0_(tppsf.T0_)
{}


Foam::inletOutletTotalTemperatureFvPatchScalarField::
inletOutletTotalTemperatureFvPatchScalarField
(
	const inletOutletTotalTemperatureFvPatchScalarField& tppsf,
	const DimensionedField<scalar, volMesh>& iF
)
:
	mixedFvPatchScalarField(tppsf, iF),
	UName_(tppsf.UName_),
	phiName_(tppsf.phiName_),
	psiName_(tppsf.psiName_),
	gamma_(tppsf.gamma_),
	T0_(tppsf.T0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inletOutletTotalTemperatureFvPatchScalarField::autoMap
(
	const fvPatchFieldMapper& m
)
{
	mixedFvPatchScalarField::autoMap(m);
	T0_.autoMap(m);
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::rmap
(
	const fvPatchScalarField& ptf,
	const labelList& addr
)
{
	mixedFvPatchScalarField::rmap(ptf, addr);

	const inletOutletTotalTemperatureFvPatchScalarField& tiptf =
		refCast<const inletOutletTotalTemperatureFvPatchScalarField>(ptf);

	T0_.rmap(tiptf.T0_, addr);
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::updateCoeffs()
{
	if (updated())
	{
		return;
	}

	const fvPatchVectorField& Up =
		lookupPatchField<volVectorField, vector>(UName_);

	const fvsPatchField<scalar>& phip =
		lookupPatchField<surfaceScalarField, scalar>(phiName_);

	const fvPatchField<scalar>& psip =
		lookupPatchField<volScalarField, scalar>(psiName_);

	scalar gM1ByG = (gamma_ - 1.0)/gamma_;

	this->refValue() =
		T0_/(1.0 + 0.5*psip*gM1ByG*(1.0 - pos(phip))*magSqr(Up));
	this->valueFraction() = 1.0 - pos(phip);

	mixedFvPatchScalarField::updateCoeffs();
}


void Foam::inletOutletTotalTemperatureFvPatchScalarField::write(Ostream& os)
const
{
	fvPatchScalarField::write(os);
	if (UName_ != "U")
	{
		os.writeKeyword("U") << UName_ << token::END_STATEMENT << nl;
	}
	if (phiName_ != "phi")
	{
		os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
	}
	if (phiName_ != "psi")
	{
		os.writeKeyword("psi") << psiName_ << token::END_STATEMENT << nl;
	}
	os.writeKeyword("gamma") << gamma_ << token::END_STATEMENT << nl;
	T0_.writeEntry("T0", os);
	writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		inletOutletTotalTemperatureFvPatchScalarField
	);
}

// ************************************************************************* //
