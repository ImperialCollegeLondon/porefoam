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

#include "pressureInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF
)
:
	fixedValueFvPatchVectorField(p, iF),
	phiName_("phi"),
	rhoName_("rho")
{}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
	const pressureInletVelocityFvPatchVectorField& ptf,
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	fixedValueFvPatchVectorField(ptf, p, iF, mapper),
	phiName_(ptf.phiName_),
	rhoName_(ptf.rhoName_)
{}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF,
	const dictionary& dict
)
:
	fixedValueFvPatchVectorField(p, iF),
	phiName_(dict.lookupOrDefault<word>("phi", "phi")),
	rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
	fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
	const pressureInletVelocityFvPatchVectorField& pivpvf
)
:
	fixedValueFvPatchVectorField(pivpvf),
	phiName_(pivpvf.phiName_),
	rhoName_(pivpvf.rhoName_)
{}


Foam::pressureInletVelocityFvPatchVectorField::
pressureInletVelocityFvPatchVectorField
(
	const pressureInletVelocityFvPatchVectorField& pivpvf,
	const DimensionedField<vector, volMesh>& iF
)
:
	fixedValueFvPatchVectorField(pivpvf, iF),
	phiName_(pivpvf.phiName_),
	rhoName_(pivpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::pressureInletVelocityFvPatchVectorField::updateCoeffs()
{
	if (updated())
	{
		return;
	}

	if (!this->db().objectRegistry::found(phiName_))
	{
		// Flux not available, do not update
		InfoIn
		(
			"void pressureInletVelocityFvPatchVectorField::updateCoeffs()"
		)   << "Flux field " << phiName_ << " not found.  "
			<< "Performing fixed value update" << endl;

		fixedValueFvPatchVectorField::updateCoeffs();

		return;
	}

	const surfaceScalarField& phi =
		db().lookupObject<surfaceScalarField>(phiName_);

	const fvsPatchScalarField& phip =
		patch().patchField<surfaceScalarField, scalar>(phi);

	vectorField n = patch().nf();
	const Field<scalar>& magS = patch().magSf();

	if (phi.dimensions() == dimVelocity*dimArea)
	{
		operator==(n*phip/magS);
	}
	else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
	{
		if (!this->db().objectRegistry::found(rhoName_))
		{
			// Rho not available, do not update
			fixedValueFvPatchVectorField::updateCoeffs();

			return;
		}

		const fvPatchField<scalar>& rhop =
			lookupPatchField<volScalarField, scalar>(rhoName_);

		operator==(n*phip/(rhop*magS));
	}
	else
	{
		FatalErrorIn("pressureInletVelocityFvPatchVectorField::updateCoeffs()")
			<< "dimensions of phi are not correct"
			<< "\n    on patch " << this->patch().name()
			<< " of field " << this->dimensionedInternalField().name()
			<< " in file " << this->dimensionedInternalField().objectPath()
			<< exit(FatalError);
	}

	fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::pressureInletVelocityFvPatchVectorField::write(Ostream& os) const
{
	fvPatchVectorField::write(os);
	if (phiName_ != "phi")
	{
		os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
	}
	if (rhoName_ != "rho")
	{
		os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
	}
	writeEntry("value", os);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

void Foam::pressureInletVelocityFvPatchVectorField::operator=
(
	const fvPatchField<vector>& pvf
)
{
	fvPatchField<vector>::operator=(patch().nf()*(patch().nf() & pvf));
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

makePatchTypeField
(
	fvPatchVectorField,
	pressureInletVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
