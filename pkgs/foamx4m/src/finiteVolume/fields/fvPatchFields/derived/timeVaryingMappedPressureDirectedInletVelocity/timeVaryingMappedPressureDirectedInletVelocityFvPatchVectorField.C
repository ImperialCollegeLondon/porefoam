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

#include "objectRegistry.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvPatchFieldMapper.H"
#include "timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
(
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF
)
:
	timeVaryingMappedFixedValueFvPatchVectorField(p, iF),
	phiName_("phi"),
	rhoName_("rho")
{}


Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
(
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF,
	const dictionary& dict
)
:
	timeVaryingMappedFixedValueFvPatchVectorField(p, iF),
	phiName_(dict.lookupOrDefault<word>("phi", "phi")),
	rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
	fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}


Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
(
	const timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField& ptf,
	const fvPatch& p,
	const DimensionedField<vector, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	timeVaryingMappedFixedValueFvPatchVectorField(ptf, p, iF, mapper),
	phiName_(ptf.phiName_),
	rhoName_(ptf.rhoName_)
{}


Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
(
	const timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField&
	pivpvf
)
:
	timeVaryingMappedFixedValueFvPatchVectorField(pivpvf),
	phiName_(pivpvf.phiName_),
	rhoName_(pivpvf.rhoName_)
{}


Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
(
	const timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField&
	pivpvf,
	const DimensionedField<vector, volMesh>& iF
)
:
	timeVaryingMappedFixedValueFvPatchVectorField(pivpvf, iF),
	phiName_(pivpvf.phiName_),
	rhoName_(pivpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
autoMap
(
	const fvPatchFieldMapper& m
)
{
	timeVaryingMappedFixedValueFvPatchVectorField::autoMap(m);
}


void Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
rmap
(
	const fvPatchVectorField& ptf,
	const labelList& addr
)
{
	timeVaryingMappedFixedValueFvPatchVectorField::rmap(ptf, addr);
}


void Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
updateCoeffs()
{
	if (updated())
	{
		return;
	}

	// Map the tabulated velocity field onto this field
	timeVaryingMappedFixedValueFvPatchVectorField::updateCoeffs();
	vectorField inletDir = *this;

	// Normalise to obtain the flow direction
	inletDir /= (mag(inletDir) + SMALL);

	const surfaceScalarField& phi =
		db().lookupObject<surfaceScalarField>(phiName_);

	const fvsPatchField<scalar>& phip =
		patch().patchField<surfaceScalarField, scalar>(phi);

	vectorField n = patch().nf();
	scalarField ndmagS = (n & inletDir)*patch().magSf();

	if (phi.dimensions() == dimVelocity*dimArea)
	{
		operator==(inletDir*phip/ndmagS);
	}
	else if (phi.dimensions() == dimDensity*dimVelocity*dimArea)
	{
		const fvPatchField<scalar>& rhop =
			lookupPatchField<volScalarField, scalar>(rhoName_);

		operator==(inletDir*phip/(rhop*ndmagS));
	}
	else
	{
		FatalErrorIn
		(
			"timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField"
			"::updateCoeffs()"
		)   << "dimensions of phi are not correct"
			<< "\n	on patch " << this->patch().name()
			<< " of field " << this->dimensionedInternalField().name()
			<< " in file " << this->dimensionedInternalField().objectPath()
			<< exit(FatalError);
	}
}


void Foam::timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField::
write
(
	Ostream& os
) const
{
	if (phiName_ != "phi")
	{
		os.writeKeyword("phi") << phiName_ << token::END_STATEMENT << nl;
	}
	if (rhoName_ != "rho")
	{
		os.writeKeyword("rho") << rhoName_ << token::END_STATEMENT << nl;
	}
	timeVaryingMappedFixedValueFvPatchVectorField::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	makePatchTypeField
	(
		fvPatchVectorField,
		timeVaryingMappedPressureDirectedInletVelocityFvPatchVectorField
	);
}


// ************************************************************************* //
