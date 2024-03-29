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

#include "fvPatchFields.H"
#include "constantAlphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::constantAlphaContactAngleFvPatchScalarField::
constantAlphaContactAngleFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF
)
:
	alphaContactAngleFvPatchScalarField(p, iF),
	theta0_(0.0)
{}


Foam::constantAlphaContactAngleFvPatchScalarField::
constantAlphaContactAngleFvPatchScalarField
(
	const constantAlphaContactAngleFvPatchScalarField& gcpsf,
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const fvPatchFieldMapper& mapper
)
:
	alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
	theta0_(gcpsf.theta0_)
{}


Foam::constantAlphaContactAngleFvPatchScalarField::
constantAlphaContactAngleFvPatchScalarField
(
	const fvPatch& p,
	const DimensionedField<scalar, volMesh>& iF,
	const dictionary& dict
)
:
	alphaContactAngleFvPatchScalarField(p, iF),
	theta0_(readScalar(dict.lookup("theta0")))
{
	evaluate();
}


Foam::constantAlphaContactAngleFvPatchScalarField::
constantAlphaContactAngleFvPatchScalarField
(
	const constantAlphaContactAngleFvPatchScalarField& gcpsf
)
:
	alphaContactAngleFvPatchScalarField(gcpsf),
	theta0_(gcpsf.theta0_)
{}


Foam::constantAlphaContactAngleFvPatchScalarField::
constantAlphaContactAngleFvPatchScalarField
(
	const constantAlphaContactAngleFvPatchScalarField& gcpsf,
	const DimensionedField<scalar, volMesh>& iF
)
:
	alphaContactAngleFvPatchScalarField(gcpsf, iF),
	theta0_(gcpsf.theta0_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::constantAlphaContactAngleFvPatchScalarField::theta
(
	const fvPatchVectorField&,
	const fvsPatchVectorField&
) const
{
	return tmp<scalarField>(new scalarField(size(), theta0_));
}


void Foam::constantAlphaContactAngleFvPatchScalarField::write
(
	Ostream& os
) const
{
	fvPatchScalarField::write(os);
	os.writeKeyword("theta0") << theta0_ << token::END_STATEMENT << nl;
	writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
	makePatchTypeField
	(
		fvPatchScalarField,
		constantAlphaContactAngleFvPatchScalarField
	);
}

// ************************************************************************* //
