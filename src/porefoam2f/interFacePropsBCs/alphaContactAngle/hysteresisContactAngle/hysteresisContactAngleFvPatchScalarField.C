/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2010 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
	This file is partly based on OpenFOAM code.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fvPatchFields.H"
#include "hysteresisContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::hysteresisContactAngleFvPatchScalarField::
hysteresisContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(p, iF),
    thetaA_(0.0),
    //uTheta_(0.0),
    thetaR_(0.0),
    thetaAOW_(0.0),
    thetaROW_(0.0),
    maxAlphaHist(iF)
{}


Foam::hysteresisContactAngleFvPatchScalarField::
hysteresisContactAngleFvPatchScalarField
(
    const hysteresisContactAngleFvPatchScalarField& gcpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, p, iF, mapper),
    thetaA_(gcpsf.thetaA_),
    //uTheta_(gcpsf.uTheta_),
    thetaR_(gcpsf.thetaR_),
    thetaAOW_(gcpsf.thetaAOW_),
    thetaROW_(gcpsf.thetaROW_),
    maxAlphaHist(iF)
{}


Foam::hysteresisContactAngleFvPatchScalarField::
hysteresisContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    alphaContactAngleFvPatchScalarField(p, iF, dict),
    thetaA_(readScalar(dict.lookup("thetaA"))),
    //uTheta_(readScalar(dict.lookup("uTheta"))),
    thetaR_(readScalar(dict.lookup("thetaR"))),
    maxAlphaHist(iF)
{
	if (dict.found("thetaAOW"))	thetaAOW_=(pTraits<scalar>(dict.lookup("thetaAOW")));
	else		thetaAOW_=thetaA_;

	if (dict.found("thetaROW"))	thetaROW_=(pTraits<scalar>(dict.lookup("thetaROW")));
	else		thetaROW_=thetaR_;

    evaluate();
}


Foam::hysteresisContactAngleFvPatchScalarField::
hysteresisContactAngleFvPatchScalarField
(
    const hysteresisContactAngleFvPatchScalarField& gcpsf
)
:
    alphaContactAngleFvPatchScalarField(gcpsf),
    thetaA_(gcpsf.thetaA_),
    //uTheta_(gcpsf.uTheta_),
    thetaR_(gcpsf.thetaR_),
    thetaAOW_(gcpsf.thetaAOW_),
    thetaROW_(gcpsf.thetaROW_),
    maxAlphaHist(gcpsf)
{}


Foam::hysteresisContactAngleFvPatchScalarField::
hysteresisContactAngleFvPatchScalarField
(
    const hysteresisContactAngleFvPatchScalarField& gcpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    alphaContactAngleFvPatchScalarField(gcpsf, iF),
    thetaA_(gcpsf.thetaA_),
    //uTheta_(gcpsf.uTheta_),
    thetaR_(gcpsf.thetaR_),
    thetaAOW_(gcpsf.thetaAOW_),
    thetaROW_(gcpsf.thetaROW_),
    maxAlphaHist(iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::scalarField>
Foam::hysteresisContactAngleFvPatchScalarField::theta
(
    const vectorField& Up,
    const vectorField & nHat,
    const vectorField & nWall
) const
{
    //if (uTheta_ < SMALL)
    //{
        //return tmp<scalarField>(new scalarField(size(), thetaA_));
    //}
    return maxAlphaHist*(0.5*(thetaROW_+thetaAOW_))+(1.0-maxAlphaHist)*(0.5*(thetaR_+thetaA_));

    //vectorField nf = patch().nf();

    // Calculated the component of the velocity parallel to the wall
    //vectorField Uwall = Up; //Up.patchInternalField() - Up;
    //Uwall -= (nWall & Uwall)*nWall;

    // Find the direction of the interface parallel to the wall
    //vectorField sWall = nHat - (nWall & nHat)*nWall;

    // Normalise sWall
    //sWall /= (mag(sWall) + SMALL);

    // Calculate Uwall resolved normal to the interface parallel to the interface
    //scalarField uwall = sWall & Uwall;

    //return thetaA_ + (thetaR_ - thetaR_);
}

void Foam::hysteresisContactAngleFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
	maxAlphaHist = max(*this, maxAlphaHist);

    alphaContactAngleFvPatchScalarField::evaluate();

}



void Foam::hysteresisContactAngleFvPatchScalarField::write(Ostream& os) const
{
    alphaContactAngleFvPatchScalarField::write(os);
    os.writeKeyword("thetaA") << thetaA_ << token::END_STATEMENT << nl;
    //os.writeKeyword("uTheta") << uTheta_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaR") << thetaR_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaAOW") << thetaAOW_ << token::END_STATEMENT << nl;
    os.writeKeyword("thetaROW") << thetaROW_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        hysteresisContactAngleFvPatchScalarField
    );
}


// ************************************************************************* //
