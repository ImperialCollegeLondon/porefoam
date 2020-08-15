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
#include "alphaContactAngleFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(alphaContactAngleFvPatchScalarField, 0);
}

bool Foam::alphaContactAngleFvPatchScalarField::reset = false;
//bool Foam::alphaContactAngleFvPatchScalarField::reset2 = false;

template<>
const char* Foam::NamedEnum
<
    Foam::alphaContactAngleFvPatchScalarField::limitControls,
    4
>::names[] =
{
    "none",
    "gradient",
    "zeroGradient",
    "alpha"
};

const Foam::NamedEnum
<
    Foam::alphaContactAngleFvPatchScalarField::limitControls,
    4
> Foam::alphaContactAngleFvPatchScalarField::limitControlNames_;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(p, iF),
    limit_(lcZeroGradient),
    transitionFactor_(0.9998)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& acpsf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(acpsf, p, iF, mapper),
    limit_(acpsf.limit_),
    transitionFactor_(acpsf.transitionFactor_)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF),
    limit_(limitControlNames_.read(dict.lookup("limit"))),
    transitionFactor_(0.9998)
{
    if (dict.found("gradient"))
    {
        gradient() = scalarField("gradient", dict, p.size());
        fixedGradientFvPatchScalarField::updateCoeffs();
        fixedGradientFvPatchScalarField::evaluate();
    }
    else
    {
        fvPatchField<scalar>::operator=(patchInternalField());
        gradient() = 0.0;
    }
    if (dict.found("transitionFactor"))
    {
			transitionFactor_=(pTraits<scalar>(dict.lookup("transitionFactor")));
         Info<<" transitionFactor:  " <<transitionFactor_<<endl;
   }
    else
    {
        	Warning<<" transitionFactor:  " <<transitionFactor_<<endl;
    }
}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& acpsf
)
:
    fixedGradientFvPatchScalarField(acpsf),
    limit_(acpsf.limit_),
    transitionFactor_(acpsf.transitionFactor_)
{}


Foam::alphaContactAngleFvPatchScalarField::alphaContactAngleFvPatchScalarField
(
    const alphaContactAngleFvPatchScalarField& acpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(acpsf, iF),
    limit_(acpsf.limit_),
    transitionFactor_(acpsf.transitionFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::alphaContactAngleFvPatchScalarField::evaluate
(
    const Pstream::commsTypes
)
{
	//if (reset2)
    //{
        //gradient() = -gradient();
        //return;
    //}
    
	if (reset)
    {
        gradient() = 0.0;
    }
    else if (limit_ == lcGradient)
    {
        gradient() =
        patch().deltaCoeffs()
       *(
           max(min
           (
               (1./transitionFactor_)*(*this + gradient()/patch().deltaCoeffs())-(0.5/transitionFactor_-0.5),
               scalar(1)), scalar(0)
           ) - *this
       );		
		//scalarField alphapintern=this->patchInternalField(); 
        //gradient() =
        //patch().deltaCoeffs()
       //*(
           //max(min
           //(
               //alphapintern + gradient()/patch().deltaCoeffs(),
               //scalar(1.0)), scalar(0.0)
           //) - alphapintern
       //);
    }
    else if (limit_ == lcZeroGradient)
    {
        gradient() = 0.0;
    }
    

    
    fixedGradientFvPatchScalarField::evaluate();

    if (limit_ == lcAlpha || limit_ == lcGradient)
    {
        scalarField::operator=(max(min( (1./transitionFactor_)*(*this)-(0.5/transitionFactor_-0.5) , scalar(1.0)), scalar(0.0)));
        //scalarField::operator=(max(min((*this), scalar(1.0)), scalar(0.0)));
    }
}


void Foam::alphaContactAngleFvPatchScalarField::write
(
    Ostream& os
) const
{
    fixedGradientFvPatchScalarField::write(os);
    os.writeKeyword("limit")
        << limitControlNames_[limit_] << token::END_STATEMENT << nl;
    os.writeKeyword("transitionFactor")
        << transitionFactor_ << token::END_STATEMENT << nl;
}


// ************************************************************************* //
