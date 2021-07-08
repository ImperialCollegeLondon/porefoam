/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is partly based on OpenFOAM code.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

\*---------------------------------------------------------------------------*/

#include "mathematicalConstants.H"
#include "volFields.H"
#include "surfaceFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "primitivePatchInterpolation.H"
#undef NoRepository
#include "fixedRelaxedMeanValueFvPatchField.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
fixedRelaxedMeanValueFvPatchField<Type>::fixedRelaxedMeanValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    meanValue1_(pTraits<Type>::zero),
    meanValue2_(pTraits<Type>::zero),
    curTimeIndex_(-1),
    relaxationFactor_(0.05),
    thicknessFactor_(0.1)
{}


template<class Type>
fixedRelaxedMeanValueFvPatchField<Type>::fixedRelaxedMeanValueFvPatchField
(
    const fixedRelaxedMeanValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    meanValue1_(ptf.meanValue1_),
    meanValue2_(ptf.meanValue2_),
    curTimeIndex_(-1),
    relaxationFactor_(0.05),
    thicknessFactor_(0.1)
{}


template<class Type>
fixedRelaxedMeanValueFvPatchField<Type>::fixedRelaxedMeanValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    meanValue1_(pTraits<Type>(dict.lookupOrDefault("meanValue1",pTraits<Type>(dict.lookupOrDefault("meanValue",pTraits<Type>::zero))))),
    meanValue2_(pTraits<Type>(dict.lookupOrDefault("meanValue1",pTraits<Type>(dict.lookupOrDefault("meanValue",pTraits<Type>::zero))))),
    curTimeIndex_(-1),
    relaxationFactor_(0.05),
    thicknessFactor_(dict.lookupOrDefault("thicknessFactor",scalar(0.1)))
{
    if (dict.found("value"))
    {
        fixedValueFvPatchField<Type>::operator==
        (
            Field<Type>("value", dict, p.size())
        );
    }
    else
    {
        fixedValueFvPatchField<Type>::operator==(0.5*meanValue1_+0.5*meanValue2_);
    }
    
    if (dict.found("relax"))
    {
		relaxationFactor_=readScalar(dict.lookup("relax"));
		Info<<" P BC relaxation factor used is  "<< relaxationFactor_<<endl ;
    }
    else    if (dict.found("relaxationFactor"))
    {
		relaxationFactor_=readScalar(dict.lookup("relaxationFactor"));
		Info<<" P BC relaxationFactor used is  "<< relaxationFactor_<<endl ;

    }
    else
    {
		Warning<<" warning relaxationFactor is not found in  " <<dict.name();
		Info<<" P BC relaxationFactor used is  "<< relaxationFactor_<<endl ;

	}
	
		Info<<" P BC mean values: "<< meanValue1_ <<"  "<<meanValue2_<<endl ;
		Info<<" P BC thicknessFactor: "<< thicknessFactor_ <<endl ;

}


template<class Type>
fixedRelaxedMeanValueFvPatchField<Type>::fixedRelaxedMeanValueFvPatchField
(
    const fixedRelaxedMeanValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    meanValue1_(ptf.meanValue1_),
    meanValue2_(ptf.meanValue2_),
    curTimeIndex_(ptf.curTimeIndex_),
    relaxationFactor_(ptf.relaxationFactor_),
    thicknessFactor_(ptf.thicknessFactor_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
template<class Type>
void fixedRelaxedMeanValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    Field<Type>::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
template<class Type>
void fixedRelaxedMeanValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);
}

#define curtailBADOFSET(a,b) (min (max(a,b),(1.-(b)))    )

// Update the coefficients associated with the patch field
template<class Type>
void fixedRelaxedMeanValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }
		//Info<<" P BC relaxationFactor used is  "<< relaxationFactor_ ;
	Field<Type>& patchField = *this;
	const Field<scalar> pMagPhi = mag(this->patch().template lookupPatchField<surfaceScalarField, scalar>("phi"))+1e-24;

	const Field<scalar> alphap = curtailBADOFSET((1./thicknessFactor_)* (this->patch().template lookupPatchField<volScalarField, scalar>("alpha1")) -(0.5/thicknessFactor_-0.5),0.); 

	//primitivePatchInterpolation pinterpolator(patch().patch());
	//alphap = 1.4*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(alphap) )-0.4*alphap;

//     if (curTimeIndex_ != this->db().time().timeIndex())  
    {
        Field<Type> pif = this->patchInternalField();

        patchField =  (alphap*meanValue1_+(1.-alphap)*meanValue2_)  + relaxationFactor_*(pif-gSum(pif*pMagPhi)/(gSum(pMagPhi))) ;
        patchField == (alphap*meanValue1_+(1.-alphap)*meanValue2_)  + relaxationFactor_*(pif-gSum(pif*pMagPhi)/(gSum(pMagPhi))) ;
        //patchField = pif  + relaxationFactor_*(meanValue_-gSum(pif*pMagPhi)/(gSum(pMagPhi))) ;
        //patchField == pif  + relaxationFactor_*(meanValue_-gSum(pif*pMagPhi)/(gSum(pMagPhi))) ;

        curTimeIndex_ = this->db().time().timeIndex();
    }

    fixedValueFvPatchField<Type>::updateCoeffs();
}


// Write
template<class Type>
void fixedRelaxedMeanValueFvPatchField<Type>::write(Ostream& os) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("meanValue1")
        << meanValue1_ << token::END_STATEMENT << nl;
    os.writeKeyword("meanValue2")
        << meanValue2_ << token::END_STATEMENT << nl;
    os.writeKeyword("relaxationFactor")
        << relaxationFactor_ << token::END_STATEMENT << nl;        
    os.writeKeyword("thicknessFactor")
        << thicknessFactor_ << token::END_STATEMENT << nl;        
    this->writeEntry("value", os); 
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
