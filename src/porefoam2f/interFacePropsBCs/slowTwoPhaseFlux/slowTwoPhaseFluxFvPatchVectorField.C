/*-------------------------------------------------------------------------*\
 Copyright (C) 2010-2020  Ali Qaseminejad Raeini 

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>.
\*-------------------------------------------------------------------------*/

//! Description:
//!   Boundary condition for two-phase flow simulations, prescribing phase flow rates



#include "slowTwoPhaseFluxFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "primitivePatchInterpolation.H"
#include "foamTime.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

slowTwoPhaseFluxFvPatchVectorField::slowTwoPhaseFluxFvPatchVectorField
(
    const fvPatch& pd,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pd, iF),
    flowRate0_(0.0),
    flowRate1_(0.0),
    gradientFactor0_(0.0),
    gradientFactor1_(0.0),    
    pdFactor_(1.0),
    pcFactor_(0.0),
    curTimeIndex_(-1)
    //relaxationFactor_(0.1)
    
{
    //refValue() = *this;
    //refGrad() = vector::zero;
    //valueFraction() = 1.0;
}

slowTwoPhaseFluxFvPatchVectorField::slowTwoPhaseFluxFvPatchVectorField
(
    const slowTwoPhaseFluxFvPatchVectorField& ptf,
    const fvPatch& pd,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, pd, iF, mapper),
    flowRate0_(ptf.flowRate0_),
    flowRate1_(ptf.flowRate1_),
    gradientFactor0_(ptf.gradientFactor0_),
    gradientFactor1_(ptf.gradientFactor1_),    
    pdFactor_(ptf.pdFactor_),
    pcFactor_(ptf.pcFactor_),
    curTimeIndex_(-1)
    //  relaxationFactor_(ptf.relaxationFactor_)
{}


slowTwoPhaseFluxFvPatchVectorField::slowTwoPhaseFluxFvPatchVectorField
(
    const fvPatch& pd,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(pd, iF),
    flowRate0_(pTraits<scalar>(dict.lookup("flowRate0"))),
    flowRate1_(pTraits<scalar>(dict.lookup("flowRate1"))),
    gradientFactor0_(pTraits<scalar>(dict.lookup("gradientFactor0"))),
    gradientFactor1_(pTraits<scalar>(dict.lookup("gradientFactor1"))),
    //pdFactor_(0.9),
    //pcFactor_(2.0),
    pdFactor_(pTraits<scalar>(dict.lookupOrDefault("pdUniformization",0.9))),
    pcFactor_(pTraits<scalar>(dict.lookupOrDefault("pcUniformization",0.1))),
    curTimeIndex_(-1)
    //,    relaxationFactor_(0.1)    
    //phiName_(dict.lookupOrDefault<word>("phi", "phi")),
    //rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{



    if (dict.found("value"))
    {
        fixedValueFvPatchVectorField::operator==
        (
            Field<vector>("value", dict, pd.size())
        );
    }
    else
    {
        fvPatchField<vector>::operator==(vector(0,0,0));
    }
    
	//if (dict.found("pdUniformization"))
	//{
		//pdFactor_=(pTraits<scalar>(dict.lookup("pdUniformization")));
		//pcFactor_=(pTraits<scalar>(dict.lookup("pcUniformization")));
	//}
	//else if (dict.found("pdUniformizationFactor"))
	//{
		//pdFactor_=(pTraits<scalar>(dict.lookup("pdUniformizationFactor")));
		//pcFactor_=(pTraits<scalar>(dict.lookup("pcUniformizationFactor")));
		//Warning<<" Warning use of obselete keywords: pdUniformizationFactor / pcUniformizationFactor  " <<dict.name();
	//}
	//else
	//{
		//Warning<<" Warning pdUniformization / pcUniformization are not found in  " <<dict.name()<<endl;
		//Warning<<" Warning pdUniformization:  " <<pdFactor_<<endl;
		//Warning<<" Warning pcUniformization:  " <<pcFactor_<<endl;
	//}
		Info<<"  pdUniformization:  " <<pdFactor_<<endl;
		Info<<"  pcUniformization:  " <<pcFactor_<<endl;
}


slowTwoPhaseFluxFvPatchVectorField::slowTwoPhaseFluxFvPatchVectorField
(
    const slowTwoPhaseFluxFvPatchVectorField& pivpvf
)
:
    fixedValueFvPatchVectorField(pivpvf),
    flowRate0_(pivpvf.flowRate0_),
    flowRate1_(pivpvf.flowRate1_),
    gradientFactor0_(pivpvf.gradientFactor0_),
    gradientFactor1_(pivpvf.gradientFactor1_),
    pdFactor_(pivpvf.pdFactor_),
    pcFactor_(pivpvf.pcFactor_),
    curTimeIndex_(-1)
    //~, relaxationFactor_(pivpvf.relaxationFactor_)
{}


slowTwoPhaseFluxFvPatchVectorField::slowTwoPhaseFluxFvPatchVectorField
(
    const slowTwoPhaseFluxFvPatchVectorField& pivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(pivpvf, iF),
    flowRate0_(pivpvf.flowRate0_),
    flowRate1_(pivpvf.flowRate1_),
    gradientFactor0_(pivpvf.gradientFactor0_),
    gradientFactor1_(pivpvf.gradientFactor1_),
    pdFactor_(pivpvf.pdFactor_),
    pcFactor_(pivpvf.pcFactor_),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void slowTwoPhaseFluxFvPatchVectorField::updateCoeffs()
{
	if (updated())   {  return; }

	if (curTimeIndex_ == this->db().time().timeIndex()) return; 
	curTimeIndex_ = this->db().time().timeIndex();

	#define  curtailBADOFSET(a,b) ( min( max(a,b), (1.0-(b)) ) )

 
   const Field<scalar>& magS = patch().magSf();
	primitivePatchInterpolation pinterpolator(patch().patch());

	scalarField UBn =  0.0*mag(this->patchInternalField());

	Info().precision(4);

	const volScalarField& alpha1 =  db().lookupObject<volScalarField>("alpha1");
	const fvPatchField<scalar>& alpha1p =  patch().patchField<volScalarField, scalar>(alpha1);
	scalarField upif =  -this->patchInternalField() &  patch().Sf()/magS;
	upif = 0.99*upif+0.01*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(upif) );




	if (flowRate0_>1e-27)
	{
	 scalarField alphapSharp=curtailBADOFSET(10*(1.0-alpha1p.patchInternalField())-7.0,0.0);
	 scalar sumMagSa = gSum(magS*alphapSharp);
	 if (sumMagSa>1e-32)
	 {
		//scalar sign=(flowRate0_/sumMagS>0.0)*2.0-1.0;
		scalarField  phipn0=upif/(flowRate0_/sumMagSa);

		phipn0=phipn0*alphapSharp;
		

		phipn0=max(phipn0,1e-24);

		scalarField _logPatchValues = 0.95*log(phipn0)-0.05;//////////////   --init--   //////////////////////////////////////////

		const fvPatchField<scalar>&   pdp = patch().patchField<volScalarField, scalar>(db().lookupObject<volScalarField>("pd"));
		const fvPatchField<scalar>&   pcp = patch().patchField<volScalarField, scalar>(db().lookupObject<volScalarField>("pc"));


		Info<<" "<<patch().name()<<"0: alfa:"<<int(sumMagSa/gSum(magS)*100.0)/100.0<<"  u["<<gMin( phipn0 )<<" "<<gMax( phipn0 )<<"], ";	 
		Info<<" :P["<< gMax(pdp.patchInternalField())<< "  " <<gMin(pdp.patchInternalField()) 
		  <<"]  logU["<<int(gMin( _logPatchValues )*10.0)/10.0<<" "<<int(gMax(_logPatchValues)*10.0)/10.0<<"] ~ un:"<<gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa<<":  ";

		scalarField  pdpif = 0.9*(pdFactor_*pdp.patchInternalField() + pcFactor_*pcp.patchInternalField());
		scalarField alphaShSmooth = (pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(alphapSharp+1e-9) ));
		pdpif = 1.4*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(pdpif*alphapSharp+1e-9) )/alphaShSmooth-0.4*pdpif;
		pdpif = 1.4*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(pdpif*alphapSharp+1e-9) )/alphaShSmooth-0.4*pdpif;
		pdpif += 0.1*(pdFactor_*pdp.patchInternalField() + pcFactor_*pcp.patchInternalField());
		pdpif = alphapSharp*pdpif + (1.0-alphapSharp)*gSum(pdpif*alphapSharp)/gSum(alphapSharp);
		pdpif-= gMin(pdpif) - 10.0;//
		pdpif/=gMax(pdpif) + 20.0 ;//

		_logPatchValues +=  pdFactor_*(curtailBADOFSET(pdpif,0.2)-pdpif);

		scalar  uAvg=gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa +1e-8;
		Info<<"Pcor u="<<	 uAvg <<"; ";


		if (log(uAvg)>0.1)	 	   _logPatchValues += (0.0 - 0.5*(log(uAvg)) )            -0.09;	
		else if (log(uAvg)< -0.1)  _logPatchValues += (0.0 - 0.25*(log(uAvg)*(1.5-pdpif)) )+0.09;	/// mostly active	
		else					   _logPatchValues += (0.0 - 0.9*(log(uAvg)) );	 

		uAvg=gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa;
		Info<<"Qcor u="<<	 uAvg <<"; ";
		
		_logPatchValues = min(log(10.0)+min(log(uAvg+1e-1),0.0) , _logPatchValues);

		Info<<"cut u="<<gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa<<"\n";

		UBn -= flowRate0_/sumMagSa*alphapSharp*exp(_logPatchValues);
	 }
	}





	if (flowRate1_>1e-27)
	{
	 scalarField alphapSharp=curtailBADOFSET(10.0*(alpha1p.patchInternalField())-7.0,0.0);
	 scalar sumMagSa = gSum(magS*alphapSharp);
	 if (sumMagSa>1e-32)
	 {
		//scalar sign=(flowRate1_/sumMagS>0.0)*2.0-1.0;
		scalarField  phipn0=upif/(flowRate1_/sumMagSa);

		phipn0=phipn0*alphapSharp;
		

		phipn0=max(phipn0,1e-24);

		scalarField _logPatchValues = 0.95*log(phipn0)-0.05;//////////////   --init--   //////////////////////////////////////////

		const fvPatchField<scalar>&   pdp = patch().patchField<volScalarField, scalar>(db().lookupObject<volScalarField>("pd"));
		const fvPatchField<scalar>&   pcp = patch().patchField<volScalarField, scalar>(db().lookupObject<volScalarField>("pc"));



		Info<<" "<<patch().name()<<"1: alfa:"<<int(sumMagSa/gSum(magS)*100.0)/100.0<<"  u["<<gMin( phipn0 )<<" "<<gMax( phipn0 )<<"], ";	 
		Info<<" :P["<< gMax(pdp.patchInternalField())<< "  " <<gMin(pdp.patchInternalField()) 
		  <<"]  logU["<<int(gMin( _logPatchValues )*10.0)/10.0<<" "<<int(gMax(_logPatchValues)*10.0)/10.0<<"] ~ un:"<<gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa<<":  ";

		scalarField  pdpif = 0.9*(pdFactor_*pdp.patchInternalField() + pcFactor_*pcp.patchInternalField());
		scalarField alphaShSmooth = (pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(alphapSharp+1e-9) ));
		pdpif = 1.4*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(pdpif*alphapSharp+1e-9) )/alphaShSmooth-0.4*pdpif;
		pdpif = 1.4*pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(pdpif*alphapSharp+1e-9) )/alphaShSmooth-0.4*pdpif;
		pdpif += 0.1*(pdFactor_*pdp.patchInternalField() + pcFactor_*pcp.patchInternalField());
		pdpif = alphapSharp*pdpif + (1.0-alphapSharp)*gSum(pdpif*alphapSharp)/gSum(alphapSharp);
		pdpif-= gMin(pdpif) - 10.0;//
		pdpif/=gMax(pdpif) + 20.0 ;//

		_logPatchValues +=  pdFactor_*(curtailBADOFSET(pdpif,0.2)-pdpif);

		scalar  uAvg=gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa +1e-8;
		Info<<"Pcor u="<<	 uAvg <<"; ";


		if (log(uAvg)>0.1)	 	   _logPatchValues += (0.0 - 0.5*(log(uAvg)) )            -0.09;	
		else if (log(uAvg)< -0.1)  _logPatchValues += (0.0 - 0.25*(log(uAvg)*(1.5-pdpif)) )+0.09;	/// mostly active	
		else					   _logPatchValues += (0.0 - 0.9*(log(uAvg)) );	 

		uAvg=gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa ;
		Info<<"Qcor u="<<	 uAvg <<"; ";
		
		_logPatchValues = min(log(10.0)+min(log(uAvg+1e-1),0.0) , _logPatchValues);

		Info<<"cut u="<<gSum(alphapSharp*exp(_logPatchValues)*magS)/sumMagSa<<"\n";

		UBn -= flowRate1_/sumMagSa*alphapSharp*exp(_logPatchValues);
	 }
	}








	*this == UBn*patch().nf();

   fixedValueFvPatchVectorField::updateCoeffs();
}


void slowTwoPhaseFluxFvPatchVectorField::
write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("flowRate0")
        << flowRate0_ << token::END_STATEMENT << nl;
    os.writeKeyword("flowRate1")
        << flowRate1_ << token::END_STATEMENT << nl;
    os.writeKeyword("gradientFactor0")
        << gradientFactor0_ << token::END_STATEMENT << nl;
    os.writeKeyword("gradientFactor1")
        << gradientFactor1_ << token::END_STATEMENT << nl;
    os.writeKeyword("pdUniformization")
        << pdFactor_ << token::END_STATEMENT << nl;        
    os.writeKeyword("pcUniformization")
        << pcFactor_ << token::END_STATEMENT << nl;        
    writeEntry("value", os); 
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    slowTwoPhaseFluxFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
