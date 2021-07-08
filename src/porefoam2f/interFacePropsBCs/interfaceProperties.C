/*-------------------------------------------------------------------------*\
 Interface force and capillary pressure computation library

 Copyright (C) 2014-2020  Mosayeb Shams
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


#include "interfaceProperties.H"
#include "alphaContactAngleFvPatchScalarField.H"
//#include "sharpInOutAlphaFvPatchScalarField.H"
#include "mathematicalConstants.H"
#include "primitivePatchInterpolation.H"
#include "surfaceInterpolate.H"
#include "fvcDiv.H"
#include "fvcGrad.H"
#include "fvcSnGrad.H"
#include "fvcAverage.H"
#include "fvmLaplacian.H"
#include "fvCFD.H"
#include "OFstream.H"
#include "syncTools.H"

//#include "alphaf.H" // TODO this needs a lot of work to be usable

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define curtailBADOFSET(a,b) ( min( max(a,b), (1.-(b)) ) )



// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::mathematicalConstant::pi/180.;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
#ifndef interfaceProperties_interfaceProperties // folding hint
    alpha1_(alpha1), // mesh() requires this to be initialized first
    alpha1S_(IOobject( "alpha1S",  timeName(),  mesh()), alpha1),

    U_(U),

    pc_(IOobject("pc", timeName(), mesh(),  IOobject::MUST_READ, IOobject::AUTO_WRITE), mesh() ),

    deltaS_(IOobject( "deltaS", timeName(), mesh()), fvc::snGrad(alpha1, "uncorrected")),
    alphaSh_(IOobject("alphaSh", timeName(), mesh()),  mesh(),  dimensionedScalar("alphaSh", dimless, 0.),  pc_.boundaryField().types()),    

    //K_
    //(
        //IOobject( "K",  timeName(),  mesh()   ),
        //mesh(),        dimensionedScalar("K", dimless/dimLength, 0.),
        //pc_.boundaryField().types()
    //),

    nHatf_
    (
        IOobject(  "nHatf",  timeName(),  mesh()   ),
        mesh(),        dimensionedScalar("nHatf", dimArea, 0.)
    ),
    gPc_
    (
        IOobject( "gPc",  timeName(),  mesh(),  IOobject::READ_IF_PRESENT   ),//,  IOobject::AUTO_WRITE
        mesh(),   dimensionedVector("gPc", dimPressure/dimLength, vector(0.,0.,0.)),
        pc_.boundaryField().types()
    ),
    sgPc_
    (
        IOobject( "sgPc",  timeName(),  mesh(),  IOobject::READ_IF_PRESENT   ), //! NO_WRITE is not the most accurate when restarting
        linearInterpolate(gPc_) & mesh().Sf()
    ),
    alpha1f_
    (
        IOobject( "alpha1f_",  timeName(),  mesh(),  IOobject::READ_IF_PRESENT   ), //! NO_WRITE is not the most accurate when restarting
        linearInterpolate(alpha1_)
    ),


    //nS_
    //(
        //IOobject( "nS",  timeName(),  mesh() ),
        //mesh(),       dimensionedVector("nS", dimless, vector(0.,0.,0.))//,
    //),
    nw_
    (
        IOobject( "nw",  timeName(),  mesh() ),
        mesh(),     dimensionedVector("nw", dimless, vector(0.,0.,0.))
    ),
    
    edgemarks_(mesh().edges().size(),0),
    
    Internalfaces1_
    (
        IOobject( "Internalfaces1",  timeName(), mesh() ),
        mesh(),  dimensionedScalar("Internalfaces1", dimless, 1.)
    ),
    BInternalfs_
    (
        IOobject( "BInternalfs",  timeName(), mesh() ),
        mesh(),  dimensionedScalar("BInternalfs", dimless, 0.)
    ),
    AvgInternFaces1_
    (
        IOobject( "AvgInternFaces1",  timeName(), mesh() ),
        mesh(),  dimensionedScalar("AvgInternFaces1", dimless, 1.)
    ),
    IsRefCandid_(mesh().cells().size(),1),

    sgPcErr_
    (  IOobject( "sgPce", timeName(), mesh(), IOobject::READ_IF_PRESENT ), //! NO_WRITE means restarting releases the filters
       mesh(),    dimensionedScalar("sgPce", dimPressure/dimLength*dimArea, 0.)
    )  ,

    sgPcErrn_
    ( IOobject( "sgPcen", timeName(), mesh() ),
      mesh(),
      dimensionedScalar("sgPcen", dimPressure/dimLength*dimArea, 0.)
    ),

  pMesh_(mesh()),
  vpi_(mesh()),
                   	 //pvi_(pMesh_,mesh())
	//interfPointsOld_(mesh().points().size(),0),
  distPointInterface_
  ( IOobject( "distPointInterface", timeName(), mesh()  ,  IOobject::NO_READ ),//!  WRITE for visualization
	 pMesh_,	 dimensionedVector("distPointInterface",dimLength,vector::zero)
  ),


    //transportPropertiesDict_(dict),
    pcThicknessFactor_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("pcThicknessFactor") ) ),
    fcCorrectTangent_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("fcCorrectTangent") ) ),
    fcCorrectTangentRelax_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("fcCorrectTangentRelax") ) ),
    fcdFilter_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("fcdFilter") ) ),
    nPcNonOrthCorr_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("nPcNonOrthCorr") ) ),
    pcRefCellOrig_( readLabel( mesh().solutionDict().subDict("PIMPLE").lookup("pcRefCell") ) ),
    pcRefValueOrig_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("pcRefValue") ) ),
    smoothingKernel_( readLabel( mesh().solutionDict().subDict("PIMPLE").lookup("smoothingKernel") ) ),
    smoothingRelaxFactor_( readScalar( mesh().solutionDict().subDict("PIMPLE").lookup("smoothingRelaxFactor") ) ),
    wallSmoothingKernel_( readLabel( mesh().solutionDict().subDict("PIMPLE").lookup("wallSmoothingKernel") ) ),
    sigma_(dict.lookup("sigma")),
    deltaN_  ("deltaN", 1e-12/pow(average(mesh().V()), 1./3.) )
#endif

{
    const fvMesh& msh = mesh();

	Info <<"setting up filtered surface force model"<<endl;

	//SET_REF_(pcRefCell,pcRefValue,Orig_,alpha1_);

    //setRefCell(pc_, msh.solutionDict().subDict("PIMPLE"), pcRefCell, pcRefValue);

    //Info<< "cAlpha_" <<cAlpha_<<endl;
    Info<< "pcThicknessFactor_" <<pcThicknessFactor_<<endl;
    Info<< "fcCorrectTangent_" <<fcCorrectTangent_<<endl;
    Info<< "fcCorrectTangentRelax_" <<fcCorrectTangentRelax_<<endl;
    Info<< "fcdFilter_" <<fcdFilter_<<endl;
    Info<< "nPcNonOrthCorr_" <<nPcNonOrthCorr_<<endl;
    Info<< "pcRefCellOrig_" <<pcRefCellOrig_<<endl;
    Info<< "pcRefValueOrig_" <<pcRefValueOrig_<<endl;
    Info<< "smoothingKernel_" <<smoothingKernel_<<endl;
    Info<< "smoothingRelaxFactor_" <<smoothingRelaxFactor_<<endl;
    Info<< "wallSmoothingKernel_" <<wallSmoothingKernel_<<endl;

    Info<< "alpha1S_.boundaryField().types():" <<alpha1S_.boundaryField().types()<<endl;
    Info<< "alpha1_.boundaryField().types():" <<alpha1_.boundaryField().types()<<endl;

	const fvBoundaryMesh& boundary = msh.boundary();

	{ ///. nw

		forAll(boundary, bi)
		{
			if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[bi]))
			{
				nw_.boundaryField()[bi]=boundary[bi].nf();
				nw_.boundaryField()[bi]==boundary[bi].nf();//tomakesure
			}
		}	


		forAll(boundary, bi)
		{
			if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[bi]))
			{
				primitivePatchInterpolation pinterpolator(msh.boundaryMesh()[bi]);
				for (int i=0;i<wallSmoothingKernel_;i++)
				{
					nw_.boundaryField()[bi]==
					pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(
					nw_.boundaryField()[bi]));
				}
			}
			if (isA<alphaContactAngleFvPatchScalarField>(alpha1S_.boundaryField()[bi]))
					 const_cast<alphaContactAngleFvPatchScalarField&>
					 ( refCast<const alphaContactAngleFvPatchScalarField>  (alpha1S_.boundaryField()[bi]) ).noLimit();
			//if (isA<sharpInOutAlphaFvPatchScalarField>(alpha1S_.boundaryField()[bi]))
			//{
					 //(const_cast<sharpInOutAlphaFvPatchScalarField&>
					 //( refCast<const sharpInOutAlphaFvPatchScalarField>  (alpha1S_.boundaryField()[bi]) )).noLimit();
					//Info<<"\n\n\n sdffdfdf\\n\n\n "<<endl;
			//}
		}

		nw_.internalField() *= 0.;
		nw_.internalField() = fvc::interpolate(fvc::average(nw_))->internalField();		
		nw_ /=  mag(nw_) + 1e-15;



	}







{// collect faces neighbour to boundary patches boundaryInternalFaces_
	const labelListList & faceEdges=msh.faceEdges();
	//const edgeList & edges=msh.edges();
	const fvBoundaryMesh& patches = msh.boundary(); 
	forAll(patches, bI) if (!patches[bI].coupled())
	{
		Internalfaces1_.boundaryField()[bI] ==0.000000001;
		const scalarField & CApfs = alpha1_.boundaryField()[bI];
		forAll(CApfs, pfI)
		{
			const labelList& fes = faceEdges[patches[bI].patch().start()+pfI];
			forAll(fes, eI)  edgemarks_[fes[eI]] = 1.; 
		}
	}
	AvgInternFaces1_.internalField()=fvc::average(Internalfaces1_);

	//BInternalfs_=0.*nHatf_;
	forAll(BInternalfs_, fI)
	{
		const labelList& fes = faceEdges[fI];	  
		forAll(fes, eI)
		{
			BInternalfs_[fI] = max( BInternalfs_[fI], 1.*edgemarks_[fes[eI]] );
		}
	}

	forAll(patches, bI) 
		for(auto cI:patches[bI].faceCells()) IsRefCandid_[cI]=0;
	//const auto& neis=mesh().neighbour();
	//const auto& owns=mesh().owner();
	//forAll(neis,fI) {
		//if     (!IsRefCandid_[owns[fI]]) IsRefCandid_[neis[fI]]=0;
		//else if(!IsRefCandid_[neis[fI]]) IsRefCandid_[owns[fI]]=0;
	//}
	const auto& Vols=mesh().V();
	scalar Vavg=average(mesh().V()).value();
	forAll(Vols,cI) if(Vols[cI]<Vavg*0.99) IsRefCandid_[cI]=0;
	Info<<"avg(IsRefCandid): "<<gSum(IsRefCandid_)<<" / "<<returnReduce<scalar>(scalar(IsRefCandid_.size()), sumOp<scalar>())<<endl;

	label bIFssize=0;
	forAll(BInternalfs_, fI)
		if (BInternalfs_[fI]>0.5) ++bIFssize;

	boundaryInternalFaces_.resize(bIFssize);
	label bIFsi=-1;
	forAll(BInternalfs_, fI)
		if (BInternalfs_[fI]>0.5) boundaryInternalFaces_[++bIFsi]=fI;

}


/* { /// correction to handle zero-gradient boundary condition in smoothing for filtering

	volScalarField corrS1=fvc::average(fvc::interpolate(boundaryCorr_,"limitedScheme"));

	forAll(boundary, bI)
    {
        if (!boundary[bI].coupled())
        {
				boundaryCorr_.boundaryField()[bI] = 0.;
				boundaryCorr_.boundaryField()[bI] == 0.;
        }
    }
	volScalarField corrS2=fvc::average(fvc::interpolate(boundaryCorr_,"limitedScheme"));
    forAll(boundary, bI)
    {
        if (!boundary[bI].coupled())
        {
				
				boundaryCorr_.boundaryField()[bI] =

					( (corrS1.boundaryField()[bI].patchInternalField())-
					  (corrS2.boundaryField()[bI].patchInternalField()) )
					/(boundaryCorr_.boundaryField()[bI].patchInternalField());
				
				
				boundaryCorr_.boundaryField()[bI] ==
					( (corrS1.boundaryField()[bI].patchInternalField())-
					  (corrS2.boundaryField()[bI].patchInternalField()) )
					/(boundaryCorr_.boundaryField()[bI].patchInternalField());
        }
    }
    boundaryCorr_.internalField()=0;
    
    
    

	
}

*/



{ // not being used
		//const fvMesh& msh = mesh();
		//const volVectorField& Cc = msh.C();
		//const surfaceVectorField& Cf = msh.Cf();
		//const labelList& owner = msh.faceOwner();
		//const labelList& neighbour = msh.faceNeighbour();
		//const fvBoundaryMesh& patches = msh.boundary();
		//forAll(neighbour, facei)
			//deltaCC_[facei] = Cc[neighbour[facei]] - Cc[owner[facei]];
		//forAll(patches, bI)
		 //if(patches[bI].size())
		 //{
		    //const vectorField&  pCf = Cf.boundaryField()[bI];;
		    //vectorField&  pdeltaCC_ = deltaCC_.boundaryField()[bI];
		    //label pStart = patches[bI].patch().start();
	         //forAll(pCf, pfI)
					//pdeltaCC_[pfI] += pCf[pfI] - Cc[owner[pStart+pfI]];

		  //if (deltaCC_.boundaryField()[bI].coupled())
          //{
			 //const vectorField  pCcN = Cc.boundaryField()[bI].patchNeighbourField();
            //forAll(pCcN, pfI)
				//pdeltaCC_[pfI] +=  pCcN[pfI] - pCf[pfI];
          //}
		 //}
}


    //correct(fvc::snGrad(alpha1));
}


// ************************************************************************* //




// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::interfaceProperties::correct(scalar relaxDelS)
{

	const dictionary pimple =mesh().solutionDict().subDict("PIMPLE");
    //cAlpha_=   readScalar( pimple.lookup("cAlpha") );
    pcThicknessFactor_=   readScalar( pimple.lookup("pcThicknessFactor") );
    word surfaceForceModel(pimple.lookup("surfaceForceModel"));
    word alphafModel(pimple.lookup("alphafModel"));




	scalar stime=alpha1_.time( ).elapsedCpuTime();


    const fvMesh& msh = mesh();
    const surfaceVectorField& Sf = msh.Sf();
    const fvBoundaryMesh& boundary = msh.boundary(); 



	alphaSh_ == curtailBADOFSET((1./pcThicknessFactor_)*alpha1_-(0.5/pcThicknessFactor_-0.5),0.); alphaSh_.correctBoundaryConditions();
	if (relaxDelS<0.0001) deltaS_= fvc::snGrad(alphaSh_);
	else                  deltaS_= 0.5*(deltaS_+fvc::snGrad(alphaSh_));


	if(alphafModel=="advect")
		alpha1f_ = linearInterpolate(min( max(alpha1_,0.), 1. ));
	else
		alpha1f_ = linearInterpolate(min( max(alpha1_,0.), 1. )); // minmax is added hoping to solve solver crash on non-orthogonal meshes


	volScalarField alpha1Stmp=curtailBADOFSET(alpha1_*1.02-0.01,1e-3);

	volScalarField a1a2 = 0.95+0.1*sqrt(alpha1Stmp*(1.-alpha1Stmp));
	alpha1Stmp = a1a2*alpha1Stmp+(1.-a1a2)*fvc::average(alpha1f_*Internalfaces1_)/fvc::average(Internalfaces1_);
	alpha1Stmp.correctBoundaryConditions();

	a1a2 = 2.*sqrt(alpha1Stmp*(1.-alpha1Stmp));

	alpha1Stmp = a1a2*alpha1Stmp + (1.-a1a2)*fvc::average(fvc::interpolate(alpha1Stmp)*Internalfaces1_)/fvc::average(Internalfaces1_); 
	alpha1Stmp.correctBoundaryConditions();

	alpha1Stmp = 0.99*alpha1Stmp + 0.01*fvc::average(fvc::interpolate(alpha1Stmp)*Internalfaces1_)/fvc::average(Internalfaces1_);


	alpha1S_ = atanh(1.8*alpha1Stmp-0.9); //- SYN123433
	alpha1S_.correctBoundaryConditions();





	surfaceScalarField  stf
	(	IOobject( "stf", timeName(), msh ),
		msh,  dimensionedScalar("stf", dimPressure/dimLength*dimArea, 0.)
	);


	if (surfaceForceModel=="FSF")
	{
		Info<<surfaceForceModel<<endl;
		calcCurvatureFSF(stf,deltaS_,a1a2);
	}
	else
	{
		stf=sigma_*calcCurvatureFConservative(); 
	}


/// //////////////////////////////////////////////////////////////////////////////////////////

	#include "separatePc_filter.H"




}

