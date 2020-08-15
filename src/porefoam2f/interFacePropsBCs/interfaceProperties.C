/*-------------------------------------------------------------------------*\
This code is part of poreFOAM, a suite of codes written using OpenFOAM
for direct simulation of flow at the pore scale. 	
You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.


The code has been developed by Ali Qaseminejad Raeini as a part his PhD 
at Imperial College London, under the supervision of Branko Bijeljic 
and Martin Blunt. 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Mosayeb Shams:    m.shams14@imperial.ac.uk
Ali Q Raeini:     a.q.raeini@imperial.ac.uk
Branko Bijeljic:  b.bijeljic@imperial.ac.uk
Martin J Blunt:   m.blunt@imperial.ac.uk

Description:
	interface force and capillary pressure computation

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

#define curtailBADOFSET(a,b) (min (max(a,b),(1.0-(b)))    )



// * * * * * * * * * * * * * * * Static Member Data  * * * * * * * * * * * * //

const Foam::scalar Foam::interfaceProperties::convertToRad =
    Foam::mathematicalConstant::pi/180.0;


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::interfaceProperties::interfaceProperties
(
    const volScalarField& alpha1,
    const volVectorField& U,
    const IOdictionary& dict
)
:
#ifndef interfaceProperties_interfaceProperties
    transportPropertiesDict_(dict),
    pcThicknessFactor_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcThicknessFactor") ) ),
    fcCorrectTangent_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("fcCorrectTangent") ) ),
    fcCorrectTangentRelax_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("fcCorrectTangentRelax") ) ),
    fcdFilter_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("fcdFilter") ) ),
    nPcNonOrthCorr_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("nPcNonOrthCorr") ) ),
    pcRefCellOrig_( readLabel( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefCell") ) ),
    pcRefValueOrig_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("pcRefValue") ) ),
    smoothingKernel_( readLabel( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("smoothingKernel") ) ),
    smoothingRelaxFactor_( readScalar( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("smoothingRelaxFactor") ) ),
    wallSmoothingKernel_( readLabel( alpha1.mesh().solutionDict().subDict("PIMPLE").lookup("wallSmoothingKernel") ) ),
    sigma_(dict.lookup("sigma")),
    deltaN_  ("deltaN", 1.0e-12/pow(average(alpha1.mesh().V()), 1.0/3.0) ),
    alpha1_(alpha1),
    alpha1S_
    (
		IOobject( "alpha1S",  alpha1_.time().timeName(),  alpha1_.mesh(), IOobject::NO_READ, IOobject::NO_WRITE   ),
		alpha1
    ),

    U_(U),

    pc_
    (
        IOobject( "pc",  alpha1_.time().timeName(),  alpha1_.mesh(),  IOobject::MUST_READ,  IOobject::AUTO_WRITE   ),
        alpha1_.mesh()
    ),

    deltaS_
    (
        IOobject( "deltaS",  alpha1_.time().timeName(),  alpha1_.mesh()   ), 
        fvc::snGrad(alpha1, "uncorrected")
    ),
    alphaSh_
    (
        IOobject( "alphaSh",  alpha1_.time().timeName(),  alpha1_.mesh()   ),
        alpha1_.mesh(),        dimensionedScalar("alphaSh", dimless, 0.0),
        pc_.boundaryField().types()
    ),    

    //K_
    //(
        //IOobject( "K",  alpha1_.time().timeName(),  alpha1_.mesh()   ),
        //alpha1_.mesh(),        dimensionedScalar("K", dimless/dimLength, 0.0),
        //pc_.boundaryField().types()
    //),

    nHatf_
    (
        IOobject(  "nHatf",  alpha1_.time().timeName(),  alpha1_.mesh()   ),
        alpha1_.mesh(),        dimensionedScalar("nHatf", dimArea, 0.0)
    ),
    gPc_
    (
        IOobject( "gPc",  alpha1_.time().timeName(),  alpha1_.mesh(),  IOobject::READ_IF_PRESENT   ),//,  IOobject::AUTO_WRITE
        alpha1_.mesh(),   dimensionedVector("gPc", dimPressure/dimLength, vector(0.0,0.0,0.0)),
        pc_.boundaryField().types()
    ),
    sgPc_
    (
        IOobject( "sgPc",  alpha1_.time().timeName(),  alpha1_.mesh(),  IOobject::READ_IF_PRESENT   ), //! NO_WRITE is not the most accurate when restarting
        linearInterpolate(gPc_) & alpha1_.mesh().Sf()
    ),
    alpha1f_
    (
        IOobject( "alpha1f_",  alpha1_.time().timeName(),  alpha1_.mesh(),  IOobject::READ_IF_PRESENT   ), //! NO_WRITE is not the most accurate when restarting
        linearInterpolate(alpha1_)
    ),


    //nS_
    //(
        //IOobject( "nS",  alpha1_.time().timeName(),  alpha1_.mesh() ),
        //alpha1_.mesh(),       dimensionedVector("nS", dimless, vector(0.0,0.0,0.0))//,
    //),
    nw_
    (
        IOobject( "nw",  alpha1_.time().timeName(),  alpha1_.mesh() ),
        alpha1_.mesh(),     dimensionedVector("nw", dimless, vector(0.0,0.0,0.0))
    ),
    
    edgemarks_(alpha1_.mesh().edges().size(),0),
    
    Internalfaces1_
    (
        IOobject( "Internalfaces1",  alpha1_.time().timeName(), alpha1_.mesh() ),
        alpha1_.mesh(),  dimensionedScalar("Internalfaces1", dimless, 1.0)
    ),
    BInternalfs_
    (
        IOobject( "BInternalfs",  alpha1_.time().timeName(), alpha1_.mesh() ),
        alpha1_.mesh(),  dimensionedScalar("BInternalfs", dimless, 0.0)
    ),
    AvgInternFaces1_
    (
        IOobject( "AvgInternFaces1",  alpha1_.time().timeName(), alpha1_.mesh() ),
        alpha1_.mesh(),  dimensionedScalar("AvgInternFaces1", dimless, 1.0)
    ),

    sgPcErr_
    (  IOobject( "sgPce", alpha1_.time().timeName(), alpha1_.mesh(), IOobject::READ_IF_PRESENT ), //! NO_WRITE means restarting releases the filters
       alpha1_.mesh(),    dimensionedScalar("sgPce", dimPressure/dimLength*dimArea, 0.0)
    )  ,

    sgPcErrn_
    ( IOobject( "sgPcen", alpha1_.time().timeName(), alpha1_.mesh() ),
      alpha1_.mesh(),
      dimensionedScalar("sgPcen", dimPressure/dimLength*dimArea, 0.0)
    ),

  pMesh_(alpha1.mesh()),
  vpi_(alpha1.mesh()),
                   	 //pvi_(pMesh_,alpha1_.mesh())
	//interfPointsOld_(alpha1.mesh().points().size(),0),
  distPointInterface_
  ( IOobject( "distPointInterface", alpha1_.time().timeName(), alpha1_.mesh()  ,  IOobject::NO_READ ),//!  WRITE for visualization
	 pMesh_,	 dimensionedVector("distPointInterface",dimLength,vector::zero)
  )


#endif

{
    const fvMesh& mesh = alpha1_.mesh();

	Info <<"setting up filtered surface force model"<<endl;

	//SET_REF_(pcRefCell,pcRefValue,Orig_,alpha1_);

    //setRefCell(pc_, mesh.solutionDict().subDict("PIMPLE"), pcRefCell, pcRefValue);

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

	const fvBoundaryMesh& boundary = alpha1_.mesh().boundary();

	{ ///. nw

		forAll(boundary, patchi)
		{
			if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[patchi]))
			{
				nw_.boundaryField()[patchi]=boundary[patchi].nf();
				nw_.boundaryField()[patchi]==boundary[patchi].nf();//tomakesure
			}
		}	


		forAll(boundary, patchi)
		{
			if (isA<alphaContactAngleFvPatchScalarField>(alpha1_.boundaryField()[patchi]))
			{
				primitivePatchInterpolation pinterpolator(alpha1_.mesh().boundaryMesh()[patchi]);
				for (int i=0;i<wallSmoothingKernel_;i++)
				{
					nw_.boundaryField()[patchi]==
					pinterpolator.pointToFaceInterpolate(pinterpolator.faceToPointInterpolate(
					nw_.boundaryField()[patchi]));
				}
			}
			if (isA<alphaContactAngleFvPatchScalarField>(alpha1S_.boundaryField()[patchi]))
					 const_cast<alphaContactAngleFvPatchScalarField&>
					 ( refCast<const alphaContactAngleFvPatchScalarField>  (alpha1S_.boundaryField()[patchi]) ).noLimit();
			//if (isA<sharpInOutAlphaFvPatchScalarField>(alpha1S_.boundaryField()[patchi]))
			//{
					 //(const_cast<sharpInOutAlphaFvPatchScalarField&>
					 //( refCast<const sharpInOutAlphaFvPatchScalarField>  (alpha1S_.boundaryField()[patchi]) )).noLimit();
					//Info<<"\n\n\n sdffdfdf\\n\n\n "<<endl;
			//}

				
		}

		nw_.internalField()*=0.0;
		nw_.internalField()=
		fvc::interpolate(fvc::average(nw_))->internalField();		
		nw_ /=  mag(nw_) + 1.0e-15;



	}







{///. collect faces neibour to boundary patches boundaryInternalFaces_
	const labelListList & faceEdges=mesh.faceEdges();
	//const edgeList & edges=mesh.edges();
	const fvBoundaryMesh& patches = mesh.boundary(); 
	forAll(patches, patchI) if (!patches[patchI].coupled())
	{
		Internalfaces1_.boundaryField()[patchI] ==0.000000001;
		const scalarField & CApfs = alpha1_.boundaryField()[patchI];
		forAll(CApfs, pfI)
		{
			const labelList& fes = faceEdges[patches[patchI].patch().start()+pfI];
			forAll(fes, eI)
			{
				edgemarks_[fes[eI]] = 1.0;
			}
		}
	}
	AvgInternFaces1_.internalField()=fvc::average(Internalfaces1_);

	//BInternalfs_=0.0*nHatf_;
	forAll(BInternalfs_, fI)
	{
		const labelList& fes = faceEdges[fI];	  
		forAll(fes, eI)
		{
			BInternalfs_[fI] = max( BInternalfs_[fI], 1.0*edgemarks_[fes[eI]] );
		}
	}
	
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

	forAll(boundary, patchI)
    {
        if (!boundary[patchI].coupled())
        {
				boundaryCorr_.boundaryField()[patchI] = 0.0;
				boundaryCorr_.boundaryField()[patchI] == 0.0;
        }
    }
	volScalarField corrS2=fvc::average(fvc::interpolate(boundaryCorr_,"limitedScheme"));
    forAll(boundary, patchI)
    {
        if (!boundary[patchI].coupled())
        {
				
				boundaryCorr_.boundaryField()[patchI] =

					( (corrS1.boundaryField()[patchI].patchInternalField())-
					  (corrS2.boundaryField()[patchI].patchInternalField()) )
					/(boundaryCorr_.boundaryField()[patchI].patchInternalField());
				
				
				boundaryCorr_.boundaryField()[patchI] ==
					( (corrS1.boundaryField()[patchI].patchInternalField())-
					  (corrS2.boundaryField()[patchI].patchInternalField()) )
					/(boundaryCorr_.boundaryField()[patchI].patchInternalField());
        }
    }
    boundaryCorr_.internalField()=0;
    
    
    

	
}

*/



{ // not being used
		//const fvMesh& mesh = alpha1_.mesh();
		//const volVectorField& Cc = mesh.C();
		//const surfaceVectorField& Cf = mesh.Cf();
		//const labelList& owner = mesh.faceOwner();
		//const labelList& neighbour = mesh.faceNeighbour();
		//const fvBoundaryMesh& patches = mesh.boundary();
		//forAll(neighbour, facei)
			//deltaCC_[facei] = Cc[neighbour[facei]] - Cc[owner[facei]];
		//forAll(patches, patchI)
		 //if(patches[patchI].size())
		 //{
		    //const vectorField&  pCf = Cf.boundaryField()[patchI];;
		    //vectorField&  pdeltaCC_ = deltaCC_.boundaryField()[patchI];
		    //label pStart = patches[patchI].patch().start();
	         //forAll(pCf, pfI)
					//pdeltaCC_[pfI] += pCf[pfI] - Cc[owner[pStart+pfI]];

		  //if (deltaCC_.boundaryField()[patchI].coupled())
          //{
			 //const vectorField  pCcN = Cc.boundaryField()[patchI].patchNeighbourField();
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

	const dictionary pimple =alpha1_.mesh().solutionDict().subDict("PIMPLE");
    //cAlpha_=   readScalar( pimple.lookup("cAlpha") );
    pcThicknessFactor_=   readScalar( pimple.lookup("pcThicknessFactor") );
    word surfaceForceModel(pimple.lookup("surfaceForceModel"));
    word alphafModel(pimple.lookup("alphafModel"));




	scalar stime=alpha1_.time( ).elapsedCpuTime();


    const fvMesh& mesh = alpha1_.mesh();
    const surfaceVectorField& Sf = mesh.Sf();
    const fvBoundaryMesh& boundary = mesh.boundary(); 



	alphaSh_ == curtailBADOFSET((1./pcThicknessFactor_)*alpha1_-(0.5/pcThicknessFactor_-0.5),0.0); alphaSh_.correctBoundaryConditions();
	if (relaxDelS<0.0001) deltaS_= fvc::snGrad(alphaSh_);
	else                  deltaS_= 0.5*(deltaS_+fvc::snGrad(alphaSh_));


	if(alphafModel=="advect")
		
		alpha1f_ = linearInterpolate(alpha1_);
	else
		alpha1f_ = linearInterpolate(alpha1_);


	volScalarField alpha1Stmp=curtailBADOFSET(alpha1_*1.02-0.01,1.0e-3);

	volScalarField a1a2 = 0.95+0.1*sqrt(alpha1Stmp*(1.0-alpha1Stmp));
	alpha1Stmp = a1a2*alpha1Stmp+(1.0-a1a2)*fvc::average(alpha1f_*Internalfaces1_)/fvc::average(Internalfaces1_);
	alpha1Stmp.correctBoundaryConditions();

	a1a2 = 2.0*sqrt(alpha1Stmp*(1.0-alpha1Stmp));

	alpha1Stmp = a1a2*alpha1Stmp+(1.0-a1a2)*fvc::average(fvc::interpolate(alpha1Stmp)*Internalfaces1_)/fvc::average(Internalfaces1_); 
	alpha1Stmp.correctBoundaryConditions();

	alpha1Stmp = 0.99*alpha1Stmp+0.01*fvc::average(fvc::interpolate(alpha1Stmp)*Internalfaces1_)/fvc::average(Internalfaces1_);


	alpha1S_ = atanh(1.8*alpha1Stmp-0.9);
	alpha1S_.correctBoundaryConditions();





	surfaceScalarField  stf
	(	IOobject( "stf", alpha1_.time().timeName(), alpha1_.mesh() ),
		alpha1_.mesh(),  dimensionedScalar("stf", dimPressure/dimLength*dimArea, 0.0)
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

