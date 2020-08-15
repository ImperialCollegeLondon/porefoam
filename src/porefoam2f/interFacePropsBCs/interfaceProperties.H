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
Ali Q Raeini:    a.q.raeini@imperial.ac.uk
Mosayeb Shams:    m.shams14@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk

Description:
	interface force and capillary pressure computation

\*-------------------------------------------------------------------------*/


#ifndef interfaceProperties_H
#define interfaceProperties_H

#include "IOdictionary.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "processorFvsPatchFields.H"
#include "pointVolInterpolation.H"
#include "volPointInterpolation.H"
#include "pointMesh.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class interfaceProperties Declaration
\*---------------------------------------------------------------------------*/

class interfaceProperties
{
    // Private data

	//- Keep a reference to the transportProperties dictionary
	const dictionary& transportPropertiesDict_;

	//- Compression coefficient
	scalar pcThicknessFactor_;
	scalar fcCorrectTangent_;
	scalar fcCorrectTangentRelax_;
	scalar  fcdFilter_;
	scalar nPcNonOrthCorr_;
	label  pcRefCell, pcRefCellOrig_;
	scalar pcRefValue, pcRefValueOrig_;

	label  smoothingKernel_;
	scalar  smoothingRelaxFactor_;
	label  wallSmoothingKernel_;
	//- Surface tension
	dimensionedScalar sigma_;

	//- Stabilisation for normalisation of the interface normal
	const dimensionedScalar deltaN_;

	const volScalarField& alpha1_;
	volScalarField alpha1S_;
	const volVectorField& U_;

	volScalarField pc_;
	surfaceScalarField deltaS_;
	volScalarField alphaSh_; 
	//volScalarField K_; 

	surfaceScalarField nHatf_;

	volVectorField gPc_; //! gPc_ used for both visualization and post-processing
	surfaceScalarField sgPc_;
	surfaceScalarField alpha1f_;
	//volVectorField nS_;
	surfaceVectorField nw_;
	labelList			 edgemarks_; ///. 1 for boundaryInternalEdges
	labelList			 boundaryInternalFaces_;
	surfaceScalarField Internalfaces1_;
	surfaceScalarField BInternalfs_;
	volScalarField AvgInternFaces1_;
	surfaceScalarField sgPcErr_;
	surfaceScalarField sgPcErrn_; // not needed here, only for debugging and visualization, in case


	pointMesh pMesh_;
	volPointInterpolation vpi_;
	//pointVolInterpolation pvi_;
	//labelList interfPointsOld_;
	pointVectorField distPointInterface_;
	//pointVectorField nSpCurl_;
	//pointScalarField alpha1P;
	//surfaceVectorField deltaCC_;

#ifdef _USE_ALPHAF_

#endif



	// Private Member Functions

	//- Disallow default bitwise copy construct and assignment
	interfaceProperties(const interfaceProperties&);
	void operator=(const interfaceProperties&);

	//- Correction for the boundary condition on the unit normal nHat on
	//  walls to produce the correct contact dynamic angle
	//  calculated from the component of U parallel to the wall
	void correctContactAngle
	(
		surfaceVectorField::Boundary& nHat,
		volVectorField::Boundary& gradAlphaf,
		volVectorField::Boundary& nSb,
		volScalarField::Boundary& alphaSb
	) const;

        //- Re-calculate the interface curvature
	void calculateK();


	surfaceScalarField calcCurvatureFConservative
	(
		//surfaceVectorField& nSf,
		//const surfaceScalarField& delS
	);

	void calcCurvatureFSF
	(
		surfaceScalarField&  stf,
		const surfaceScalarField& delS,
		const volScalarField& a1a2
	);
	void smoothOverInterfPointsNotUsed
	(
		surfaceVectorField& curvatureForcefDiff, 
		const surfaceScalarField& magDelS, 
		const surfaceScalarField& magInterfaceSf, 
		label smoothingKernel, 
		scalar smoothingRelaxFactor
	);
	void smoothNSOverInterfPoints
	(
		pointVectorField& nHatSp, 
		const surfaceScalarField& magDelS, 
		const surfaceScalarField& magSf, 
		label smoothingKernel, 
		scalar smoothingRelaxFactor,
		const labelList& interfPoints
	);

		//- Calculate the curvature force based on the reconstructed interface
       // void calculateConservativeCurvatureForce(const surfaceVectorField&);


public:

    //- Conversion factor for degrees into radians
    static const scalar convertToRad;


    // Constructors

        //- Construct from volume fraction field gamma and IOdictionary
        interfaceProperties
        (
            const volScalarField& alpha1,
            const volVectorField& U,
            const IOdictionary&
        );


    // Member Functions


	scalar pcThicknessFactor() const          {  return pcThicknessFactor_;  }

	const dimensionedScalar& deltaN() const   {  return deltaN_;  }

	const surfaceScalarField& deltaS() const  {  return deltaS_;  }
	const surfaceScalarField& nHatf() const   {  return nHatf_;  }

	const surfaceScalarField& sgPc() const    {  return sgPc_;  }


	const volVectorField& gPc() const  {  return gPc_;  }
	const volScalarField& pc()         {  return pc_;  }
	//const volScalarField& K()          {  return K_;  }
	const surfaceScalarField& Internalfaces1()   {  return Internalfaces1_;  }
	const surfaceScalarField& BInternalfs()      {  return BInternalfs_;  }
	const volScalarField& AvgInternFaces1()      {  return AvgInternFaces1_;  }

	const dimensionedScalar& sigma() const       {  return sigma_;  }

	void correct(scalar deltaSrelax);


};







// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
