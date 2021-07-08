/*-------------------------------------------------------------------------*\
 Upscaling flow field / Abs. permeability computation
 Outdated, use calc_distributions instead

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

#include <fstream>
#include <assert.h>

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"


using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
#   include "setRootCase.H"
#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
    

	#   include "createNamedMesh.H"


    


	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);

			
	volVectorField U
	(
		IOobject
		(
			"U",   runTime.timeName(),
			mesh,  IOobject::MUST_READ
		),
		mesh
	);

    #include "createPhi.H"
    Info<< "Reading transportProperties\n" << endl;
    IOdictionary twoPhaseProperties
    (
        IOobject
        (
            "transportProperties",  runTime.constant(),
            runTime,    IOobject::MUST_READ
        )
    );
    
	const dimensionedScalar rho(twoPhaseProperties.lookup("rho"));
	const dimensionedScalar nu(twoPhaseProperties.lookup("nu"));
	const dimensionedScalar muEff(rho*nu);

	volScalarField p
	(
		IOobject
		(
			"p",   runTime.timeName(),
			mesh,  IOobject::MUST_READ
		),
		mesh
	);


	scalar V=sum(mesh.V()).value();
	
	#define x_ 0 
	#define y_ 1 
	#define z_ 2 

	scalar L[3];
	scalar K[3];
	vector dp;
	vector VDarcy;
	scalar dx=pow(average(mesh.V()), 1./3.).value();
	L[x_]=(gMax(mesh.points().component(0))-gMin(mesh.points().component(0)));
	L[y_]=(gMax(mesh.points().component(1))-gMin(mesh.points().component(1)));
	L[z_]=(gMax(mesh.points().component(2))-gMin(mesh.points().component(2)));
	Info<<" "<< dx <<" "<< L[x_] <<" "<< L[y_] <<" "<< L[z_] <<" "<< endl;

	scalar A[]={L[y_]*L[z_],L[z_]*L[x_],L[x_]*L[y_]};
	word LeftPs[]={"Left","Bottom","Front"};//
	word RightPs[]={"Right","Top","Back"};//
	word directions[]={"x","y","z"};//
		
	for (int i=0;i<3;i++)
	{

		label iLeft = mesh.boundaryMesh().findPatchID(LeftPs[i]);
		if (iLeft < 0)
		{
			Info	<< "Unable to find  patch " << LeftPs[i] << nl	<< endl;
			continue;
		}

		label iRight = mesh.boundaryMesh().findPatchID(RightPs[i]);
		if (iRight < 0)
		{
			Info	<< "Unable to find  patch " << RightPs[i] << nl	<< endl;
			Info	<< "Patchs "<<LeftPs[i]<<" and "<<RightPs[i]<<" should be normal to "<<i<<" direction (0:x,1:y,2:z)";
			continue;
		}

	 

		scalar fluxIn=gSum(phi.boundaryField()[iLeft]);



		scalar PLeft=gSum(p.boundaryField()[iLeft]*(phi.boundaryField()[iLeft]))/(fluxIn+1e-29);
		scalar PRight=gSum(p.boundaryField()[iRight]*(phi.boundaryField()[iRight]))/(fluxIn+1e-29);
			dp[i]=mag(PLeft-PRight);

		K[i]=mag(fluxIn)*muEff.value()/A[i]*L[i]/(dp[i]+1e-26);
		VDarcy[i]=mag(fluxIn/A[i]);
	}
	label i=findMax(dp);
	scalar Pmax_min=(max(p)-min(p)).value();
	scalar Umax=(max(mag(U))).value();
	scalar Re=rho.value()*VDarcy[i]*std::sqrt(K[i])/muEff.value() ;
	



    Info << runTime.caseName()
      <<"\tK_"<<directions[i]<<"= "<<  K[i]<<" m^2 \t"
      <<"\teffPorosity= "<< V/(L[x_]*L[y_]*L[z_]) <<" \t"
      <<"\tDarcyVelocity= "<< VDarcy[i] <<" m/s \t"
      <<"\tDelP= "<< dp[i] <<" Pa \t"
      <<"\tK: "<<  K[i]<<" m^2 \t"
      <<"\t\tK=( "<<  K[x_]<<"  "<<  K[y_]<<"  "<<  K[z_] <<" )"
      <<" Pmax - Pmin "<<  Pmax_min
      <<" Umax "<<  Umax
      <<"\t\t effPorosity = "<<  "V_pore / (L_x*L_y*L_z) = "<<  V<<" / "<<" ("<<L[x_]<<"*"<<L[y_]<<"*"<<L[z_]<<") = "<<V/(L[x_]*L[y_]*L[z_])
      <<"\t\t Re = "<<  "rho*VDarcy*sqrt(K)/mu = "<< rho.value() <<" * "<<VDarcy[i]<<" * sqrt("<<K[i]<<") / "<<muEff.value()<<" ) = "<<Re 
      << "\n";



	if (Pstream::master())
	{
		std::ofstream of("summary_calcPerm.txt",std::ios::app);
		assert(of);

		of<< runTime.caseName() 
		  <<"\tK_"<<directions[i]<<"= "<<  K[i]<<" m^2 \t"
		  <<"\teffPorosity= "<< V/(L[x_]*L[y_]*L[z_]) <<" \t"
		  <<"\tDarcyVelocity= "<< VDarcy[i] <<" m/s \t"
		  <<"\tDelP= "<< dp[i] <<" Pa \t"
		  <<"\tK: "<<  K[i]<<" m^2 \t"
		  <<"\t\tK=( "<<  K[x_]<<"  "<<  K[y_]<<"  "<<  K[z_] <<" )"
		  <<" Pmax - Pmin "<<  Pmax_min
		  <<" Umax "<<  Umax
		  <<"\t\t effPorosity = "<<  "V_pore / (L_x*L_y*L_z) = "<<  V<<" / "<<" ("<<L[x_]<<"*"<<L[y_]<<"*"<<L[z_]<<") = "<<V/(L[x_]*L[y_]*L[z_])
		  <<"\t\t Re = "<<  "rho*VDarcy*sqrt(K)/mu = "<< rho.value() <<" * "<<VDarcy[i]<<" * sqrt("<<K[i]<<") / "<<muEff.value()<<" ) = "<<Re 
		  << "\n";
		  
		
		of.close();
	}
	
	
    Info<< "end" << endl;

    return 0;
}


// ************************************************************************* //
