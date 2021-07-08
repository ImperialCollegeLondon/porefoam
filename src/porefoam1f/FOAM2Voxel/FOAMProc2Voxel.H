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
* 
Please see our website for relavant literature:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
\*-------------------------------------------------------------------------*/

#include <fstream>
#include <iostream>
#include <vector>
#include <assert.h>

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"

#include "OFstream.H"

#include "voxelImage.h"
 
using namespace Foam;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:





int procMain(int argc, char *argv[], 
  voxelImage& vximage, voxelField<float>& pVoxel, voxelField<float>& vface0, voxelField<float>& vface1, voxelField<float>& vface2, voxelField<float>& peVoxel, bool initialise)
{
	Foam::argList args(argc, argv, true, true);//, initialise
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"




	int3 nn = vximage.size3();
	dbl3 xmin=vximage.X0();
	dbl3 dx=vximage.dx();
	cout<<"\n  nn: "<<nn<<"   dx: "<<dx<<"   X0: "<<xmin<<std::endl;

	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);


	{
		volVectorField U
		(	IOobject
			(	"U", runTime.timeName(), mesh,
				IOobject::MUST_READ, IOobject::NO_WRITE
			), mesh
		);

		Info<<"\nmaxU "<<max(U)<<endl;

		const vectorField & V=U.internalField();

		volScalarField pp
		(	IOobject
			(	"p", runTime.timeName(), mesh, 
				IOobject::MUST_READ, IOobject::NO_WRITE
			), mesh
		);	
		#include "createPhi.H"

		{
			Info<<"faceVelocitiesWithBC procMain"<<endl;

			const surfaceVectorField& Cf = mesh.Cf();
			const surfaceVectorField& Sf = mesh.Sf();

			forAll(Cf,c)
			{

				int i=(Cf[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
				int j=(Cf[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
				int k=(Cf[c][2]-xmin[2]+dx[2]*0.4)/dx[2];

				vector nf=(Sf[c]/mag(Sf[c]));
				if      (mag(nf[0])>0.9)
					vface0(i,j,k)=phi[c]*nf[0]/mag(Sf[c]);     
				else if (mag(nf[1])>0.9)
					vface1(i,j,k)=phi[c]*nf[1]/mag(Sf[c]);     
				else if (mag(nf[2])>0.9)
					vface2(i,j,k)=phi[c]*nf[2]/mag(Sf[c]);     
				else 
					Info<<c<<"  face not in x nor y nor z direction, nf = "<<nf<<endl;

			}
		 
			const fvBoundaryMesh& patches = mesh.boundary(); 

			forAll(patches, patchI)
			{
					Info<<patches[patchI].name()<<Cf.boundaryField()[patchI].size()<<endl;

					const Field<scalar> & phip = phi.boundaryField()[patchI];
					const Field<vector> & Cfp = Cf.boundaryField()[patchI];
					const Field<vector> & Sfp = Sf.boundaryField()[patchI];


					forAll(Cfp, c)
					{
						int i=(Cfp[c][0]-xmin[0])/dx[0] +0.4;
						int j=(Cfp[c][1]-xmin[1])/dx[1] +0.4;
						int k=(Cfp[c][2]-xmin[2])/dx[2] +0.4;

					  //if( i<nn[0] && j<nn[1] && k<nn[2] )
					  {
							vector nf=(Sfp[c]/mag(Sfp[c]));
						if      (mag(nf[0])>0.9  &&  i<nn[0]+1 && j<nn[1] && k<nn[2] && i>=0 && j>=0 && k>=0)
							vface0(i,j,k)=phip[c]*nf[0]/mag(Sfp[c]);     
						else if (mag(nf[1])>0.9  &&  i<nn[0] && j<nn[1]+1 && k<nn[2] && i>=0 && j>=0 && k>=0)
							vface1(i,j,k)=phip[c]*nf[1]/mag(Sfp[c]);     
						else if (mag(nf[2])>0.9  &&  i<nn[0] && j<nn[1] && k<nn[2]+1 && i>=0 && j>=0 && k>=0)
							vface2(i,j,k)=phip[c]*nf[2]/mag(Sfp[c]);     
						else 
							Info<<c<<"  face not in x nor y nor z direction, or not on boundaries, nf = "<<nf<<"   i"<<i<<" j"<<j<<" k"<<k<<" v"<<vface0(i,j,k)<<" cfp"<<Cfp[c]<<endl;
					  }

					}
				}
		}

		const scalarField & p=pp.internalField();

		const vectorField & C =	mesh.C();
		Info<<"C.size(): "<<C.size()<<endl;
		Info<<"V.size(): "<<V.size()<<endl;
		forAll(C,c)
		{
			int i=(C[c][0]-xmin[0])/dx[0];
			int j=(C[c][1]-xmin[1])/dx[1];
			int k=(C[c][2]-xmin[2])/dx[2];

			pVoxel(i,j,k)=p[c];
			vximage(i,j,k)=0;
		}
	}


	if(!peVoxel.nx()) { Info<<"skipping psi, electric/diffusion potential"<<endl; }
	else
	{
		volScalarField psi
		(	IOobject
			(	"psi", runTime.timeName(), mesh, 
				IOobject::READ_IF_PRESENT, IOobject::NO_WRITE
			), mesh, dimensionedScalar("psi", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0),
			  "zeroGradient"
		);

		const scalarField & p=psi.internalField();

		const vectorField & C =	mesh.C();
		Info<<"C.size(): "<<C.size()<<"  p.size(): "<<p.size()<<"  peVoxel.size(): "<<peVoxel.nx()<<endl;
		forAll(C,c)
		{
			int i=(C[c][0]-xmin[0])/dx[0];
			int j=(C[c][1]-xmin[1])/dx[1];
			int k=(C[c][2]-xmin[2])/dx[2];
			peVoxel(i,j,k)=p[c];
		}

	}



    Info<< "end" << endl;

    return 0;
}











int procMainUx(int argc, char *argv[], voxelField<float>& vface0, int3& nn,  dbl3& xmin,  dbl3& dx, bool initialise)
{
	Foam::argList args(argc, argv, true, true);//, initialise
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"



	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);


	
	volVectorField U
	(	IOobject
		(	"U", runTime.timeName(), mesh,
			IOobject::MUST_READ, IOobject::NO_WRITE
		), mesh
	);
    
	Info<<"\nmaxU "<<max(U)<<endl;

	const vectorField & V=U.internalField();
	

	surfaceScalarField phi
	(	IOobject
		(	"phi",  runTime.timeName(),  mesh,
			IOobject::READ_IF_PRESENT,  IOobject::NO_WRITE
		), linearInterpolate(U) & mesh.Sf()
	);
	//if (outputFormat=="binary")
	{
		Info<<"faceVelocitiesWithBC"<<endl;
		
		const surfaceVectorField& Cf = mesh.Cf();
		const surfaceVectorField& Sf = mesh.Sf();


		forAll(Cf,c)
		{

			int i=(Cf[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
			int j=(Cf[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
			int k=(Cf[c][2]-xmin[2]+dx[2]*0.4)/dx[2];
	
			vector nf=(Sf[c]/mag(Sf[c]));
			if      (mag(nf[0])>0.9)
				vface0(i,j,k)=phi[c]*nf[0]/mag(Sf[c]);     
			else if (mag(nf[1])>0.9)
				;//vface1(i,j,k)=phi[c]*nf[1]/mag(Sf[c]);     
			else if (mag(nf[2])>0.9)
				;//vface2(i,j,k)=phi[c]*nf[2]/mag(Sf[c]);     
			else 
				Info<<c<<"  face not in x nor y nor z direction, nf = "<<nf<<endl;

		}

		const fvBoundaryMesh& patches = mesh.boundary(); 

		forAll(patches, patchI)
		{
				Info<<patches[patchI].name()<<Cf.boundaryField()[patchI].size()<<endl;

				const Field<scalar> & phip = phi.boundaryField()[patchI];
				const Field<vector> & Cfp = Cf.boundaryField()[patchI];
				const Field<vector> & Sfp = Sf.boundaryField()[patchI];


				forAll(Cfp, c)
				{

					int i=(Cfp[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
					int j=(Cfp[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
					int k=(Cfp[c][2]-xmin[2]+dx[2]*0.4)/dx[2];
				
				  {
						vector nf=(Sfp[c]/mag(Sfp[c]));
					if      (mag(nf[0])>0.9  &&  i<nn[0]+1 && j<nn[1] && k<nn[2] && i>=0 && j>=0 && k>=0)
						vface0(i,j,k)=phip[c]*nf[0]/mag(Sfp[c]);     
					else if (mag(nf[1])>0.9  &&  i<nn[0] && j<nn[1]+1 && k<nn[2] && i>=0 && j>=0 && k>=0)
						;//vface1(i,j,k)=phip[c]*nf[1]/mag(Sfp[c]);     
					else if (mag(nf[2])>0.9  &&  i<nn[0] && j<nn[1] && k<nn[2]+1 && i>=0 && j>=0 && k>=0)
						;//vface2(i,j,k)=phip[c]*nf[2]/mag(Sfp[c]);     
					else 
						Info<<c<<"  face not in x (nor y nor z) direction, or not on boundaries, nf = "<<nf<<"   i"<<i<<" j"<<j<<" k"<<k<<" v"<<vface0(i,j,k)<<endl;
				  }

				}
			}
	}

	//const scalarField & p=pp.internalField();


	const vectorField & C =	mesh.C();
	Info<<"C.size(): "<<C.size()<<endl;
	Info<<"V.size(): "<<V.size()<<endl;



    Info<< "end" << endl;

    return 0;
}



int procMainUy(int argc, char *argv[], voxelField<float>& vface1, int3& nn,  dbl3& xmin,  dbl3& dx, bool initialise)
{
	Foam::argList args(argc, argv, true, true);//, initialise
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"



	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);


	
	volVectorField U
	(	IOobject
		(	"U", runTime.timeName(), mesh,
			IOobject::MUST_READ, IOobject::NO_WRITE
		), mesh
	);
    
	Info<<"\nmaxU "<<max(U)<<endl;

	const vectorField & V=U.internalField();
	
	//volScalarField pp
	//(	IOobject
		//(	"p", runTime.timeName(), mesh, 
			//IOobject::MUST_READ, IOobject::NO_WRITE
		//), mesh
	//);

	surfaceScalarField phi
	(	IOobject
		(	"phi",  runTime.timeName(),  mesh,
			IOobject::READ_IF_PRESENT,  IOobject::NO_WRITE
		), linearInterpolate(U) & mesh.Sf()
	);

	//if (outputFormat=="binary")
	{
		Info<<"faceVelocitiesWithBC"<<endl;
		
		const surfaceVectorField& Cf = mesh.Cf();
		const surfaceVectorField& Sf = mesh.Sf();


		forAll(Cf,c)
		{

			int i=(Cf[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
			int j=(Cf[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
			int k=(Cf[c][2]-xmin[2]+dx[2]*0.4)/dx[2];

			vector nf=(Sf[c]/mag(Sf[c]));
			if      (mag(nf[0])>0.9)
				;//vface0(i,j,k)=phi[c]*nf[0]/mag(Sf[c]);     
			else if (mag(nf[1])>0.9)
				vface1(i,j,k)=phi[c]*nf[1]/mag(Sf[c]);     
			else if (mag(nf[2])>0.9)
				;//vface2(i,j,k)=phi[c]*nf[2]/mag(Sf[c]);     
			else 
				Info<<c<<"  face not in x nor y nor z direction, nf = "<<nf<<endl;

		}
	 
		const fvBoundaryMesh& patches = mesh.boundary(); 

		forAll(patches, patchI)
		{
				Info<<patches[patchI].name()<<Cf.boundaryField()[patchI].size()<<endl;

				const Field<scalar> & phip = phi.boundaryField()[patchI];
				const Field<vector> & Cfp = Cf.boundaryField()[patchI];
				const Field<vector> & Sfp = Sf.boundaryField()[patchI];


				forAll(Cfp, c)
				{

					int i=(Cfp[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
					int j=(Cfp[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
					int k=(Cfp[c][2]-xmin[2]+dx[2]*0.4)/dx[2];
				
				  {
						vector nf=(Sfp[c]/mag(Sfp[c]));
					if      (mag(nf[0])>0.9  &&  i<nn[0]+1 && j<nn[1] && k<nn[2] && i>=0 && j>=0 && k>=0)
						;//vface0(i,j,k)=phip[c]*nf[0]/mag(Sfp[c]);     
					else if (mag(nf[1])>0.9  &&  i<nn[0] && j<nn[1]+1 && k<nn[2] && i>=0 && j>=0 && k>=0)
						vface1(i,j,k)=phip[c]*nf[1]/mag(Sfp[c]);     
					else if (mag(nf[2])>0.9  &&  i<nn[0] && j<nn[1] && k<nn[2]+1 && i>=0 && j>=0 && k>=0)
						;//vface2(i,j,k)=phip[c]*nf[2]/mag(Sfp[c]);     
					else 
						Info<<c<<"  face not in y (nor x nor z) direction, or not on boundaries, nf = "<<nf<<"   i"<<i<<" j"<<j<<" k"<<k<<" v"<<vface1(i,j,k)<<endl;
				  }

				}
			}
	}

	//const scalarField & p=pp.internalField();


	const vectorField & C =	mesh.C();
	Info<<"C.size(): "<<C.size()<<endl;
	Info<<"V.size(): "<<V.size()<<endl;


    Info<< "end" << endl;

    return 0;
}



int procMainUz(int argc, char *argv[], voxelField<float>& vface2, int3& nn,  dbl3& xmin,  dbl3& dx, bool initialise)
{
	Foam::argList args(argc, argv, true, true);//, initialise
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"



	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);


	
	volVectorField U
	(	IOobject
		(	"U", runTime.timeName(), mesh,
			IOobject::MUST_READ, IOobject::NO_WRITE
		), mesh
	);
    
	Info<<"\nmaxU "<<max(U)<<endl;

	const vectorField & V=U.internalField();

	
	surfaceScalarField phi
	(	IOobject
		(	"phi",  runTime.timeName(),  mesh,
			IOobject::READ_IF_PRESENT,  IOobject::NO_WRITE
		), linearInterpolate(U) & mesh.Sf()
	);
	//if (outputFormat=="binary")
	{
		Info<<"faceVelocitiesWithBC"<<endl;
		
		const surfaceVectorField& Cf = mesh.Cf();
		const surfaceVectorField& Sf = mesh.Sf();


		forAll(Cf,c)
		{

			int i=(Cf[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
			int j=(Cf[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
			int k=(Cf[c][2]-xmin[2]+dx[2]*0.4)/dx[2];

			vector nf=(Sf[c]/mag(Sf[c]));
			if      (mag(nf[0])>0.9)
				;//vface0(i,j,k)=phi[c]*nf[0]/mag(Sf[c]);     
			else if (mag(nf[1])>0.9)
				;//vface1(i,j,k)=phi[c]*nf[1]/mag(Sf[c]);     
			else if (mag(nf[2])>0.9)
				vface2(i,j,k)=phi[c]*nf[2]/mag(Sf[c]);     
			else 
				Info<<c<<"  face not in x nor y nor z direction, nf = "<<nf<<endl;

		}
	 
		const fvBoundaryMesh& patches = mesh.boundary(); 

		forAll(patches, patchI)
		{
				Info<<patches[patchI].name()<<Cf.boundaryField()[patchI].size()<<endl;

				const Field<scalar> & phip = phi.boundaryField()[patchI];
				const Field<vector> & Cfp = Cf.boundaryField()[patchI];
				const Field<vector> & Sfp = Sf.boundaryField()[patchI];


				forAll(Cfp, c)
				{

					int i=(Cfp[c][0]-xmin[0]+dx[0]*0.4)/dx[0];
					int j=(Cfp[c][1]-xmin[1]+dx[1]*0.4)/dx[1];
					int k=(Cfp[c][2]-xmin[2]+dx[2]*0.4)/dx[2];
				
				  {
						vector nf=(Sfp[c]/mag(Sfp[c]));
					if      (mag(nf[0])>0.9  &&  i<nn[0]+1 && j<nn[1] && k<nn[2] && i>=0 && j>=0 && k>=0)
						;//vface0(i,j,k)=phip[c]*nf[0]/mag(Sfp[c]);     
					else if (mag(nf[1])>0.9  &&  i<nn[0] && j<nn[1]+1 && k<nn[2] && i>=0 && j>=0 && k>=0)
						;//vface1(i,j,k)=phip[c]*nf[1]/mag(Sfp[c]);     
					else if (mag(nf[2])>0.9  &&  i<nn[0] && j<nn[1] && k<nn[2]+1 && i>=0 && j>=0 && k>=0)
						vface2(i,j,k)=phip[c]*nf[2]/mag(Sfp[c]);     
					else 
						Info<<c<<"  face not in z (nor x nor y) direction, or not on boundaries, nf = "<<nf<<"   i"<<i<<" j"<<j<<" k"<<k<<" v"<<vface2(i,j,k)<<endl;
				  }

				}
			}
	}


	const vectorField & C =	mesh.C();
	Info<<"C.size(): "<<C.size()<<endl;
	Info<<"V.size(): "<<V.size()<<endl;




    Info<< "end" << endl;

    return 0;
}




int procMainP(int argc, char *argv[], voxelField<float>& pVoxel, int3& nn,  dbl3& xmin,  dbl3& dx,  bool initialise)
{
	Foam::argList args(argc, argv, true, true);//, initialise
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"





	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);

	volScalarField pp
	(	IOobject
		(	"p", runTime.timeName(), mesh, 
			IOobject::MUST_READ, IOobject::NO_WRITE
		), mesh
	);

	const scalarField & p=pp.internalField();

	const vectorField & C =	mesh.C();
	Info<<"C.size(): "<<C.size()<<endl;
	forAll(C,c)
	{
		int i=(C[c][0]-xmin[0])/dx[0];
		int j=(C[c][1]-xmin[1])/dx[1];
		int k=(C[c][2]-xmin[2])/dx[2];
		pVoxel(i,j,k)=p[c];
	}

   Info<< "end" << endl;

   return 0;
}


///. should be run first,
int procMainV(int argc, char *argv[], voxelImage& vximage, int3& nn,  dbl3& xmin,  dbl3& dx, bool initialise)
{
	Foam::argList args(argc, argv, true, true);//, initialise
	if (!args.checkRootCase())
	{
		Foam::FatalError.exit();
	}

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"



	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);



	const vectorField & C =	mesh.C();
	Info<<"C.size(): "<<C.size()<<endl;
	forAll(C,c)
	{

		int i=(C[c][0]-xmin[0])/dx[0];
		int j=(C[c][1]-xmin[1])/dx[1];
		int k=(C[c][2]-xmin[2])/dx[2];
		vximage(i,j,k)=0;
	}

    Info<< "end" << endl;

    return 0;
}

int procMainPE(int argc, char *argv[], voxelField<float>& peVoxel, int3& nn,  dbl3& xmin,  dbl3& dx,  bool initialise)
{
	if(!peVoxel.nx()) { Info<<"skipping psi, electric/diffusion potential"<<endl;  return 0; }
	Foam::argList args(argc, argv, true, true);//, initialise
    if (!args.checkRootCase())
    {
        Foam::FatalError.exit();
    }

#   include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);
#   include "createNamedMesh.H"





	runTime.setTime(timeDirs[timeDirs.size()-1], timeDirs.size()-1);

	volScalarField psi
	(	IOobject
		(	"psi", runTime.timeName(), mesh, 
			IOobject::READ_IF_PRESENT, IOobject::NO_WRITE
		), mesh, dimensionedScalar("psi", dimensionSet(0, 2, -2, 0, 0, 0, 0), 0),
        "zeroGradient"
	);

	const scalarField & p=psi.internalField();

	const vectorField & C =	mesh.C();
	Info<<"C.size(): "<<C.size()<<endl;
	forAll(C,c)
	{
		int i=(C[c][0]-xmin[0])/dx[0];
		int j=(C[c][1]-xmin[1])/dx[1];
		int k=(C[c][2]-xmin[2])/dx[2];
		peVoxel(i,j,k)=p[c];
	}

   Info<< "end" << endl;

   return 0;
}






// ************************************************************************* //
