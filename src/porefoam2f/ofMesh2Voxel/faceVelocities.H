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
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
\*-------------------------------------------------------------------------*/

    
	voxelField<double> vface0(n,0.0);
	voxelField<double> vface1(n,0.0);
	voxelField<double> vface2(n,0.0);
    const faceList& faces = mesh.faces();
    const surfaceVectorField& Cf = mesh.Cf();
    const surfaceVectorField& Sf = mesh.Sf();
    //const surfaceVectorField nf = mesh.Sf()/mag(mesh.Sf());

    #include "createPhi.H"
	//Info<<"phi.size(): "<<phi.size()<<endl;
	//Info<<"faces.size(): "<<faces.size()<<endl;
	//Info<<"Cf.size(): "<<Cf.size()<<endl;
	//Info<<"Cf.size(): "<<Cf.boundaryField().size()<<endl;

	forAll(Cf,c)
	{

		int i=(Cf[c][0]+dx[0]/4.0)/dx[0];
		int j=(Cf[c][1]+dx[1]/4.0)/dx[1];
		int k=(Cf[c][2]+dx[2]/4.0)/dx[2];
                
        vector nf=(Sf[c]/mag(Sf[c]));
        if      (mag(nf[0])>0.95)
            vface0[k][j][i]=phi[c]*nf[0]/mag(Sf[c]);     
        else if (mag(nf[1])>0.95)
            vface1[k][j][i]=phi[c]*nf[1]/mag(Sf[c]);     
        else if (mag(nf[2])>0.95)
            vface2[k][j][i]=phi[c]*nf[2]/mag(Sf[c]);     
        else 
            Info<<c<<"  face not in x nor y nor z direction, nf = "<<nf<<endl;

	}
 
    const fvBoundaryMesh& patches = mesh.boundary(); 

    forAll(patches, patchI)
    {

            const Field<scalar> & phip = phi.boundaryField()[patchI];
            const Field<vector> & Cfp = Cf.boundaryField()[patchI];
            //const Field<vector> & nfp = nf.boundaryField()[patchI];
            const Field<vector> & Sfp = Sf.boundaryField()[patchI];

            const faceList& faces = mesh.faces();

            forAll(phip, c)
            {
                //label faceI=patches[patchI].patch().start()+i;
                //const point pi = Cf[faceI];
                int i=(Cfp[c][0]+dx[0]/4.0)/dx[0];
                int j=(Cfp[c][1]+dx[1]/4.0)/dx[1];
                int k=(Cfp[c][2]+dx[2]/4.0)/dx[2];
            
              if( i<n[0] && j<n[1] && k<n[2] )
              {
                    vector nf=(Sfp[c]/mag(Sfp[c]));
                if      (mag(nf[0])>0.95)
                    vface0[k][j][i]=phip[c]*nf[0]/mag(Sfp[c]);     
                else if (mag(nf[1])>0.95)
                    vface1[k][j][i]=phip[c]*nf[1]/mag(Sfp[c]);     
                else if (mag(nf[2])>0.95)
                    vface2[k][j][i]=phip[c]*nf[2]/mag(Sfp[c]);     
                else 
                    Info<<c<<"  face not in x nor y nor z direction, nf = "<<nf<<endl;
              }

            }
        }

	 vface0.write(runTime.constant()+"/../Ux.dat");
	 vface1.write(runTime.constant()+"/../Uy.dat");	
	 vface2.write(runTime.constant()+"/../Uz.dat");	

