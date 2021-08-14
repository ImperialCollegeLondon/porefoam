    
 template<class Type>
 tmp<GeometricField<Type, fvsPatchField, surfaceMesh> >
 fvcDelta
 (
     const GeometricField<Type, fvPatchField, volMesh>& vf
 )
 {
     const fvMesh& mesh = vf.mesh();
 
     // construct GeometricField<Type, fvsPatchField, surfaceMesh>
     tmp<GeometricField<Type, fvsPatchField, surfaceMesh> > tssf(
         new GeometricField<Type, fvsPatchField, surfaceMesh>(
             IOobject("delta("+vf.name()+')', vf.instance(), vf.mesh()),
             mesh,  vf.dimensions()
         )
     );
     GeometricField<Type, fvsPatchField, surfaceMesh>& ssf = tssf.ref();
 

     // owner/neighbour addressing
     const labelList& owner = mesh.owner();
     const labelList& neighbour = mesh.neighbour();
 
 
     forAll(owner, faceI)
         ssf[faceI] = (vf[neighbour[faceI]] - vf[owner[faceI]]);
 
     forAll(vf.boundaryField(), bi)
         ssf.boundaryFieldRef()[bi] = 
              vf.boundaryField()[bi] - vf.boundaryField()[bi].patchInternalField();
 
     return tssf;
 }     
 
 
 //#include "mylib.h"


scalar volIntegrate(const volScalarField& vf)  {
	return (sum(vf*vf.mesh().V())).value();
}

//scalar volAverage(const volScalarField & vf,const volScalarField & clipV)
//{
	//return (sum(vf*clipV)/sum(clipV)).value();
//}

scalar surfIntegrate(const surfaceScalarField& sf)  {
	return (sum(sf*sf.mesh().magSf())).value();
}

scalar surfAverage(const surfaceScalarField& sf,const surfaceScalarField& deltaMagSf)  {
	return (sum(sf*deltaMagSf)/sum(deltaMagSf)).value();
}


