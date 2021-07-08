/*-------------------------------------------------------------------------*\

You can redistribute this code and/or modify this code under the
terms of the GNU General Public License (GPL) as published by the
Free Software Foundation, either version 3 of the License, or (at
your option) any later version. see <http://www.gnu.org/licenses/>.

Please see our website for relavant literature making use of this code:
https://www.imperial.ac.uk/earth-science/research/research-groups/pore-scale-modelling/

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/


//!\brief Converts 3D Image files into openfoam format for computing flow fields.



#include <fstream>
#include <assert.h>
#include <algorithm>
#include <vector>
#include <array>
#include <valarray>
#include <iostream>
#include <string>
#include <sstream>

using namespace std;


#include "voxelImage.h"
#include "voxelRegions.h"

#include "voxelImage.cpp"

int usage()  {
	cout<<
		"converts micro-CT images to OpenFOAM serial or simple parallel meshes\n"
		"usage:\n"
		"  voxelToFoamPar image.mhd nProcY nProcY nProcY  F/T:resetX0 F/T:keepBCs\n"
		"  examples:\n"
		"    voxelToFoamPar imageName.mhd 3 2 2 \n"
		"    voxelToFoamPar imageName.mhd 1 1 1 \n"
		"    voxelToFoamPar imageName.mhd 1 1 1 T "<<endl;
	return 1;
}


void toFoam(voxelImage& vxlImg, int nVVs, const voxelField<int>& procIsijk, int iPrc, int jPrc, int kPrc);

int main(int argc, char** argv)  {  //!- reads image, adds a boundary layers around it
	int irg = 0;
	if(argc<5)		return usage();
	string hdr(argv[++irg]); // mhd header name
	if(hdr.size()<4 || ( hdr.compare(hdr.size()-4,4,".mhd")!=0 && hdr.compare(hdr.size()-4,4,".tif")!=0) )
		return usage();
	int3 nPar(1,1,1); 	nPar.x= atoi(argv[++irg]);  nPar.y= atoi(argv[++irg]);  nPar.z= atoi(argv[++irg]);
	char resetX0 = argc>1+irg ? argv[++irg][0] : 'F';
	char keepBCs = argc>1+irg ? argv[++irg][0] : 'F';
	if (nPar.z*nPar.y*nPar.x<1) { cout<<"\nError: nPar.z*nPar.y*nPar.x<1\n"; return usage(); }

	cout<<"voxelToFoamPar  header        nPar       resetX0   keepBCs "<<endl;
	cout<<"voxelToFoamPar "<<hdr<<"  "<<nPar<<"  "<<resetX0<<"    "<<keepBCs<<endl;

	voxelImage vimage;  readConvertFromHeader(vimage, hdr);
	int3 n = vimage.size3();

	if(resetX0=='T' || resetX0=='t')
		vimage.X0Ch()=dbl3(0.,0.,0.);

	vimage.printInfo();

	vimage.writeHeader("vxlImage.mhd");

	cout <<"finding connected parts of the image "<<endl;
	if(keepBCs!='T' && keepBCs!='t')
		vimage.threshold101(0,0);
	cout <<" X0: "<<vimage.X0()<<endl;
	vimage.cropD(int3(0,0,0),n,1,1);	//		 XXXXXXXXXXXXXXXXXXXXXXXXXXXX
	cout <<" X0: "<<vimage.X0()<<endl;

	keepLargest0(vimage); //! CtrlF:isolated=254
	cout <<" X0: "<<vimage.X0()<<endl;

	vimage.printInfo();
	cout <<"Converting to OpenFOAM format "<<endl;



	//! solid phase should have the highest vv TODO: add multilabel
	int nVVs=0;   forAllvv_seq(vimage) if(vv<240) nVVs =max(nVVs,int(vv));
	cout<<" maxVxlvalue:"<<nVVs<<endl;
	forAllvp_(vimage) if(*vp>nVVs) *vp=nVVs; //! set isolateds to last vv, CtrlF:isolated=254
	++nVVs;
	cout<<"nVVs:"<<int(nVVs)<<endl;


const int
	Left  =nVVs+0,
	Right =nVVs+1,
	Bottom=nVVs+2,
	Top   =nVVs+3,
	Back  =nVVs+4,
	Front =nVVs+5;
	vimage.setSlice('i',0,    0+1*Left  );
	vimage.setSlice('i',n.x+1,0+1*Right );
	vimage.setSlice('j',0,    0+1*Bottom);
	vimage.setSlice('j',n.y+1,0+1*Top   );
	vimage.setSlice('k',0,    0+1*Back  );
	vimage.setSlice('k',n.z+1,0+1*Front );


	voxelField<int> procIsijk(nPar.x+2,nPar.y+2,nPar.z+2,-1);
	if (nPar.z*nPar.y*nPar.x==1)
		toFoam(vimage,nVVs, procIsijk, 1, 1, 1);
	else if (nPar.z*nPar.y*nPar.x>1)  {
		voxelField<voxelImage>  vimages(nPar.x,nPar.y,nPar.z,voxelImage());

		vector<int> iBs(nPar.x+1,n.x);
		vector<int> jBs(nPar.y+1,n.y);
		vector<int> kBs(nPar.z+1,n.z);
		for (int iz=0;iz<nPar.z;iz++)	kBs[iz]=int(n.z/nPar.z)*iz;
		for (int iy=0;iy<nPar.y;iy++)	jBs[iy]=int(n.y/nPar.y)*iy;
		for (int ix=0;ix<nPar.x;ix++)	iBs[ix]=int(n.x/nPar.x)*ix;


		cout<<"iBs: "<<*iBs.begin()<<" ... "<<*iBs.rbegin()<<endl;
		cout<<"jBs: "<<*jBs.begin()<<" ... "<<*jBs.rbegin()<<endl;
		cout<<"kBs: "<<*kBs.begin()<<" ... "<<*kBs.rbegin()<<endl;

			OMPragma("omp parallel for")
			for (int ixyz=0;ixyz<nPar.x*nPar.y*nPar.z;ixyz++)  {
				int ix= ixyz%nPar.x;
				int iz= ixyz/(nPar.x*nPar.y);
				int iy=(ixyz/nPar.x)%nPar.y;
				vimages(ix,iy,iz).reset(iBs[ix+1]-iBs[ix]+2, jBs[iy+1]-jBs[iy]+2, kBs[iz+1]-kBs[iz]+2,0);
				vimages(ix,iy,iz).setFrom(vimage, iBs[ix], jBs[iy], kBs[iz]);
				(cout<<"poro_"+_s(ix)+_s(iy)+_s(iz)+":"+_s(vimages(ix,iy,iz).volFraction(0,0))+"  iBx"+_s(iBs[ix+1]-iBs[ix]+2)+"\n").flush();
			};
			int iPrc=-1;
			forAllkji(vimages)
				if(vimages(i,j,k).volFraction(0,0)>1e-12)		procIsijk(i+1,j+1,k+1)=++iPrc;

			cout<<"\n************** generating meshes ************"<<endl;
			vimage.reset(0,0,0,0);
			OMPragma("omp parallel for")
			for (int ixyz=0;ixyz<nPar.x*nPar.y*nPar.z;ixyz++)  {
			 int ix= ixyz%nPar.x;
			 int iz= ixyz/(nPar.x*nPar.y);
			 int iy=(ixyz/nPar.x)%nPar.y;
			 if(procIsijk(ix+1,iy+1,iz+1)>=0)  {
				(cout<<"************* processor: "+_s(ix)+_s(iy)+_s(iz)+", Phi="+_s(100*vimages(ix,iy,iz).volFraction(0,0))+" *************\n").flush();
				toFoam(vimages(ix,iy,iz),nVVs, procIsijk, ix+1, iy+1, iz+1);
				cout<<endl;
			 }
			}
	} else cout<<"!!! npx x npy x npz = "<<nPar.z*nPar.y*nPar.x<<endl;
	cout<<":/"<<endl;
   return 0;
}

void toFoam(voxelImage& vxlImg, int nVVs, const voxelField<int>& procIsijk, int iPrc, int jPrc, int kPrc)  { //!- generate an openfaom (processor) mesh similar to voxelToFoamPar.cpp
	int myprocI=procIsijk(iPrc,jPrc,kPrc);
	int3 n=vxlImg.size3();n.x-=2;n.y-=2;n.z-=2;
	dbl3 X0=vxlImg.X0();
	dbl3 dx=vxlImg.dx();
	X0+=dx;



	string Folder = myprocI>=0 ?  "processor"+_s(myprocI)  :  ".";
	(cout<<"procIs, size: "+_s(procIsijk.size3())+",   ijk: "+_s(int3(iPrc,jPrc,kPrc))+"  >> "+Folder+"  N: "+_s(n)+"  X0: "+_s(X0)+"\n").flush();
	mkdirs(Folder);	mkdirs(Folder+"/constant");
	Folder=Folder+"/constant/polyMesh";  	cout<<"creating folder "+Folder+ _TRY_(mkdirs(Folder)) <<endl;


//=======================================================================

	voxelField<int> point_mapper(n.x+1,n.y+1,n.z+1,-1);




 {

	int iPoints=-1;
	forAllkji_1(vxlImg)
		if (!vxlImg(i,j,k))  {
				int*
				dd=&point_mapper(i-1,j-1,k-1);
				      if (*dd<0) *dd=++iPoints;
				++dd; if (*dd<0) *dd=++iPoints;
				dd+=n.x;
						if (*dd<0) *dd=++iPoints;
				++dd; if (*dd<0) *dd=++iPoints;

				dd=&point_mapper(i-1,j-1,k);
				      if (*dd<0) *dd=++iPoints;
				++dd; if (*dd<0) *dd=++iPoints;
				dd+=n.x;
						if (*dd<0) *dd=++iPoints;
				++dd; if (*dd<0) *dd=++iPoints;
		}


	cout<<"writing  points; nPoints: "+_s(iPoints+1)+";\n";cout.flush();


	ofstream pointsf(Folder+"/points");   assert(pointsf);

	pointsf<<
		"FoamFile\n"
		"{\n"
		"	version	 2.0;\n"
		"	format	  ascii;\n"
		"	class	   vectorField;\n"
		"	location	\""<<Folder+"\";\n"
		"	object	  points;\n"
		"}\n\n";

	pointsf<<iPoints+1<<"\n("<<endl;
	pointsf.precision(8);
	iPoints=-1;
	for (int iz=0;iz<point_mapper.nz();++iz)  {	double z=iz*dx[2]+X0[2];
		for (int iy=0;iy<point_mapper.ny();iy++)  {	double y=iy*dx[1]+X0[1];
			for (int ix=0;ix<point_mapper.nx();ix++)  {
				if(point_mapper(ix,iy,iz)>=0)  {
					point_mapper(ix,iy,iz)=++iPoints;//. sort point_mapper
					double x=ix*dx[0]+X0[0];
					pointsf<< "("<<x<< ' '<<y<<' '<<z<<")\n";
				}
			}
		}
	}

	pointsf<<"\n)"<<endl;

 }
	cout <<" :/"<<endl;
//=======================================================================




	size_t nCells=0;

	const int
	Internal     = 0,
	Grainwalls   = 1,
	Left  =nVVs+0,
	Right =nVVs+1,
	Bottom=nVVs+2,
	Top   =nVVs+3,
	Back  =nVVs+4,
	Front =nVVs+5;

	array<string,255> B_nams;
	B_nams[Internal]="Internal";
	B_nams[Grainwalls]="Grainwalls";
	int nBoundaries=5+nVVs;
	cout<<"ib:Nvvs: ";
	for(int ib=2;ib<nVVs;++ib)  { (cout<<ib<<" ").flush();  B_nams[ib]="VV"+_s(ib)+"B"; }
	B_nams[Left]="Left";
	B_nams[Right]="Right";
	B_nams[Bottom]="Bottom";
	B_nams[Top]="Top";
	B_nams[Back]="Back";
	B_nams[Front]="Front";

	array<size_t,255> nFaces; nFaces.fill(0);


	for (int iz=1;iz<=n.z;++iz)
	 for (int iy=1;iy<=n.y;++iy)
	  for (int ix=1;ix<=n.x;++ix)
		if (!vxlImg(ix,iy,iz))  {
			nCells++;

			unsigned char
			neiv=vxlImg(ix-1,iy,iz);
			if (ix!=1)  {
			  if (neiv)     ++nFaces[neiv]; }
			else if (neiv)  ++nFaces[neiv];
			else            ++nFaces[Left];

			neiv=vxlImg(ix+1,iy,iz);
			if (ix!=n.x)  {
			  if (neiv)     ++nFaces[neiv];
			  else          ++nFaces[Internal]; }
			else if (neiv)  ++nFaces[neiv];
			else            ++nFaces[Right];

			neiv=vxlImg(ix,iy-1,iz);
			if (iy!=1)  {
			  if (neiv)     ++nFaces[neiv]; }
			else if (neiv)  ++nFaces[neiv];
			else            ++nFaces[Bottom];

			neiv=vxlImg(ix,iy+1,iz);
			if (iy!=n.y)  {
			  if (neiv)     ++nFaces[neiv];
			  else          ++nFaces[Internal]; }
			else if (neiv)  ++nFaces[neiv];
			else            ++nFaces[Top];


			neiv=vxlImg(ix,iy,iz-1);
			if (iz!=1)  {
			  if (neiv)     ++nFaces[neiv]; }
			else if (neiv)  ++nFaces[neiv];
			else            ++nFaces[Back];

			neiv=vxlImg(ix,iy,iz+1);
			if (iz!=n.z)  {
			  if (neiv)     ++nFaces[neiv];
			  else          ++nFaces[Internal]; }
			else if (neiv)  ++nFaces[neiv];
			else            ++nFaces[Front];
		}


	(cout<<",  nCells: "+_s(nCells)+",    B:nFaces: ").flush();
	for(int ib=0;ib<255;++ib)  if(nFaces[ib])  cout<< " "+_s(ib)+":"<<nFaces[ib]<<", ";
	cout<<endl;

	nBoundaries+=	int(procIsijk(iPrc-1,jPrc,kPrc)>=0 && nFaces[Left])+
						int(procIsijk(iPrc+1,jPrc,kPrc)>=0 && nFaces[Right])+
						int(procIsijk(iPrc,jPrc-1,kPrc)>=0 && nFaces[Bottom])+
						int(procIsijk(iPrc,jPrc+1,kPrc)>=0 && nFaces[Top])+
						int(procIsijk(iPrc,jPrc,kPrc-1)>=0 && nFaces[Back])+
						int(procIsijk(iPrc,jPrc,kPrc+1)>=0 && nFaces[Front]);


	array<size_t,255> iStartFaces; iStartFaces.fill(0);





	{	ofstream boundary((Folder+"/boundary"));
		assert(boundary);
		boundary<<
		"FoamFile\n"
		"{\n"
		"	version	 2.0;\n"
		"	format	  ascii;\n"
		"	class	   polyBoundaryMesh;\n"
		"	location	\""<<Folder+"\";\n"
		"	object	  boundary;\n"
		"}\n\n";


		#define write_boundary(nwrite, name,ptype)          {          \
			iStartFaces[name] = iLastFace;										\
			boundary<<		                                                \
			"	"<<  B_nams[name]							 <<endl<<			   \
			"	{"												 <<endl<<		         \
			"		type			"<<ptype<<";"				  <<endl<<			   \
			"		nFaces		  "<<(nwrite)* nFaces[name]<<';'<<endl<<	\
			"		startFace	   "<<iLastFace<<';'   <<endl<<  			 \
			"	}"												 <<endl; 					\
			iLastFace += (nwrite)* nFaces[name];			}

		#define write_procBoundary(name,myId,neiId)     {              \
			iStartFaces[name] = iLastFace;										\
			if (nFaces[name] >0)                                       \
			boundary<<		                                                \
			"	"<<  "processor"+	_s(myId)+"to"+_s(neiId)	<<endl<<	   \
			"	{"												 <<endl<<		         \
			"		type			processor;"				  <<endl<<			      \
			"		inGroups        1(processor);"	  <<endl<<			      \
			"		nFaces		  "<<nFaces[name]<<';'	   <<endl<<	      \
			"		startFace	   "<<iLastFace<<';'   <<endl<<   \
			"		matchTolerance  0.0001;"	  <<endl<<			            \
			"		transform       unknown;"	  <<endl<<			            \
			"		myProcNo        "<<myId<<";"	  <<endl<<			         \
			"		neighbProcNo    "<<neiId<<";"	  <<endl<<			         \
			"	}"												 <<endl;						\
			iLastFace += nFaces[name];					  }


		boundary<< nBoundaries <<endl	<<'('<<endl;






		int iLastFace = nFaces[Internal];

		//write_boundary(true, Grainwalls,"patch");
		for(int ib=1;ib<nVVs;++ib)  write_boundary(1, ib,"patch");

		write_boundary(procIsijk(iPrc-1,jPrc,kPrc)<0, Left,"patch");
		write_boundary(procIsijk(iPrc+1,jPrc,kPrc)<0, Right,"patch");
		write_boundary(procIsijk(iPrc,jPrc-1,kPrc)<0, Bottom,"patch");
		write_boundary(procIsijk(iPrc,jPrc+1,kPrc)<0, Top,"patch");
	 #ifdef _2D_
		write_boundary(true, Back,"empty");
		write_boundary(true, Front,"empty");
	 #else
		write_boundary(procIsijk(iPrc,jPrc,kPrc-1)<0, Back,"patch");
		write_boundary(procIsijk(iPrc,jPrc,kPrc+1)<0, Front,"patch");
	 #endif


		if (procIsijk(iPrc-1,jPrc,kPrc)>=0) write_procBoundary(Left,myprocI, procIsijk(iPrc-1,jPrc,kPrc));
		if (procIsijk(iPrc+1,jPrc,kPrc)>=0) write_procBoundary(Right,myprocI, procIsijk(iPrc+1,jPrc,kPrc));
		if (procIsijk(iPrc,jPrc-1,kPrc)>=0) write_procBoundary(Bottom,myprocI, procIsijk(iPrc,jPrc-1,kPrc));
		if (procIsijk(iPrc,jPrc+1,kPrc)>=0) write_procBoundary(Top,myprocI, procIsijk(iPrc,jPrc+1,kPrc));
		if (procIsijk(iPrc,jPrc,kPrc-1)>=0) write_procBoundary(Back,myprocI, procIsijk(iPrc,jPrc,kPrc-1));
		if (procIsijk(iPrc,jPrc,kPrc+1)>=0) write_procBoundary(Front,myprocI, procIsijk(iPrc,jPrc,kPrc+1));

		boundary<<")"   <<endl;
		boundary.close();
	}







	cout<<"creating faces  "<<endl;


	array<vector<array<int,6> >,255> faces_bs;
	size_t sumnFaces=0;
	for(int ib=0;ib<255;++ib)  if(nFaces[ib])  {
		sumnFaces+=nFaces[ib];
		faces_bs[ib].resize(nFaces[ib]);
		fill(faces_bs[ib].begin(),faces_bs[ib].end(), array<int,6>{{-1,-1,-1,-1,-1,-1}});
	}

	voxelField<int3> ownerMapper(n.x+1,n.y+1,n.z+1,int3(-1,-1,-1));


	cout<<"collecting faces  "<<endl;


#define recordF_m( l10,l11,l20,l21,l30,l31,dir,ii,jj,kk,type )  {			\
	if (ownerMapper(ii,jj,kk)[dir]<0)  {											 \
		++iFaces[type] ;																		 \
		ownerMapper(ii,jj,kk)[dir]=iFaces[type] + iStartFaces[type];	              \
		faces_bs[type][iFaces[type]][4]=iCells;	/* cell number (for owners) */         \
		faces_bs[type][iFaces[type]][0]=point_mapper(ii,jj,kk);								  \
		faces_bs[type][iFaces[type]][1]=point_mapper(ii+l11,jj+l21,kk+l31);					\
		faces_bs[type][iFaces[type]][2]=point_mapper(ii+l10+l11,jj+l20+l21,kk+l30+l31);  \
		faces_bs[type][iFaces[type]][3]=point_mapper(ii+l10,jj+l20,kk+l30);					\
	}																			   \
	else {																		   \
		faces_bs[type][ownerMapper(ii,jj,kk)[dir]][5]=iCells; /* cell number (for neighbours) */ \
	}																			   \
  }


	#define  kclockwiseRecordF(type) recordF_m( 1,0,0,1,0,0, 2, ix-1,iy-1,iz-1, type)
	#define kuclockwiseRecordF(type) recordF_m( 0,1,1,0,0,0, 2, ix-1,iy-1,iz  , type)
	#define  jclockwiseRecordF(type) recordF_m( 0,1,0,0,1,0, 1, ix-1,iy-1,iz-1, type)
	#define juclockwiseRecordF(type) recordF_m( 1,0,0,0,0,1, 1, ix-1,iy  ,iz-1, type)
	#define  iclockwiseRecordF(type) recordF_m( 0,0,1,0,0,1, 0, ix-1,iy-1,iz-1, type)
	#define iuclockwiseRecordF(type) recordF_m( 0,0,0,1,1,0, 0, ix  ,iy-1,iz-1, type)




	int iCells=-1;


	array<int,255> iFaces; iFaces.fill(-1);

	for (int iz=1;iz<=n.z;iz++)  {  cout<<(iz%80 ? '.' : '\n');cout.flush();
	 for (int iy=1;iy<=n.y;iy++)
	  for (int ix=1;ix<=n.x;ix++)
		if (!vxlImg(ix,iy,iz))  {
			iCells++;

			unsigned char
			neiv=vxlImg(ix-1,iy,iz);
			if (ix!=1)  {
			  if (neiv)       iclockwiseRecordF(neiv)
			  else            iclockwiseRecordF(Internal)    }
			else if (neiv)   iclockwiseRecordF(neiv)
			else              iclockwiseRecordF(Left)

			neiv=vxlImg(ix+1,iy,iz);
			if (ix!=n.x)  {
			  if (neiv)       iuclockwiseRecordF(neiv)
			  else            iuclockwiseRecordF(Internal)   }
			else if (neiv)    iuclockwiseRecordF(neiv)
			else              iuclockwiseRecordF(Right)


			neiv=vxlImg(ix,iy-1,iz);
			if (iy!=1)  {
			  if (neiv)       jclockwiseRecordF(neiv)
			  else            jclockwiseRecordF(Internal)   }
			else if (neiv)    jclockwiseRecordF(neiv)
			else              jclockwiseRecordF(Bottom)

			neiv=vxlImg(ix,iy+1,iz);
			if (iy!=n.y)  {
			  if (neiv)       juclockwiseRecordF(neiv)
			  else            juclockwiseRecordF(Internal)  }
			else if (neiv)    juclockwiseRecordF(neiv)
			else              juclockwiseRecordF(Top)
			

			neiv=vxlImg(ix,iy,iz-1);
			if (iz!=1)  {
			  if (neiv)       kclockwiseRecordF(neiv)
			  else            kclockwiseRecordF(Internal)  }
			else if (neiv)    kclockwiseRecordF(neiv)
			else              kclockwiseRecordF(Back)

			neiv=vxlImg(ix,iy,iz+1);
			if (iz!=n.z)  {
			  if (neiv)       kuclockwiseRecordF(neiv)
			  else            kuclockwiseRecordF(Internal)   }
			else if (neiv)    kuclockwiseRecordF(neiv)
			else              kuclockwiseRecordF(Front)

		}
	}




	point_mapper.reset(0,0,0,0);



  {
	ofstream faces(Folder+"/faces");	assert(faces);
	faces<<
		"FoamFile\n"
		"{\n"
		"	version	 2.0;\n"
		"	format	  ascii;\n"
		"	class	   faceList;\n"
		"	location	\""<<Folder+"\";\n"
		"	object	  faces;\n"
		"}\n\n"<<sumnFaces<<"\n("<<endl;


	ofstream owner(Folder+"/owner");	assert(owner);
	owner<<
		"FoamFile\n"
		"{\n"
		"	version	 2.0;\n"
		"	format	  ascii;\n"
		"	class	   labelList;\n"
		"	location	\""<<Folder+"\";\n"
		"	object	  owner;\n"
		"}\n\n"<<sumnFaces<<"\n("<<endl;

	#define write_faces_owners(ipp) \
	 for (const auto& ff:faces_bs[ipp])  { \
		faces<<'('<<ff[0]<<' '<<ff[1]<<' '<<ff[2]<<' '<<ff[3]<<")\n"; \
		owner<<ff[4]<<"\n";  \
	}
	//write_faces_owners(Internal)
	for(int ib=0;ib<nVVs;++ib)  	write_faces_owners(ib)
	if (procIsijk(iPrc-1,jPrc,kPrc)<0)   write_faces_owners(Left)
	if (procIsijk(iPrc+1,jPrc,kPrc)<0)   write_faces_owners(Right)
	if (procIsijk(iPrc,jPrc-1,kPrc)<0)   write_faces_owners(Bottom)
	if (procIsijk(iPrc,jPrc+1,kPrc)<0)   write_faces_owners(Top)
	if (procIsijk(iPrc,jPrc,kPrc-1)<0)   write_faces_owners(Back)
	if (procIsijk(iPrc,jPrc,kPrc+1)<0)   write_faces_owners(Front)

	if (procIsijk(iPrc-1,jPrc,kPrc)>=0)  write_faces_owners(Left)
	if (procIsijk(iPrc+1,jPrc,kPrc)>=0)  write_faces_owners(Right)
	if (procIsijk(iPrc,jPrc-1,kPrc)>=0)  write_faces_owners(Bottom)
	if (procIsijk(iPrc,jPrc+1,kPrc)>=0)  write_faces_owners(Top)
	if (procIsijk(iPrc,jPrc,kPrc-1)>=0)  write_faces_owners(Back)
	if (procIsijk(iPrc,jPrc,kPrc+1)>=0)  write_faces_owners(Front)


	faces<<")"<<endl;
	owner<<")"<<endl;
  }



	ofstream neighbour(Folder+"/neighbour");	assert(neighbour);
	neighbour<<
		"FoamFile\n"
		"{\n"
		"	version	 2.0;\n"
		"	format	  ascii;\n"
		"	class	   labelList;\n"
		"	location	\""<<Folder+"\";\n"
		"	object	  neighbour;\n"
		"}\n\n"
		<<nFaces[Internal]<<"\n("<<endl;
	for (const auto& ff:faces_bs[Internal])  neighbour<<ff[5]<<"\n";
	neighbour<<")"<<endl;

	cout<<" :/ "<<endl;

}



