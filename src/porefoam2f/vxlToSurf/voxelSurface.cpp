// conversion code from voxel images to 3D surfaces

// Copyright (C) 2020  Ali Qaseminejad Raeini 

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

// For any queries, contact Ali Q. Raeini: email aliqasemi@gmail.com


#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>



#include "voxelImage.h"
#include "InputFile.h"
//#include "happly.h"



using namespace std;



#include "surfUtils.h"




inline int collectManifoldFaces(label pI, label connectingFace ,DynamicList<label> & group1, bool handlNonManiflEdges,
const facePiece& faces, const piece<point>& points, const labelListList& pPoints, const labelListList& pFaces )
{


	const bool selectMin=false;

	int nClctd=0;
	nClctd+=appendUnique(group1,connectingFace);

	const labelListList eFaces = edgeFaces(pPoints[pI],faces,pFaces[pI]);
	for_(pPoints[pI], eI)
	{

		const ints & myEFs=eFaces[eI];
		if (myEFs.size()==1)   Info<<"\n Error point "<<points[pI]<<" conected to edge of only one face   " <<points[pI]<<"  " <<points[pPoints[pI][eI]]<<"	face: " <<eFaces[eI] << "\n\n";

		else if (myEFs.size()==2)
		{
			if      (myEFs[0]==connectingFace) nClctd+=appendUnique(group1,myEFs[1]);
			else if (myEFs[1]==connectingFace) nClctd+=appendUnique(group1,myEFs[0]);
		}
		else if (myEFs.size()>2 && handlNonManiflEdges)
		{
			bool isMyEdge =false;
			for_(myEFs, fI)
				if( myEFs[fI]==connectingFace )
					isMyEdge=true;

			if(isMyEdge)
			{


				dbls closeNess(myEFs.size(), selectMin ? 1000.0: -1000.0 );
				dbl3 masterNormal=normal(faces[connectingFace],points);
				dbl3 Ce=0.5*(points[pI]+points[pPoints[pI][eI]]);
				dbl3 tmf=centre(faces[connectingFace],points)-Ce;
				tmf/=mag(tmf)+1.0e-15;

				for_(myEFs, fI) if ( !(myEFs[fI]==connectingFace) )
				{
					dbl3 tf=centre(faces[myEFs[fI]],points)-Ce;
					tf/=mag(tf)+1.0e-15;
					scalar sin=tf&masterNormal;
					scalar cos=tf&tmf;
					const double PI=3.14159265;
					double angle=std::atan2 (sin,cos) * 180 / PI;

					if ( angle<0.0) angle=360+angle;
					closeNess[fI]=angle;
				}

				if (!selectMin)        closeNess=-closeNess;

				label nei=std::distance(closeNess.begin(), std::min_element(closeNess.begin(), closeNess.end()));
				if (appendUnique(group1,myEFs[nei])==1)   ++nClctd;

			}
		}

	}

  return nClctd;
}


//=============================================================================================
//=============================================================================================
//=============================================================================================
//=============================================================================================


void correct( facePiece & faces, DynamicList<point> & points, bool handlNonManiflEdges )
{
	Info<<"	"<<points.size()<<"  points and  "<<faces.size()<<"  faces in surface, looking for errors:"<<endl;		/*Info.flush()*/;

	// new points and changed faces
	DynamicList<point> addedPoints; addedPoints.reserve(points.size()/100+1);
	//DynamicList<face>  modifiedFaces(faces.size()/100+1);
	std::map<label,face>  changedFaces;//(faces.size()/100+1);
	const auto& meshPoints = points;

	label nProblemPoints = 0;



	const labelListList pFaces  = pointFaces(points.size(),faces);
	const labelListList pPoints = getPointPoints(points.size(),faces);


	size_t nSkipped = 0;
	label iLastPoint=points.size()-1;
	for_(points, pI)
	{
		const ints& myFaces = pFaces[pI];
		if (myFaces.size()>5)
		{

			DynamicList<label> group1; group1.reserve(myFaces.size());

			group1.append(myFaces[0]);
			int nClctdTotal=1;

			//int gfI=0;
			{int nClctd; do { nClctd=0;
				 for_(group1, gfI)  nClctd += collectManifoldFaces(pI, group1[gfI] ,group1,handlNonManiflEdges, faces, points, pPoints, pFaces);
				 nClctdTotal+=nClctd;
			} while (nClctd>0);}
			for_(group1, gfI)  
			nClctdTotal+= collectManifoldFaces(pI, group1[gfI] ,group1,handlNonManiflEdges, faces, points, pPoints, pFaces);

			if(nClctdTotal<myFaces.size())
			{
				  bool PreviouslyModified=false;
				  for_(group1,gfI)   if (changedFaces.find(group1[gfI])!=changedFaces.end())	PreviouslyModified=true;

				
				if( (nClctdTotal<myFaces.size()-2) &&  (nClctdTotal>2) && !PreviouslyModified)
				{	++nProblemPoints;

					addedPoints.append(points[pI]); ++iLastPoint;  ///. clone the point, the clone will go to the end
					ensure(group1.size());
					for_(group1,gfI)
					{
						face modifedF=faces[group1[gfI]];        ///. get the face
						label index=distance(modifedF.begin(), find(modifedF.begin(),modifedF.end(),pI));

						auto modifOld=changedFaces.find(group1[gfI]);
						if (modifOld!=changedFaces.end())	modifedF=modifOld->second;

						if (index>=0)	modifedF[index]=iLastPoint; ///. change the face // TODO change inline
						else			Info<<index<<"Error in collecting connected faces : negative array index "<<endl;

						changedFaces.insert({group1[gfI],modifedF});
						//modifiedFaces.append( modifedF );  ///. save the changed face
					}
				}
				else  ++nSkipped; //Info<<"Point "<<pI<<" skipped, as this will cause singly connected edges,  collected " <<  nClctdTotal<<" faces out of "<<myFaces.size()<<endl;
			}

		}
		//else if(myFaces.size()<3)  Info<<myFaces.size()<<" "<<pI<<":wrong point:("<<points[pI]<<"), nf:"<<myFaces.size()<<" \n";


	}
	
	if(nSkipped)  Info<<"  "<<nSkipped<<" points skipped, as this will cause singly connected edges"<<endl;


	Info<< nProblemPoints<< " multiply shared (nonmanifold) edges,  ";
	Info<<"addedPoints: "<<addedPoints.size() <<"  changedFaces: "<<changedFaces.size() <<endl;


	//points.pback(piece<point>(&addedPoints[0], &addedPoints.back()));
	points.insert( points.end(), addedPoints.begin(), addedPoints.end() );


	for(const auto& chngdF:changedFaces)	faces[chngdF.first]= chngdF.second;

	Info<<faces.size()<<" faces "<<points.size()<<" points "<<endl;

}







void correctbioti( facePiece & faces, labelList& fMarks, dbl3s& points, int stage )
{
	Info<<points.size()<<" points,  "<<faces.size()<<" faces,  stage:"<<stage<<"   correcting: ";  cout.flush();

	std::set<label>  deletedFaceIs;//(faces.size()/100+1);

	label nProblemPoints = 0;
	label nbads = 0;


	facePiece & Sfaces=faces;
	ints pMarks(points.size(),0);
	std::vector<labelList> pFaces(points.size());
	for_(Sfaces, fI)
	{	const face& f = Sfaces[fI];
		for_(f,i) ++pMarks[f[i]];
	}
	for_(pFaces, pI)	pFaces[pI].resize(pMarks[pI]);
	pMarks=-1;
	for_(Sfaces, fI)
	{	const face& f = Sfaces[fI];
		for_(f,i) pFaces[f[i]][++pMarks[f[i]]]=fI;
	}

	labelList pointPointmap(points.size(),-1);




	for_(pFaces, pI)
	{
		label pII = pI;///surf1.meshPoints()[pI];
		label n3fs(0),delFI(-1);
		const labelList& myFaces = pFaces[pI];
		if (myFaces.size()>4)
		{
			label myMark(fMarks[myFaces[0]]); bool CLine(false);
			for_(myFaces,i)
			{
				const face& f=Sfaces[myFaces[i]];
				auto fItr=find(f.begin(), f.end(),pII);  //find(list.begin(), list.end(),val)
				if (fItr==f.end())
					++nbads;
				else
					if (pFaces[*nextCircIter(f,fItr)].size() ==3 && pFaces[*prevCircIter(f,fItr)].size() ==3) {++n3fs; delFI=myFaces[i];}
					if(myMark!=fMarks[myFaces[i]]) CLine=true;
			}
			if (CLine && n3fs==0 && stage==2)
			{
			  for_(myFaces,i)
			  {
				const face& f=Sfaces[myFaces[i]];
				auto fItr=find(f.begin(), f.end(),pII);  //find(list.begin(), list.end(),val)
				if (fItr==f.end())
				{
					++nbads;//Info<<"  bad:"<<pI<<" "<<pII<<"  ";
				}
				else
				{
					auto opospItr= nextCircIter(f,nextCircIter(f,fItr)); //f[(mePIinf+2)%f.size()];
					if (pFaces[*opospItr].size() ==3 && pFaces[*prevCircIter(f, fItr)].size()>6 && pFaces[*nextCircIter(f,fItr)].size()>6)
					{
						++n3fs; delFI=myFaces[i]; pII=*prevCircIter(f, fItr);
						//~ Info<<"!"<<(mePIinf+2)%f.size()<<" "<<delFI<<"    "<<pI<<" "<<opospI<<" "<<pII<<" "<<endl;
						//~ exit(-1);
					}
				}
			  }
			}
		}

		if ( (n3fs==1 ||  (stage>=1 && n3fs>=1)) && (deletedFaceIs.find(delFI)==deletedFaceIs.end())  && pointPointmap[pII]==-1  && pointPointmap[pI]==-1 )
		{
			nProblemPoints++;

			const face& f=Sfaces[delFI];
			auto fItr=find(f.begin(), f.end(),pII);  //find(list.begin(), list.end(),val)
			if (fItr==f.end())
				Info<<"  bad:"<<pI<<" "<<pII<<"  ";
			else
			{
				int keptP,delP;
				if(*nextCircIter(f, fItr) < *prevCircIter(f, fItr))
					{keptP=*nextCircIter(f, fItr);  delP=*prevCircIter(f, fItr);}
				else{delP=*nextCircIter(f, fItr);  keptP=*prevCircIter(f, fItr);}
				
				if (pointPointmap[delP]==-1 && pointPointmap[keptP]==-1)
				{
					pointPointmap[delP]=keptP;
					deletedFaceIs.insert(delFI);
					//~ mergedPointIndices.append(keptP);
					//~ mergedPointIndices.append(delP);
				}
			}
		}
		//~ else if(myFaces.size()<3)  Info<<pI<<": wrong point : "<<points[pI]<<endl;

	}



	Info<< nProblemPoints<< " dbl 3-nei points in face. bads:"<<nbads<<" *  ";


	{
		label iLastPoint=-1;
		for_(points,i)
		{
			if (pointPointmap[i]<0)
			{
				points[++iLastPoint]=points[i];
				pointPointmap[i]=iLastPoint;
			}
			else
			{
				if (pointPointmap[i]>=i) Info<<" Errorsddsfdf "<<endl;
				pointPointmap[i]=pointPointmap[pointPointmap[i]];
				points[pointPointmap[i]]=0.5*(points[i]+points[pointPointmap[i]]);
			}
		}
		points.resize(iLastPoint+1);
	}
	Info<<"    "<<points.size()<<" points "<<" ";
	{
		for(auto delF:deletedFaceIs)	faces[ delF ]=face({-1,-1,-1,-1});
		label iLastFace=-1;
		for_(Sfaces,i)
		{
			if(faces[i].size())
			{
				face f=faces[i];
				for_(f,ii)  f[ii]=pointPointmap[f[ii]];
				faces[++iLastFace]=f;
				fMarks[iLastFace]=fMarks[i];
			}
		}
		//faces.resize(iLastFace+1);
		//fMarks.resize(iLastFace+1); 
		int Errr, Fixme;
	}
   Info<<faces.size()<<" faces "<<endl;


}






surfMsh createSurface(InputFile& meshingDict, const voxelImage & vxlImg, const int nVVs, const std::string& outputSurface)
{
	surfMsh srfmsh;
	facePieceList& faces_bs = srfmsh.faces_bs;
	DynamicField<point>& points = srfmsh.points;

	Info<<__FUNCTION__<<": "<<endl;
	int3 n=vxlImg.size3();n[0]-=2;n[1]-=2;n[2]-=2;
	dbl3 X0=vxlImg.X0();
	dbl3 dx=vxlImg.dx();
	X0+=dx;


	std::array<size_t,256> nFaces; nFaces.fill(0);

	for (int iz=1;iz<=n[2];++iz)
	 for (int iy=1;iy<=n[1];++iy)
	  for (int ix=1;ix<=n[0];++ix)
	  {
		//if (!vxlImg(ix,iy,iz))
	    {
			unsigned char vv=vxlImg(ix,iy,iz);
			unsigned char
			neiv=vxlImg(ix-1,iy,iz);
			if (neiv>vv) ++nFaces[vv];

			neiv=vxlImg(ix+1,iy,iz);
			if (neiv>vv) ++nFaces[vv];

			neiv=vxlImg(ix,iy-1,iz);
			if (neiv>vv) ++nFaces[vv];

			neiv=vxlImg(ix,iy+1,iz);
			if (neiv>vv) ++nFaces[vv];

			neiv=vxlImg(ix,iy,iz-1);
			if (neiv>vv) ++nFaces[vv];

			neiv=vxlImg(ix,iy,iz+1);
			if (neiv>vv) ++nFaces[vv];
	    }
	  }

	for(int ib=0;ib<256;++ib)  if(nFaces[ib])  cout<< " Faces_"<<ib<<" "<<nFaces[ib]<<endl;
	cout<<endl;




	cout<<"creating faces"<<endl;

	{ 	size_t nAll=0; 
		for(int ib=0;ib<256;++ib)  nAll+=nFaces[ib]; 
		srfmsh.faces.resize(nAll); 
		nAll=0; 
		for(int ib=0;ib<256;++ib)  { faces_bs[ib].reset(facePiece(&srfmsh.faces[nAll], nFaces[ib]));  nAll+=nFaces[ib]; }
	}

	voxelField<int> point_mapper(n[0]+1,n[1]+1,n[2]+1,-1);




	cout<<"collecting faces"<<endl;
	int iPoints=-1;
	auto point_mapper_insert = [&point_mapper, &points,  &dx,  &X0, &iPoints](int iix, int iiy, int iiz) ->
	 int {
		if (point_mapper(iix,iiy,iiz)<0)
		{
			points.append(point(dx[0]*iix+X0[0],dx[1]*iiy+X0[1],dx[2]*iiz+X0[2]));
			point_mapper(iix,iiy,iiz)=++iPoints;
			return	iPoints;
		}
		return point_mapper(iix,iiy,iiz);
	 };


	#define recordF_m( l10,l11,l20,l21,l30,l31, ii,jj,kk,type )					  \
	  {                                                                                  \
		   ++iFaces[type] ;														\
			/*faces_bs[type][iFaces[type]].setSize(4+1);*/								  \
			faces_bs[type][iFaces[type]][0]=point_mapper_insert(ii,jj,kk);								  \
			faces_bs[type][iFaces[type]][1]=point_mapper_insert(ii+l11,jj+l21,kk+l31);					\
			faces_bs[type][iFaces[type]][2]=point_mapper_insert(ii+l10+l11,jj+l20+l21,kk+l30+l31);  \
			faces_bs[type][iFaces[type]][3]=point_mapper_insert(ii+l10,jj+l20,kk+l30);					\
			faces_bs[type][iFaces[type]].zone=neiv;					\
	  }


	#define  kclockwiserecordF(type) recordF_m( 1,0,0,1,0,0,   ix-1,iy-1,iz-1, type)
	#define kuclockwiserecordF(type) recordF_m( 0,1,1,0,0,0,   ix-1,iy-1,iz  , type)
	#define  jclockwiserecordF(type) recordF_m( 0,1,0,0,1,0,   ix-1,iy-1,iz-1, type)
	#define juclockwiserecordF(type) recordF_m( 1,0,0,0,0,1,   ix-1,iy  ,iz-1, type)
	#define  iclockwiserecordF(type) recordF_m( 0,0,1,0,0,1,   ix-1,iy-1,iz-1, type)
	#define iuclockwiserecordF(type) recordF_m( 0,0,0,1,1,0,   ix  ,iy-1,iz-1, type)





	std::array<int,256> iFaces; iFaces.fill(-1);

	for (int iz=1;iz<=n[2];iz++)
	{  cout<<(iz%100 ? '.' : '\n');cout.flush();
	  for (int iy=1;iy<=n[1];iy++)
	    for (int ix=1;ix<=n[0];ix++)
		  //if (!vxlImg(ix,iy,iz))
		  {
			unsigned char vv=vxlImg(ix,iy,iz);

				unsigned char
				neiv=vxlImg(ix-1,iy,iz);
				if (neiv>vv)  {iclockwiserecordF( vv)}

				neiv=vxlImg(ix+1,iy,iz);
				if (neiv>vv)  {iuclockwiserecordF( vv)}


				neiv=vxlImg(ix,iy-1,iz);
				if (neiv>vv)  {jclockwiserecordF( vv)}

				neiv=vxlImg(ix,iy+1,iz);
				if (neiv>vv)  {juclockwiserecordF( vv)}


				neiv=vxlImg(ix,iy,iz-1);
				if (neiv>vv)  {kclockwiserecordF( vv)}

				neiv=vxlImg(ix,iy,iz+1);
				if (neiv>vv)  {kuclockwiserecordF( vv)}

		}
	}


//______________________________________________

	point_mapper.reset(0,0,0,0);


	Info<<"nPoints: "<<points.size()<<"  "<<endl;

	if (meshingDict.getOr(string("true"),"correctSurface")[0]=='t')
	 for(auto& facezs:faces_bs) if(facezs.size())
	 {
		correct(facezs,points,true);
		correct(facezs,points,true);
		correct(facezs,points,true);
		correct(facezs,points,true);
		correct(facezs,points,true);
		correct(facezs,points,false);
		correct(facezs,points,true);
		correct(facezs,points,false);
		correct(facezs,points,false);
		correct(facezs,points,false);
	 }

	return srfmsh;

}




