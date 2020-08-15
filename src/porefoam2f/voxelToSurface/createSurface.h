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
Ali Q Raeini:	a.q.raeini@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk

 Description:
	creates a surface between the pore and the solid from a 3D rock image
\*-------------------------------------------------------------------------*/
#include <array>

void correct( faceList & faces, DynamicField<point> & points, bool handlemultipliConnectedEdges  );



void writeSTLBINARY( const voxelImage & vxlImg, std::string outputSurface)
{


	Info<<"writeSTLBINARY: "<<endl;
	int3 n=vxlImg.size3();n[0]-=2;n[1]-=2;n[2]-=2;
	dbl3 X0=vxlImg.X0(); 
	dbl3 dx=vxlImg.dx();
	X0+=dx;



//=======================================================================



	const int nVVs=2;
	const int

	Internal     = 0, 
	Grainwalls   = 1, 
	Left  =nVVs+0,
	Right =nVVs+1,
	Bottom=nVVs+2,
	Top   =nVVs+3,
	Back  =nVVs+4,
	Front =nVVs+5;

	std::array<std::string,255> B_nams;
	B_nams[Internal]="Internal";
	B_nams[Grainwalls]="Grainwalls";
	int nBoundaries=5+nVVs;
	for(int ib=2;ib<nVVs;++ib)  {B_nams[ib]="VV"+_s(ib)+"B"; }
	B_nams[Left]="Left";
	B_nams[Right]="Right";
	B_nams[Bottom]="Bottom";
	B_nams[Top]="Top";
	B_nams[Back]="Back";
	B_nams[Front]="Front";

	std::array<size_t,255> nFaces; nFaces.fill(0);


	for (int iz=1;iz<=n[2];iz++)
	 for (int iy=1;iy<=n[1];iy++)
	  for (int ix=1;ix<=n[0];ix++)
		{
			if (!vxlImg(ix,iy,iz))
			{
				//nCells++;

				unsigned char 
				neiv=vxlImg(ix-1,iy,iz);
				if (ix!=1)
				{
					if (neiv) ++nFaces[neiv];
				}else if (neiv) ++nFaces[neiv];

				neiv=vxlImg(ix+1,iy,iz);
				if (ix!=n[0])
				{
				  if (neiv)     ++nFaces[neiv];
				}else if (neiv) ++nFaces[neiv];

				neiv=vxlImg(ix,iy-1,iz);
				if (iy!=1)
				{
				  if (neiv)     ++nFaces[neiv];
				}else if (neiv) ++nFaces[neiv];

				neiv=vxlImg(ix,iy+1,iz);
				if (iy!=n[1])
				{
				  if (neiv)     ++nFaces[neiv];
				}else if (neiv) ++nFaces[neiv];

				neiv=vxlImg(ix,iy,iz-1);
				if (iz!=1)
				{
				  if (neiv)     ++nFaces[neiv];
				}else if (neiv) ++nFaces[neiv];

				neiv=vxlImg(ix,iy,iz+1);
				if (iz!=n[2])
				{
				  if (neiv)     ++nFaces[neiv];
				}else if (neiv) ++nFaces[neiv];
		}
	}

	//cout<<"nCells: "<<nCells<<endl;
	for(int ib=0;ib<255;++ib)  if(nFaces[ib])  cout<< " Faces_"<<ib<<" "<<nFaces[ib]<<endl;


	nBoundaries+=	( 1 && nFaces[Left])+ 
						( 1 && nFaces[Right])+
						( 1 && nFaces[Bottom])+ 
						( 1 && nFaces[Top])+
						( 1 && nFaces[Back])+ 
						( 1 && nFaces[Front]);


	std::array<size_t,255> iStartFaces; iStartFaces.fill(0);


	cout<<"creating faces"<<endl;
 

	std::array<faceList,255> faces_bs;
	size_t sumnFaces=0;
	for(int ib=0;ib<255;++ib)  if(nFaces[ib])
	{
		sumnFaces+=nFaces[ib];
		faces_bs[ib].resize(nFaces[ib]);
		//std::fill(faces_bs[ib].begin(),faces_bs[ib].end(), face(4));
	}

voxelField<int> point_mapper(n[0]+1,n[1]+1,n[2]+1,-1);

	DynamicField<point> points;



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
  {																				  \
		   ++iFaces[type] ;														\
			faces_bs[type][iFaces[type]].setSize(4);								  \
			faces_bs[type][iFaces[type]][0]=point_mapper_insert(ii,jj,kk);								  \
			faces_bs[type][iFaces[type]][1]=point_mapper_insert(ii+l11,jj+l21,kk+l31);					\
			faces_bs[type][iFaces[type]][2]=point_mapper_insert(ii+l10+l11,jj+l20+l21,kk+l30+l31);  \
			faces_bs[type][iFaces[type]][3]=point_mapper_insert(ii+l10,jj+l20,kk+l30);					\
 }


#define  kclockwiserecordF(type) recordF_m( 1,0,0,1,0,0,   ix-1,iy-1,iz-1, type)
#define kuclockwiserecordF(type) recordF_m( 0,1,1,0,0,0,   ix-1,iy-1,iz  , type)
#define  jclockwiserecordF(type) recordF_m( 0,1,0,0,1,0,   ix-1,iy-1,iz-1, type)
#define juclockwiserecordF(type) recordF_m( 1,0,0,0,0,1,   ix-1,iy  ,iz-1, type)
#define  iclockwiserecordF(type) recordF_m( 0,0,1,0,0,1,   ix-1,iy-1,iz-1, type)
#define iuclockwiserecordF(type) recordF_m( 0,0,0,1,1,0,   ix  ,iy-1,iz-1, type)




	int iCells=-1;


	std::array<int,255> iFaces; iFaces.fill(-1);

	for (int iz=1;iz<=n[2];iz++)
	{  cout<<(iz%50 ? '.' : '\n');cout.flush();
		for (int iy=1;iy<=n[1];iy++)
		 for (int ix=1;ix<=n[0];ix++)
		  if (!vxlImg(ix,iy,iz))
			{

				iCells++;

				unsigned char 
				neiv=vxlImg(ix-1,iy,iz);
				if (ix!=1)
				{
				  if (neiv)    {iclockwiserecordF( neiv)}
				}else if (neiv) {iclockwiserecordF( neiv)}

				neiv=vxlImg(ix+1,iy,iz);
				if (ix!=n[0])
				{
				  if (neiv)    {iuclockwiserecordF( neiv)}
				}else if (neiv) {iuclockwiserecordF( neiv)}


				neiv=vxlImg(ix,iy-1,iz);
				if (iy!=1)
				{
				  if (neiv)    {jclockwiserecordF( neiv)}
				}else if (neiv) {jclockwiserecordF( neiv)}

				neiv=vxlImg(ix,iy+1,iz);
				if (iy!=n[1])
				{
				  if (neiv)    {juclockwiserecordF( neiv)} 
				}else if (neiv) {juclockwiserecordF( neiv)} 


				neiv=vxlImg(ix,iy,iz-1);
				if (iz!=1)
				{
				  if (neiv)    {kclockwiserecordF(neiv)}
				}else if (neiv) {kclockwiserecordF( neiv)}

				neiv=vxlImg(ix,iy,iz+1);
				if (iz!=n[2])
				{ 
				  if (neiv)    {kuclockwiserecordF( neiv)}
				}else if (neiv) {kuclockwiserecordF( neiv)}


		}
	}


//______________________________________________

point_mapper.reset(0,0,0,0);


	Info<<"nPoints: "<<points.size()<<"	  "<<endl;		/*Info.flush()*/;



	{

		correct(faces_bs[1],points,true);
		correct(faces_bs[1],points,true);
		correct(faces_bs[1],points,false);
		correct(faces_bs[1],points,false);
		correct(faces_bs[1],points,false);
		correct(faces_bs[1],points,false);
	}
	{
		Xfer<List<face> > facesFer(faces_bs[1],true);
		Field<point> & pointsSF=points;
		Xfer<Field<point> > pointsFer(pointsSF,true);
		meshedSurface surf1(pointsFer, facesFer);

		//surf.triangulate ();

		surf1.write(outputSurface);
	}

}

 
