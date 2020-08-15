/*-------------------------------------------------------------------------*\
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

//! Description:
//!   creates a surface between the pore and the solid from a 3D rock image


    #include <fstream>
    #include <iostream>
    #include <vector>

    #include <assert.h>

#include "fvCFD.H"

#include "argList.H"
#include "timeSelector.H"
#include "graph.H"
#include "mathematicalConstants.H"
#include "labelVector.H"

#include "OFstream.H"
#include "triFaceList.H"
#include "triSurface.H"

#include "DynamicField.H"
#include "MeshedSurfaces.H"

#include "voxelImage.h"
#include "createSurface.h"

using namespace Foam;

//----------------------------------------------------------------------
// Main program:



int main(int argc, char *argv[])
{





#   include "setRootCase.H"
#   include "createTime.H"




    IOdictionary meshingDict
    (
        IOobject
        (
            "meshingDict",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );
    word headerName(meshingDict.lookup("headerName"));
    //word inputName(meshingDict.lookup("inputName"));
    word outputName(meshingDict.lookup("outputName"));
    word outputSurface(meshingDict.lookup("outputSurface"));
    //word trimName(meshingDict.lookup("trimName"));




	voxelImage vximage(headerName);
	vximage.threshold101(0,0);
	vximage.rotate('z');

	int3 n=vximage.size3();


	label nAddLyrs(readLabel(meshingDict.lookup("nAddLayers"))) ;
	label nCopyLyrsX(readLabel(meshingDict.lookup("nCopyLayersX"))) ;
	label nCopyLyrsYZ(readLabel(meshingDict.lookup("nCopyLayersYZ"))) ;
	Info<<"nCopyLayersX: "<<nCopyLyrsX<<endl;
	Info<<"nCopyLayersYZ: "<<nCopyLyrsYZ<<endl;

   vximage.cropD(int3(nCopyLyrsYZ, nCopyLyrsYZ,nCopyLyrsX),n-int3(nCopyLyrsYZ,nCopyLyrsYZ,nCopyLyrsX));//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX
	vximage.growBox(nAddLyrs +nCopyLyrsYZ+nCopyLyrsX);

	labelVector nSmoothBegin(meshingDict.lookup("nSmoothBegin")) ;
	labelVector nSmoothEnd((meshingDict.lookup("nSmoothEnd"))) ;

	n = vximage.size3();

	vximage.cropD(0,n, 2,1);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX

	voxelImage SmoothedVoxels=vximage;

	vximage.cropD( int3(&nSmoothBegin[0])+2 ,n-int3(&nSmoothEnd[0])-1, 0, 0 );    


	scalar perforationDepth(readScalar(meshingDict.lookup("perforationDepth"))) ;
	if (perforationDepth>0.00000001)
    {

		scalar perforationWidth(readScalar(meshingDict.lookup("perforationWidth"))) ;
		scalar perforationLengthFraction(readScalar(meshingDict.lookup("perforationLengthFraction"))) ;
		scalar perforationLengthFraction0=perforationLengthFraction-(1.0*nSmoothBegin[0]/n[0]);
		scalar perforationLengthFraction1=perforationLengthFraction-(1.0*nSmoothBegin[1]/n[1]);
		//vector IntactFaction=( vector(1.0,1.0,1.0) - smoothFaction  )/2.0;
		//scalar perforationWidth=(perforationWidth+perforationWidth*n[1])/2.0;
		label beginSkip_0=(1.0-perforationLengthFraction0)*0.5*n[0];
		label beginSkip_1=(1.0-perforationLengthFraction1)*0.5*n[1];
		label begin2_0=n[0]-beginSkip_0-perforationWidth;
		label begin2_1=n[1]-beginSkip_1-perforationWidth;
		
 
		voxelImage tmp1(perforationLengthFraction0*n[0],perforationWidth,perforationDepth,0);
		SmoothedVoxels.setBlock	( 1+beginSkip_0, 1+beginSkip_1, 2, tmp1 );
		SmoothedVoxels.setBlock	( 1+beginSkip_0, 1+begin2_1, 2, tmp1 );
		
		voxelImage tmp2(perforationWidth,perforationLengthFraction1*n[1],perforationDepth,0);
		SmoothedVoxels.setBlock	( 1+beginSkip_0, 1+beginSkip_1, 2, tmp2 );
		SmoothedVoxels.setBlock	( 1+begin2_0, 1+beginSkip_1, 2, tmp2 );

	}
	label nSmoothing(   readLabel(meshingDict.lookup("nSmoothing")  ) );
	Info<<"\nsmoothing  outside, nIterations:" << nSmoothing<<endl;

	for ( label i=0; i<nSmoothing ; i++ )
	{
		SmoothedVoxels.PointMedian032(18,18,0,1);
		SmoothedVoxels.FaceMedian06(2,4);
		SmoothedVoxels.setBlock(nSmoothBegin[0]+2, nSmoothBegin[1]+2, nSmoothBegin[2]+2, vximage);
	}
	vximage.reset(0,0,0,0);

	Info<<"Smoothing internals"<<endl;

	SmoothedVoxels.PointMedian032(20,20,0,1);
	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(22,22,0,1);
	SmoothedVoxels.FaceMedian06(1,5);
	SmoothedVoxels.FaceMedian06(1,5);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(1,5);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);


	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(22,22,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(1,5);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.PointMedian032(22,22,0,1);
	SmoothedVoxels.PointMedian032(23,23,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(1,5);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);

	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);

	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);


	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);


	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);

	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);

	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);

	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);


	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);




	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);



	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);



	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);



	SmoothedVoxels.PointMedian032(21,21,0,1);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);
	SmoothedVoxels.FaceMedian06(2,4);



	SmoothedVoxels.rotate('z');///. rotate back, see above
	n = SmoothedVoxels.size3();
	SmoothedVoxels.printInfo();



	writeSTLBINARY(SmoothedVoxels, outputSurface);//         XXXXXXXXXXXXXXXXXXXXXXXXXXXX

	int nAddedX =nAddLyrs + nCopyLyrsYZ+2;
	int nAddedYZ=nAddLyrs + nCopyLyrsX+2;
	SmoothedVoxels.cropD({nAddedX,nAddedYZ,nAddedYZ},n-int3(nAddedX,nAddedYZ,nAddedYZ));
	SmoothedVoxels.write(outputName);
	SmoothedVoxels.writeAConnectedPoreVoxel(headerName+"_aPoreVoxel");

	Info<<"finished  writeSTLBINARY, outputFileName: "<<outputSurface<<endl;

	Info<< "end" << endl;

	return 0;
}












int appendUnique(DynamicList<label>& dynList, label value)
{
	forAll(dynList,i) if(dynList[i]==value) return 0;
	dynList.append(value);
	return 1;
}

int getMinPos(const UList<scalar>& array)
{
	scalar min=array[0];
	int minPos=0;
	for (int i=1; i<array.size(); i++)   if(min>array[i])  { min=array[i];  minPos=i; }
	return minPos;
}

int findPos(const UList<label>& list,label value)
{
	forAll(list,i)	if (value==list[i])	return i;
	return -1;
}


inline int collectManifoldFaces(label meshPI, label connectingFace ,DynamicList<label> & group1, const meshedSurface& surf1, bool handlemultipliConnectedEdges)
{

	const labelListList& pEdges = surf1.pointEdges();
	const labelListList& eFaces = surf1.edgeFaces();
	const List<edge>& edges = surf1.edges();
	const List<face>& faces = surf1.surfFaces();
	const pointField& points = surf1.points();
	const labelList& meshPoints = surf1.meshPoints();

	const bool selectMin=false;

	int nCollected=0;
	nCollected+=appendUnique(group1,connectingFace);

	const labelList & myEdges=pEdges[meshPI];
	forAll(myEdges, myeI)
	{
		label eI=myEdges[myeI];

		const labelList & myEFs=eFaces[eI];
		if (myEFs.size()==1)   Info<<"\n Error point "<<meshPoints[meshPI]<<" conected to edge of only one face   " <<meshPoints[edges[eI][0]]<<"  " <<meshPoints[edges[eI][1]]<<"	face: " <<eFaces[eI] << "\n\n";

		else if (myEFs.size()==2)
		{
			if	  (myEFs[0]==connectingFace) nCollected+=appendUnique(group1,myEFs[1]);
			else if (myEFs[1]==connectingFace) nCollected+=appendUnique(group1,myEFs[0]);
		}
		else if (myEFs.size()>2 && handlemultipliConnectedEdges)
		{
			bool isMyEdge =false;
			forAll(myEFs, fI)
				if( myEFs[fI]==connectingFace )
				{
					isMyEdge=true;
				}

			if(isMyEdge)
			{


				SortableList<scalar> closeNess(myEFs.size(), selectMin ? 1000.0: -1000.0 );
				vector masterNormal=faces[connectingFace].normal(points);
				vector Ce=0.5*(points[meshPoints[edges[eI][0]]]+points[meshPoints[edges[eI][1]]]);
				vector tmf=faces[connectingFace].centre(points)-Ce;
				tmf/=mag(tmf)+1.0e-15;

				forAll(myEFs, fI) if ( !(myEFs[fI]==connectingFace) )
				{
					vector tf=faces[myEFs[fI]].centre(points)-Ce;
					tf/=mag(tf)+1.0e-15;
					scalar sin=tf&masterNormal;
					scalar cos=tf&tmf;
					const double PI=3.14159265;
					double angle=std::atan2 (sin,cos) * 180 / PI;

					if ( angle<0.0) angle=360+angle;
					closeNess[fI]=angle;
				}

				if (!selectMin)		closeNess=-closeNess;

				label nei=getMinPos(closeNess);
				if (appendUnique(group1,myEFs[nei])==1)   nCollected++;

			}
		}

	}

	return nCollected;
}


//=============================================================================================
//=============================================================================================
//=============================================================================================
//=============================================================================================


void correct( faceList & faces, DynamicField<point> & points, bool handlemultipliConnectedEdges )
{
	Info<<"	"<<points.size()<<"  points and  "<<faces.size()<<"  faces in surface, looking for errors:"<<endl;		/*Info.flush()*/;

	// new points and changed faces
	DynamicList<point> addedPoints(points.size()/100+1);
	DynamicList<face>  modifiedFaces(faces.size()/100+1);
	DynamicList<label>  modifiedFacesIndices(faces.size()/100+1);

	label nProblemPoints = 0;


	Xfer<List<face> > facesFer(faces,false);
	//Field<point> & pointsSF=points;
	Xfer<Field<point> > pointsFer(points,false);
	MeshedSurface<face> surf1(pointsFer, facesFer);

	const labelListList& pFaces = surf1.pointFaces();
	const labelList& meshPoints = surf1.meshPoints();

	const List<face>& Sfaces = surf1.surfFaces();



	label iLastPoint=points.size()-1;
	forAll(meshPoints, meshPI)
	{
		label pI=(meshPoints[meshPI]);
		const labelList& myFaces = pFaces[meshPI];
		if (myFaces.size()>5)
		{

			DynamicList<label> group1(myFaces.size());

			group1.append(myFaces[0]);
			int nCollectedTotal=1;


			int nCollected;
			do { nCollected=0;
				 forAll(group1, gfI)	nCollected += collectManifoldFaces(meshPI, group1[gfI] ,group1,surf1,handlemultipliConnectedEdges);
				 nCollectedTotal+=nCollected;
			} while (nCollected>0);

			forAll(group1, gfI)   nCollectedTotal+= collectManifoldFaces(meshPI, group1[gfI] ,group1,surf1,handlemultipliConnectedEdges);

			if(nCollectedTotal<myFaces.size())
			{
				bool PreviouslyModified=false;
				forAll(group1,gfI)   if (findPos(modifiedFacesIndices,group1[gfI])>=0)	PreviouslyModified=true;

				
				if( (nCollectedTotal<myFaces.size()-2) &&  (nCollectedTotal>2) && !PreviouslyModified)
				{   nProblemPoints++;

					addedPoints.append(points[pI]); iLastPoint++;  ///. duplicate the point, it will go to the end

					forAll(group1,gfI)
					{
						face modifiedFace=Sfaces[group1[gfI]];		///. get the face

						label iModFace=findPos(modifiedFacesIndices,group1[gfI]);
						if (iModFace>=0)	modifiedFace=modifiedFaces[iModFace];

						modifiedFacesIndices.append(group1[gfI]);
						label index=findPos(Sfaces[group1[gfI]],pI);
						if (index>=0)	modifiedFace[index]=iLastPoint; ///. change the face
						else			Info<<index<<"Error in collecting connected faces : negative array index "<<endl;

						modifiedFaces.append( modifiedFace );  ///. save the changed face
					}

				}
				else  Info<<"Point "<<pI<<" skipped, as this will cause singly connected edges,  collected " <<  nCollectedTotal<<" faces out of "<<myFaces.size()<<endl;
			}

		}
		else if(myFaces.size()<3)  Info<<myFaces.size()<<" "<<pI<<":wrong point:("<<points[pI]<<")"<<myFaces.size()<<" \n";


	}



	Info<<"found "<< nProblemPoints<< " problem points - points shared by edges that connected to more than 2 face"<<endl;
	Info<<"addedPoints: "<<addedPoints.size() <<"  modifiedFacesIndices: "<<modifiedFacesIndices.size() <<"  modifiedFaces: "<<modifiedFaces.size()<<endl;

	points.append(addedPoints);
	forAll(modifiedFacesIndices,i)	faces[ modifiedFacesIndices[i] ]= modifiedFaces[i];

	Info<<faces.size()<<" faces "<<points.size()<<" points "<<endl;


}


//----------------------------------------------------------------------
