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
//!   Add a face layer between internal faces and boundary faces such that
//!     no two internal faces touch a single boundary edge, retired, thanks to cfMesh



#include <fstream>
#include <cassert>
#include <vector>

#include "argList.H"
#include "timeSelector.H"
#include "Time.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "transformGeometricField.H"
#include "IOobjectList.H"


#include "syncTools.H"

using namespace Foam;


using Foam::label;
using Foam::labelList;
using Foam::Info;
using Foam::face;
using Foam::edge;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    //argList::addNote
    //(
        //"adds a face layer between internal faces and boundary ensuring each face connects to only one boundary face."
    //);

    //Foam::timeSelector::addOptions();

    //argList::validArgs.append("relax");
    //argList::validArgs.append("nIter");
    //argList::validArgs.append("strictOnBoundary");
//
    //scalar relax = args.argRead<scalar>(1);
//
    //label nIter = args.argRead<label>(2);
    //bool strictOnBoundary = args.argRead<bool>(3);
//
	//Info<<"  relax:"<<relax <<"  nIter:"<<nIter <<endl;


#   include "setRootCase.H"
#   include "createTime.H"





#   include "createMesh.H"





	int nInFaces = mesh.nInternalFaces();
	Foam::faceList    faces = mesh.faces();
	const Foam::pointField &   points = mesh.points();
	const Foam::pointField &   CCntres = mesh.C();

	const Foam::labelListList &   edgesFs = mesh.edgeFaces();
	const Foam::labelListList &   edgesCs = mesh.edgeCells();
	label iGrainWalls = mesh.boundaryMesh().findPatchID("Grainwalls");
	const Foam::polyPatch& gwPatch= mesh.boundaryMesh()[iGrainWalls];
	const Foam::edgeList&  	edges = mesh.edges();




	const Foam::labelListList & fEdgesStupidOpenFOAM = mesh.faceEdges();
	Foam::labelListList fEdges = fEdgesStupidOpenFOAM;
	forAll(faces, fI)
	{
		const face& fps = faces[fI];
		const labelList& fesStupid = fEdgesStupidOpenFOAM[fI];
		labelList& fes = fEdges[fI];
		forAll(fps,fpi)
		{
			label p0=fps[fpi];
			label p1=fps.nextLabel(fpi);
			bool found = false;
			forAll(fesStupid, fei) ///. find edge
			{
				const edge& e = edges[fesStupid[fei]];
				if ((e[0] == p0 && e[1] == p1) || (e[0] == p1 && e[1] == p0))
				{
					 found=true;   fes[fpi] = fesStupid[fei];   break;
				}
			}
			if (!found) Info<<"badege"<<Foam::endl;
			if (fes[fpi]!=fesStupid[fpi]) Info<<"fixtege"<<Foam::endl;
		}
	}

	Foam::labelList     owns = mesh.faceOwner();
	Foam::labelList     neis = mesh.neighbour();



	//forAll(edges, ei) { edges[ei].start()=grainpts[edges[ei].start()];  edges[ei].end()=grainpts[edges[ei].end()]; }

	labelList grainpts = gwPatch.meshPoints();
	labelList pointLbls(points.size(), 0);
	forAll(grainpts, ip)		pointLbls[grainpts[ip]]=2^1;

	//forAll(edges, ie)		e1s[ie]=edges[ie].fi

	//labelList meshEdges(gwPatch.meshEdges(mesh.edges(), mesh.pointEdges()));
	labelList 	cellNNeutral(mesh.cells().size(),0);
	labelList 	cellNBoundary(mesh.cells().size(),0);
	labelList 	edgeNeutC(edges.size(),-1);
	labelList 	edgeNiuNei(edges.size(),-1);
	labelList 	pointNiuFac0(points.size(),-1);
	labelList 	pointNiuFac1(points.size(),-1);
	labelList 	pointNiuFac2(points.size(),-1);
	labelList 	BedgesNInternalFaces(edges.size(),0);
	std::vector<int> 	pointNewPoints(points.size(),-2);
	std::vector<int> 	pointNNewFaces(points.size(),0);

	Info<<"  finding Internal edge faces" <<Foam::endl;

	forAll(edgesFs, ei)
	{	const labelList& efs = edgesFs[ei];
		int nedgeNBoundaryFaces = 0;
		forAll(efs, fI)
			if (efs[fI]<nInFaces)
				++BedgesNInternalFaces[ei];
			else
				++nedgeNBoundaryFaces;
		if(nedgeNBoundaryFaces==0) BedgesNInternalFaces[ei]=0;
	}


	int iNewP=-1;
	int iNewF=-1;
	forAll(BedgesNInternalFaces, ei)	if (BedgesNInternalFaces[ei]>1)
	{
		++iNewF;
		edge ej=edges[ei];
		if (pointNewPoints[ej.end()  ]<0)
		{	pointNewPoints[ej.end()  ] = points.size()+(++iNewP);  ++pointNNewFaces[ej.end()  ];}
		else                                                      ++pointNNewFaces[ej.end()  ];
		if (pointNewPoints[ej.start()]<0)
		{	pointNewPoints[ej.start()] = points.size()+(++iNewP);  ++pointNNewFaces[ej.start()  ];}
		else                                                      ++pointNNewFaces[ej.start()  ];
	}///\\\///\\\///\\\///\\\///

	Info<<points.size()<<"  +  "<<iNewP+1<<" new points"<<Foam::endl;
	Info<<faces.size()<<"  +  "<<iNewF+1<<" new faces"<<Foam::endl;

	std::vector<face>  newFaces(iNewF+1, face(4));
	std::vector<int>  newOwners(iNewF+1, -2);
	std::vector<int>  newNeibrs(iNewF+1, -2);
	//std::vector<int>  dblIEdgeEnd(points.size(), 0); /// begging of double-internalFace edges
	//std::vector<int>  dblIEdgeBgn(points.size(), 0); /// begging of double-internalFace edges
	iNewF=-1;
	forAll(BedgesNInternalFaces, ei)	if (BedgesNInternalFaces[ei]>1)
	{
		++iNewF;
		label ownnei[2], iownnei(-1);
		const labelList efs = edgesFs[ei];
		forAll(efs, efi)
		{	label fI=efs[efi];
			if (fI>=nInFaces)
			{
				ownnei[++iownnei]=owns[fI]; ++cellNBoundary[owns[fI]];
			}
			//else {++cellNNeutral[owns[fI]] ; ++cellNNeutral[neis[fI]] ; };
		}
		const labelList ecs = edgesCs[ei];
		Foam::vector CCmid = Foam::vector::zero;
		Foam::scalar sumWeightCCMid = 1e-15;
		forAll(ecs, eci)
		{	label cI=ecs[eci];
			double weight=0.5;
			if (cI!=ownnei[0] && cI!=ownnei[1])
			{
				edgeNeutC[ei]=cI;
				++cellNNeutral[cI];
				weight=2.;

			}
			//else
			{
				CCmid+=weight*CCntres[cI];
				sumWeightCCMid+=weight;
			}
		}
		CCmid/=sumWeightCCMid;

		edge eg=edges[ei];


		if (pointNiuFac0[eg[0]]<0) pointNiuFac0[eg[0]]=iNewF;
		else if (pointNiuFac1[eg[0]]<0) pointNiuFac1[eg[0]]=iNewF;
		else {Info<<" 3Niu "; pointNiuFac2[eg[0]]=iNewF; }

		if (pointNiuFac0[eg[1]]<0) pointNiuFac0[eg[1]]=iNewF;
		else if (pointNiuFac1[eg[1]]<0) pointNiuFac1[eg[1]]=iNewF;
		else {Info<<" 3Niu "; pointNiuFac2[eg[1]]=iNewF; }

		if(ownnei[0]>ownnei[1]) {label tmp=ownnei[0]; ownnei[0]=ownnei[1]; ownnei[1]=tmp; }
		newOwners[iNewF]=ownnei[0];
		newNeibrs[iNewF]=ownnei[1];
		//if(ownnei[0]<ownnei[1])
		Foam::vector CC = CCntres[ownnei[1]]-CCntres[ownnei[0]];
		Foam::vector PP = points[eg.end()]-points[eg.start()];
		Foam::vector CPmid = CCmid - 0.5*(points[eg.end()]+points[eg.start()]);
		if((CC&(PP^CPmid))>=0.)
		{
			newFaces[iNewF][0]=eg.start();
			newFaces[iNewF][1]=eg.end();
			newFaces[iNewF][2]=pointNewPoints[eg.end()];
			newFaces[iNewF][3]=pointNewPoints[eg.start()];



		}
		else
		{
			newFaces[iNewF][1]=eg.start();
			newFaces[iNewF][0]=eg.end();
			newFaces[iNewF][3]=pointNewPoints[eg.end()];
			newFaces[iNewF][2]=pointNewPoints[eg.start()];

			//newOwners[iNewF]=ownnei[1];
			//newNeibrs[iNewF]=ownnei[0];

		}


	}


Info <<" max(cellNNeutral)  "<< max(cellNNeutral)<<Foam::endl;
   ////Returns label of edge nEdges away from startEdge (in the direction of startVertI)
  //Foam::label Foam::meshTools::walkFace (const primitiveMesh& mesh, label facei, label startEdgeI, label startVertI,  label nEdges )
  // {   const labelList& fEdges = mesh.faceEdges(facei);		label edgeI = startEdgeI,  vertI = startVertI;
		//for (label iter = 0; iter < nEdges; iter++)  {  edgeI = otherEdge(mesh, fEdges, edgeI, vertI);   vertI = mesh.edges()[edgeI].otherVertex(vertI);  }
		//return edgeI; }

	forAll(faces, fI)
	{
		const face& ff = faces[fI];
		face ffnew = faces[fI];

		//face sharedps=ff*0-1;
		label ipn(0);
		for(int pi=0; pi<ff.size();++pi,++ipn)  if( pointNewPoints[ff[pi]]>=0 )
		{
			label pI=ff[pi];
			label pIPrv=ff[ff.rcIndex(pi)];
			label pINxt=ff[ff.fcIndex(pi)];

			label niuF0 = pointNiuFac0[pI];
			label niuF1 = pointNiuFac1[pI];
			label niuF2 = pointNiuFac2[pI];
			label ownr = owns[fI];
			label nibr = fI<nInFaces ? neis[fI] : -7;

			if ( !( ( newNeibrs[niuF1]==nibr ) || newNeibrs[niuF1]==ownr || newOwners[niuF1]==ownr  ) )
			 if ( niuF2>=0 && ( ( newNeibrs[niuF2]==nibr ) || newNeibrs[niuF2]==ownr || newOwners[niuF2]==ownr  ) )
				{	label tmp = niuF1; niuF1 = niuF2; niuF2 = tmp;  }
			if ( !( ( newNeibrs[niuF0]==nibr ) || newNeibrs[niuF0]==ownr || newOwners[niuF0]==ownr  ) )
			 if ( niuF1>=0 && ( ( newNeibrs[niuF1]==nibr ) || newNeibrs[niuF1]==ownr || newOwners[niuF1]==ownr  ) )
				{	label tmp = niuF0; niuF0 = niuF1; niuF1 = tmp;  }
			if ( !( ( newNeibrs[niuF1]==nibr ) || newNeibrs[niuF1]==ownr || newOwners[niuF1]==ownr  ) )
			 if ( niuF2>=0 && ( ( newNeibrs[niuF2]==nibr ) || newNeibrs[niuF2]==ownr || newOwners[niuF2]==ownr  ) )
				{	label tmp = niuF1; niuF1 = niuF2; niuF2 = tmp;  }

			//if   (pointNewPoints[pI]>=0)
				//sharedps[++iownnei] = ff[pi];
		//}
		//for(int nn=0;nn<=iownnei;++nn)
		//{
			if( pointNewPoints[pIPrv]>=0 && BedgesNInternalFaces[fEdges[fI][ff.rcIndex(pi)]]>1 && BedgesNInternalFaces[fEdges[fI][pi            ]]<=1)
			{
				//if (1) cout<<1;
				//if ( fI>=nInFaces && ( pointNNewFaces[pI]==3 ) && ( (pointNNewFaces[pINxt]==0)) )
				//{
					 //ffnew[ipn] = pointNewPoints[pI];
				//}
			}
			else
			if( BedgesNInternalFaces[fEdges[fI][pi]]>1 ) //fcIndex
			{
				if(fI<nInFaces)
				{
						//ffnew[ipn] = pointNewPoints[pI];
						//ffnew[ffnew.fcIndex(ipn)] = pointNewPoints[ff.nextLabel(pi)];
						//++ipn;
						//++pi;
					//if(edges[fEdges[fI][pi]][0]==pI)                  ///. the same just to check for now
					//else if(edges[fEdges[fI][pi]][1]==pI) ///. the same just to check for now
					//else Info<<"badeege"<<Foam::endl;

					//if (newNeibrs[niuF0]>0)
					//{  face ftmp2=ffnew;
						//forAll(ffnew,pi)   ffnew[ffnew.size()-pi-1]=ftmp2[pi];
						//label tmp = ownr;	ownr = nibr; nibr = tmp;
					//}
					label niuF0Nxt = pointNiuFac0[pINxt];
					label niuF1Nxt = pointNiuFac1[pINxt];
					label niuF2Nxt = pointNiuFac1[pINxt];
					if ( !( ( newNeibrs[niuF1Nxt]==nibr ) || newNeibrs[niuF1Nxt]==ownr || newOwners[niuF1Nxt]==ownr  ) )
					 if ( niuF2Nxt>=0 && ( ( newNeibrs[niuF2Nxt]==nibr ) || newNeibrs[niuF2Nxt]==ownr || newOwners[niuF2Nxt]==ownr  ) )
						{	label tmp = niuF1Nxt; niuF1Nxt = niuF2Nxt; niuF2Nxt = tmp;  }
					if ( !( ( newNeibrs[niuF0Nxt]==nibr ) || newNeibrs[niuF0Nxt]==ownr || newOwners[niuF0Nxt]==ownr  ) )
					 if ( niuF1Nxt>=0 && ( ( newNeibrs[niuF1Nxt]==nibr ) || newNeibrs[niuF1Nxt]==ownr || newOwners[niuF1Nxt]==ownr  ) )
						{	label tmp = niuF0Nxt; niuF0Nxt = niuF1Nxt; niuF1Nxt = tmp;  }
					if ( !( ( newNeibrs[niuF1Nxt]==nibr ) || newNeibrs[niuF1Nxt]==ownr || newOwners[niuF1Nxt]==ownr  ) )
					 if ( niuF2Nxt>=0 && ( ( newNeibrs[niuF2Nxt]==nibr ) || newNeibrs[niuF2Nxt]==ownr || newOwners[niuF2Nxt]==ownr  ) )
						{	label tmp = niuF1Nxt; niuF1Nxt = niuF2Nxt; niuF2Nxt = tmp;  }

					if (   ( (BedgesNInternalFaces[fEdges[fI][ff.rcIndex(pi)]]<=1 && BedgesNInternalFaces[fEdges[fI][ff.fcIndex(pi)]]<=1) /*&& ((BedgesNInternalFaces[fEdges[fI][ff.rcIndex(pi)]]==1 || BedgesNInternalFaces[fEdges[fI][ff.fcIndex(pi)]]==1)  )*/)
						 //&& ( (pointNNewFaces[pI]!=2) || (pointNNewFaces[pINxt]!=2) )
						 //&& ( ( pointNNewFaces[pI]>1 && ( newOwners[niuF0]==ownr || newNeibrs[niuF0]==nibr /*|| newOwners[niuF0]==nibr || newNeibrs[niuF0]==ownr*/ )   //)
														  //&& ( /*newOwners[niuF1]==ownr || newNeibrs[niuF1]==nibr ||*/ newOwners[niuF1]==nibr || newNeibrs[niuF1]==ownr
														  //|| (niuF2>0 && ( /*newOwners[niuF2]==ownr || newNeibrs[niuF2]==nibr ||*/ newOwners[niuF2]==nibr || newNeibrs[niuF2]==ownr ) )
														  //)     //|| newOwners[niuF1]==nibr || newNeibrs[niuF1]==ownr )
							//)
						 //|| ( pointNNewFaces[pINxt]>1 && ( newOwners[niuF0Nxt]==ownr || newNeibrs[niuF0Nxt]==nibr /*|| newOwners[niuF0Nxt]==nibr || newNeibrs[niuF0Nxt]==ownr*/ )   //)
														     //&& ( newOwners[niuF1Nxt]==ownr || newNeibrs[niuF1Nxt]==nibr /*|| newOwners[niuF1Nxt]==nibr || newNeibrs[niuF1Nxt]==ownr  */
														     //|| (niuF2Nxt>0 && ( newOwners[niuF2Nxt]==ownr || newNeibrs[niuF2Nxt]==nibr /*|| newOwners[niuF2Nxt]==nibr || newNeibrs[niuF2Nxt]==ownr*/ ) )
														     //)  //|| newOwners[niuF1Nxt]==nibr || newNeibrs[niuF1Nxt]==ownr)
							//)
							 //)
						 &&( ( pointNNewFaces[pI]>1
							&& ( newOwners[niuF0]==ownr || newNeibrs[niuF0]==ownr || newOwners[niuF1]==ownr || newNeibrs[niuF1]==ownr || (niuF2>0 && ( newOwners[niuF2]==ownr || newNeibrs[niuF2]==ownr ) ) )
							&& ( newOwners[niuF0]==nibr || newNeibrs[niuF0]==nibr || newOwners[niuF1]==nibr || newNeibrs[niuF1]==nibr || (niuF2>0 && ( newOwners[niuF2]==nibr || newNeibrs[niuF2]==nibr ) ) )
							  )
						 ||  ( pointNNewFaces[pINxt]>1
							&& ( newOwners[niuF0Nxt]==ownr || newNeibrs[niuF0Nxt]==ownr || newOwners[niuF1Nxt]==ownr || newNeibrs[niuF1Nxt]==ownr || (niuF2Nxt>0 && ( newOwners[niuF2Nxt]==ownr || newNeibrs[niuF2Nxt]==ownr ) ) )
							&& ( newOwners[niuF0Nxt]==nibr || newNeibrs[niuF0Nxt]==nibr || newOwners[niuF1Nxt]==nibr || newNeibrs[niuF1Nxt]==nibr || (niuF2Nxt>0 && ( newOwners[niuF2Nxt]==nibr || newNeibrs[niuF2Nxt]==nibr ) ) )
							  )
						  )
						)
					{
						face fftmp =ffnew;
						ffnew.append(-1);
						for(int pii=ffnew.size()-1; pii>ipn;--pii)  ffnew[pii] = ffnew[pii-1];

						bool cond =  ( pointNNewFaces[pINxt]>1
										&& ( newOwners[niuF0Nxt]==ownr || newNeibrs[niuF0Nxt]==ownr || newOwners[niuF1Nxt]==ownr || newNeibrs[niuF1Nxt]==ownr || (niuF2Nxt>0 && ( newOwners[niuF2Nxt]==ownr || newNeibrs[niuF2Nxt]==ownr ) ) )
										&& ( newOwners[niuF0Nxt]==nibr || newNeibrs[niuF0Nxt]==nibr || newOwners[niuF1Nxt]==nibr || newNeibrs[niuF1Nxt]==nibr || (niuF2Nxt>0 && ( newOwners[niuF2Nxt]==nibr || newNeibrs[niuF2Nxt]==nibr ) ) )
										 )
										;
						//bool cond =  ( ( newOwners[niuF0]==ownr || newNeibrs[niuF0]==nibr ) )
										//&& (niuF1<0 || newOwners[niuF1]==ownr || newNeibrs[niuF1]==nibr )
										//&& (niuF2<0 || newOwners[niuF2]==ownr || newNeibrs[niuF2]==nibr )
										//;
						if(cond)
						{
							ffnew[ipn] = pointNewPoints[pI];
							ffnew[++ipn] = pointNewPoints[pINxt];
							ffnew[ffnew.fcIndex(ipn)] = pINxt;
						}
						else
						{
							ffnew[ipn++] = pI;
							ffnew[ipn] = pointNewPoints[pI];
							ffnew[ffnew.fcIndex(ipn)] = pointNewPoints[pINxt];
						}
					}
					else
					{
						ffnew[ipn] = pointNewPoints[pI];
						ffnew[ffnew.fcIndex(ipn)] = pointNewPoints[ff.nextLabel(pi)];
					}
					//if (newOwners[niuF0]<0)
					//{  face ftmp2=ffnew;
						//forAll(ffnew,pi)   ffnew[ffnew.size()-pi-1]=ftmp2[pi];
						//label tmp = owns[fI];	owns[fI] = neis[fI]; neis[fI] = tmp;
					//}
				}
				//else if ( ( pointNNewFaces[pI]==3 ) && ( (pointNNewFaces[pINxt]==0)) )
				//{
						//ffnew[ffnew.fcIndex(ipn)] = pointNewPoints[ff.nextLabel(pi)];
				//}
					//++ipn; ++pi;
			}
			else if( ( fI<nInFaces &&  ( newOwners[niuF0]==nibr || newNeibrs[niuF0]==nibr
													|| (niuF1>0 && ( newOwners[niuF1]==nibr || newNeibrs[niuF1]==nibr ) )
													|| (niuF2>0 && ( newOwners[niuF2]==nibr || newNeibrs[niuF2]==nibr ) )
												)
						)
						||
						( (fI<nInFaces || pointNNewFaces[pI]==1) &&
														 (
															               (  ( newOwners[niuF0]==ownr || newNeibrs[niuF0]==ownr ) )
															||  (niuF1>0 && ( newOwners[niuF1]==ownr || newNeibrs[niuF1]==ownr ) )
															||  (niuF2>0 && ( newOwners[niuF2]==ownr || newNeibrs[niuF2]==ownr ) )
														 )
						)

					) //pointNNewFaces[pI]==1 &&
			{
				face fftmp =ffnew;
				ffnew.append(-1);
				for(int pii=ffnew.size()-1; pii>ipn;--pii)  ffnew[pii] = ffnew[pii-1];

					//bool cond = ( fI<nInFaces &&            (newOwners[niuF0]==nibr  || newNeibrs[niuF0]==nibr || (niuF1<0 && newNeibrs[niuF1]==nibr) ) )
				           //||  ( pointNNewFaces[pI]==1 && ( (newOwners[niuF0]==owns[fI] && dblIEdgeEnd[pI]) || (newNeibrs[niuF0]==owns[fI] && !dblIEdgeEnd[pI]) ) );

					bool cond =  ( fI<nInFaces && BedgesNInternalFaces[fEdges[fI][pi]]==0 )
									//( fI<nInFaces && (( (newOwners[niuF0]==nibr && dblIEdgeEnd[pI])  || (newNeibrs[niuF0]==nibr&& !dblIEdgeEnd[pI] ) )
									//|| ( niuF1>0 && ( (newOwners[niuF1]==nibr && !dblIEdgeEnd[pI]) || (newNeibrs[niuF1]==nibr&& dblIEdgeEnd[pI] ))  ) ))
						        || ( fI>=nInFaces && BedgesNInternalFaces[fEdges[fI][pi]]==1 ) ; //( (newOwners[niuF0]==owns[fI] && dblIEdgeEnd[pI]) || (newNeibrs[niuF0]==owns[fI] && !dblIEdgeEnd[pI]) ) );

					if(cond)	ffnew[++ipn] = pointNewPoints[pI];
					else		ffnew[ipn++] = pointNewPoints[pI];

			}
			else if ( (fI<nInFaces || pointNNewFaces[pI]==1 ) )
			{
					 ffnew[ipn] = pointNewPoints[pI];
			}

		}




		//if (fI<nInFaces && owns[fI]>nibr)
		//{  face ftmp2=ffnew;
			//forAll(ffnew,pi)   ffnew[ffnew.size()-pi-1]=ftmp2[pi];
			//label tmp = owns[fI];	owns[fI] = nibr; nibr = tmp;
		//}


		faces[fI].resize(ffnew.size()); faces[fI]=ffnew;
	}///\\\///\\\///\\\///\\\///

	Info<<"-----------"<<Foam::endl;


//{

	std::cout <<"writing points";std::cout.flush();
	std::ofstream pointFile("constant/polyMesh/points");
	std::vector<int>  newPointsMap(iNewP+1, -1);
	cout<<"nNewPoints: "<<newPointsMap.size()<<std::endl;
	for (unsigned int pi=0;pi<pointNewPoints.size();++pi) 	if (pointNewPoints[pi]>=0)
		newPointsMap[pointNewPoints[pi]-points.size()]=pi;;


	assert(pointFile);
	pointFile<<
	"FoamFile\n{\n	version	 2.0;\n	format	  ascii;\n"
	"	class	   vectorField;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  points;\n"
	"}\n\n";
	pointFile<<newPointsMap.size()+points.size()<<"\n("<<std::endl;

	forAll(points, pi) pointFile<< "("<<points[pi][0]<< ' '<<points[pi][1]<<' '<<points[pi][2]<<")\n";
	//for(const auto pi :  newPointsMap) { pointFile<< "("<<0<< ' '<<0<<' '<<0<<")\n"; }
	for(const auto pi :  newPointsMap) { pointFile<< "("<<points[pi][0]<< ' '<<points[pi][1]<<' '<<points[pi][2]<<")\n"; }

	pointFile<<"\n)"<<std::endl;

	pointFile.close();

//}

	std::cout <<"cxcxc";std::cout.flush();


{
	std::ofstream boundary("constant/polyMesh/boundary");	assert(boundary);
	boundary<<
	"FoamFile\n{\n	version	 2.0;\n	format	  ascii;\n"
	"	class	   polyBoundaryMesh;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  boundary;\n"
	"}\n\n";


	boundary<<mesh.boundary().size()  <<"\n("<<std::endl;
	forAll(mesh.boundary(), bpi)
	{
		const Foam::polyPatch&  bi = mesh.boundary()[bpi].patch();
		boundary<<
		"	"<<  bi.name()							<<std::endl<<
		"	{"													<<std::endl<<
		"		type			"<<"patch"/*bi.physicalType()*/<<";\n"	<<
		"		nFaces		  "<<bi.size()<<";\n"	  	<<
		"		startFace	   "<<bi.start()+newFaces.size()<<";\n"<<
		"	}"												 <<std::endl;
	}
	boundary<<")"   <<std::endl;
	boundary.close();

}


{


//______________________________________________________________________________________________________________________

	std::ofstream facfile("constant/polyMesh/faces");	assert(facfile);
	facfile<<
	"FoamFile\n{\n	version	 2.0;\n	format	  ascii;\n"
	"	class	   faceList;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  faces;\n"
	"}\n\n";


	std::ofstream owner("constant/polyMesh/owner");	assert(owner);
	owner<<
	"FoamFile\n{\n	version	 2.0;\n	format	  ascii;\n"
	"	class	   labelList;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  owner;\n"
	"}\n\n";

	std::ofstream neighbour("constant/polyMesh/neighbour");	assert(neighbour);
	neighbour<<
	"FoamFile\n{\n	version	 2.0;\n	format	  ascii;\n"
	"	class	   labelList;\n"
	"	location	\"constant/polyMesh\";\n"
	"	object	  neighbour;\n"
	"}\n\n";



//______________________________________________________________________________________________________________________


	facfile<<faces.size()+newOwners.size()<<"\n("<<std::endl;
	for (int fI=0;fI<nInFaces;++fI)
	{  const face& ff = faces[fI];
	   facfile<<'(';  forAll(ff,i) { facfile<<ff[i]<<" "; }  facfile<<")\n";  }
	for (unsigned int fI=0;fI<newFaces.size();++fI)
	{  const face& ff = newFaces[fI];
	   facfile<<'(';  forAll(ff,i) { facfile<<ff[i]<<" "; }  facfile<<")\n";  }
	for (int fI=nInFaces;fI<faces.size();++fI)
	{  const face& ff = faces[fI];
	   facfile<<'(';  forAll(ff,i) { facfile<<ff[i]<<" "; }  facfile<<")\n";  }
	facfile<<")"<<std::endl;
	facfile.close();

	owner<<faces.size()+newOwners.size()<<"\n("<<std::endl;
	for (int fI=0;fI<nInFaces;++fI)
		owner<<owns[fI]<<"\n";
	for (unsigned int fI=0;fI<newOwners.size();++fI)
		owner<<newOwners[fI]<<"\n";
	for (int fI=nInFaces;fI<faces.size();++fI)
		owner<<owns[fI]<<"\n";
	owner<<")"<<std::endl;
	owner.close();


	neighbour<<nInFaces+newNeibrs.size()<<"\n("<<std::endl;
	for (int fI=0;fI<nInFaces;++fI)
		neighbour<<neis[fI]<<"\n";
	for (unsigned int fI=0;fI<newNeibrs.size();++fI)
		neighbour<<newNeibrs[fI]<<"\n";
	neighbour<<")"<<std::endl;
	neighbour.close();

}


    Info<< "\nEnd.\n" << Foam::endl;

    return 0;

}
// ************************************************************************* //
