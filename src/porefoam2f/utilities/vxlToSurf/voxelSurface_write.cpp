/*-------------------------------------------------------------------------*\
 This is part of surfLib, a library for working with surface files and data

 Copyright (C) 2018-2020  Ali Qaseminejad Raeini 

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

using namespace std;



#include "surfUtils.h"


//#include "yocto_utils.h"
dbl3s mapFacesGetPoints(std::vector<std::vector<face> >& facess, const piece<point> pointsAll)  {
	ints pValidIs(pointsAll.size(),-1); int pInd=-1;
	for(auto&faces:facess) for(auto&fac:faces) for(auto& pi:fac) if(pValidIs[pi]==-1) pValidIs[pi]=++pInd;
	dbl3s points(pInd+1);
	for(auto&faces:facess) for(auto&fac:faces) for(auto& pi:fac) pi=pValidIs[pi];
	for_i(pValidIs) if(pValidIs[i]!=-1) { points[pValidIs[i]]=pointsAll[i]; }
	return std::move(points);
}

std::vector<std::vector<face> > getbsoleteOrderedFaces(const facePieceList& facezsz, int vv)  {
	std::vector<std::vector<face> > facess(256);
	{
		ints nFacs(256,0);
		const facePiece& facezs=facezsz[vv];
		for(const auto&fac:facezs) { ensure(0<=fac.zone && fac.zone<255);  ++nFacs[fac.zone];}
		for(int i=0;i<256;++i) if(nFacs[i]) { facess[i].resize(nFacs[i], face({-1,-1,-1,-1})); nFacs[i]=-1; }
		for(const auto&fac:facezs) { facess[fac.zone][++nFacs[fac.zone]]=fac; ensure(facess[fac.zone][nFacs[fac.zone]][0]>=0);}
	}


	auto invertFace = [](const face& f) { return face({f[3],f[2],f[1],f[0]}); };
	for(int v2=0;v2<vv;++v2) if(facezsz[v2].size()) // when searching for faces only faces toward higher labels are extracted, here we also add those toward lower voxel-labels
	{	/// beak each voxel value surfacce manifold into zones
		ints nFacs(256,0);
		const facePiece& facezs=facezsz[v2];
		for(const auto&fac:facezs) { if(fac.zone==vv)  ++nFacs[v2];}
		for(int iz=0;iz<vv;++iz) if(nFacs[iz]) { facess[iz].resize(nFacs[iz], face({-1,-1,-1,-1})); nFacs[iz]=-1; }
		for(const auto&fac:facezs)  if(fac.zone==vv) { facess[v2][++nFacs[v2]]=invertFace(fac);  ensure(facess[v2][nFacs[v2]][0]>=0);} // wrong zone
	}
	return std::move(facess);
}

void writeSurfaceFiles(const facePieceList& facezsz, const piece<point> pointsAll, const std::string& fileNames)  {

 for(size_t vv=0;vv<facezsz.size();++vv) if(facezsz[vv].size())  {
	string ext = fileNames.substr(fileNames.size()-4);
	string basNam = fileNames.substr(0,fileNames.size()-4)+_s(vv);

	cout<<"writing surface "<<basNam+ext<<" for voxels values "<<vv<<endl;
	/// beak each voxel value surface manifold into zones

	std::vector<std::vector<face> > facess = getbsoleteOrderedFaces(facezsz, vv);

	dbl3s points = mapFacesGetPoints(facess,pointsAll);
	/// get rid of unused points

	int nFaces=0, nZones=0;
	for(const auto& rg:facess) { nFaces+=rg.size(); nZones+=rg.size()>0; }

	if(ext==".obj")  {
		ofstream fil(basNam+ext);
		ensure(fil.good());

		fil	<< "o " << basNam<<"\n"
			<< "# n_points : " << points.size() << "\n"
			<< "# n_faces  : " << nFaces << "\n"
			<< "# n_zones  : " << nZones << "\n";

		// Print zone names as comment
		for_(facess, zoneI)
		  if(facess[zoneI].size())
			fil  << "#	" << zoneI << "  zone" << zoneI
				 << "  (nFaces: " << facess[zoneI].size() << ")" << "\n";

		// Write vertex coords
		for_(points, ptI)  {
			const point& pt = points[ptI];
			fil  << "v " << pt.x << ' '  << pt.y << ' '  << pt.z << "\n";
		}

		fil  << endl;

		for_(facess, zoneI)  {
			const facePiece& zone = facess[zoneI];
			fil << "vt  " << zoneI/255.0 <<' '<< zoneI/255.0 << endl;

			if (zone.size())  {
				fil << "g  zone" << zoneI << endl;
				fil << "#  zone"<<zoneI<<" size: " << zone.size() << endl;
				{
					for_(zone, zfI)  {
						const face& f = zone[zfI];

						fil << 'f';
						for_(f, fp)  fil << ' ' << f[fp] + 1<<'/'<<zoneI;//+1
						fil << "\n";
					}
				}
			}
		}
		fil << endl;
	}
	else 
	{	if(ext!=".vtk") 		cerr << "outputSurface format "<<fileNames<<"not supported switching to vtk";

		ofstream fil(basNam+".vtk");
		ensure(fil.good());

		fil << "# vtk DataFile Version 2.0\n"
			<< "written by voxelImage library by Ali Q. Raeini ...\n"
			<< "ASCII\n\n"
			<< "DATASET POLYDATA"<<endl;

		// Write vertex coords
		fil  << "POINTS " << points.size() << " float" << '\n';
		for_(points, ptI)  {
			const point& pt = points[ptI];

			fil  << pt.x << ' ' << pt.y << ' ' << pt.z << '\n';
		}

		int nFaces=0, nZones=0;
		for(const auto& rg:facess) { nFaces+=rg.size(); nZones+=rg.size()>0; }

		label nNodes = nFaces*4;
		fil << '\n'
			<< "POLYGONS " << nFaces << ' ' << nFaces + nNodes << '\n';


		for_(facess, zoneI)  {
			const facePiece& zone = facess[zoneI];

			for_(zone, zfI)  {
				const face& f = zone[zfI];
				fil << f.size();
				for_(f, fp)  fil << ' ' << f[fp];
				fil << '\n';
			}
		}


		fil  << '\n'
			<< "CELL_DATA " << nFaces << '\n'
			<< "FIELD attributes 1\n"
			<< "zone 1 " << nFaces << " int\n" ;

		for_(facess, zoneI)  {
			const facePiece& zone = facess[zoneI];

			for_(zone, zfI)  {
				fil << zoneI;//+1
				fil << '\n';
			}
		}
	}


 }
}



void writeMergeSurfaceFile(const facePieceList& facezsz, const piece<point> points, const std::string& fileNames) //,std::vector<string> rockTypes
{
	string ext = fileNames.substr(fileNames.size()-4);
	string basNam = fileNames.substr(0,fileNames.size()-4);
	(cout<<"writing surface "<<basNam+ext<<" ").flush();

	std::vector<std::vector<face> > facess(256*256);

	for(size_t vv=0;vv<facezsz.size();++vv) if(facezsz[vv].size())  {
		(cout<<" "<<vv).flush();
		/// beak  each voxel value surfacce manifold into zones
		const facePiece& facezs=facezsz[vv];
		ints nFacs(256,0);
		for(const auto&fac:facezs) { ensure(0<=fac.zone && fac.zone<255);  ++nFacs[fac.zone];}
		for(int iz=0;iz<256;++iz) if(nFacs[iz]) { facess[iz+vv*256].resize(nFacs[iz], face({-1,-1,-1,-1})); nFacs[iz]=-1; }
		for(const auto&fac:facezs) { facess[fac.zone+vv*256][++nFacs[fac.zone]]=fac;  ensure(facess[fac.zone+vv*256][nFacs[fac.zone]][0]>=0);}
		(cout<<".").flush();
	}
	(cout<<" .").flush();

	int nFaces=0, nZones=0;
	for(const auto& rg:facess) { nFaces+=rg.size(); nZones+=rg.size()>0; }

	if(ext==".obj")  {
		ofstream fil(basNam+ext);
		ensure(fil.good());

		fil	<< "o " << basNam<<"\n"
			<< "# n_points : " << points.size() << "\n"
			<< "# n_faces  : " << nFaces << "\n"
			<< "# n_zones  : " << nZones << "\n";

		// Print zone names as comment
		for_(facess, zoneI)
		  if(facess[zoneI].size())
			fil  << "#	" << zoneI << "  v"<< zoneI/256<<"_to_v"<<zoneI%256
				 << "  (nFaces: " << facess[zoneI].size() << ")" << "\n";

		// Write vertex coords
		for_(points, ptI)  {
			const point& pt = points[ptI];
			fil  << "v " << pt.x << ' '  << pt.y << ' '  << pt.z << "\n";
		}

		fil  << endl;

		for_(facess, zoneI)  {
			const facePiece& zone = facess[zoneI];
			fil << "vt  " << zoneI/255.0 <<' '<< zoneI/255.0 << endl;

			if (zone.size())  {
				fil << "g  v" << zoneI/256<<"_to_v"<<zoneI%256<< endl;
				fil << "#  v" << zoneI/256<<"_to_v"<<zoneI%256<<" size: " << zone.size() << endl;
				{
					for_(zone, zfI)  {
						const face& f = zone[zfI];
						fil << 'f';
						for_(f, fp)  fil << ' ' << f[fp] + 1<<'/'<<zoneI;//+1
						fil << "\n";
					}
				}
			}
		}
		fil << endl;
	}
	else 
	{	if(ext!=".vtk") 		cerr << "outputSurface format "<<ext<<" not supported switching to .vtk";

		ofstream fil(basNam+".vtk");
		ensure(fil.good());

		fil << "# vtk DataFile Version 2.0\n"
			<< "written by voxelImage library by Ali Q. Raeini ...\n"
			<< "ASCII\n\n"
			<< "DATASET POLYDATA"<<endl;

		// Write vertex coords
		fil  << "POINTS " << points.size() << " float" << '\n';
		for_(points, ptI)   fil  << points[ptI] << '\n';


		int nFaces=0;
		for(const auto& rg:facess) { nFaces+=rg.size(); }

		fil << "\n""POLYGONS " << nFaces << ' ' << nFaces + nFaces*4 << '\n';


		for_(facess, zoneI)  {
			const facePiece& zone = facess[zoneI];

			for_(zone, zfI)  {
				const face& f = zone[zfI];
				fil << f.size();
				for_(f, fp)  fil << ' ' << f[fp];
				fil << '\n';
			}
		}


		fil  << '\n'
			<< "CELL_DATA " << nFaces << '\n'
			<< "FIELD attributes 1\n"
			<< "zone 1 " << nFaces << " int\n" ;

		for_(facess, zoneI)  {
			const facePiece& zone = facess[zoneI];

			for_(zone, zfI)   fil << zoneI*256+zone[zfI].zone << '\n';
		}
	}


}



