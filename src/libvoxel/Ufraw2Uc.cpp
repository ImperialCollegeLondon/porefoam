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


#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>
#include "voxelImage.h"

#include "voxelImage.cpp"

using namespace std;

int usage()  {
	std::cout<<"Ufraw2Uc  "<<std::endl;
		std::cout<<
		  " Converts face centred velocities (Uf*s) to cell centred velocities (Uc*s)\n"
		  " Writes velocity magnitide if third argument (optional) is Umag or UmagOnly\n"
		  " Converts also file format if third arg starts with '.' (followed by output format)\n\n"
		  "Usage examples, type:\n"
		  " #cd PATH/TO/Ufx.*etc first\n"
		  " Ufraw2Uc raw  vxlImage.mhd          # generate Uc*s\n"
		  " Ufraw2Uc raw  vxlImage.mhd  Umag     # write mag(U) as well\n"
		  " Ufraw2Uc raw  vxlImage.mhd  UmagOnly # write mag(U) only\n"
		  " Ufraw2Uc dat  vxlImage.mhd          # write Uc*s in ascii format\n"
		  " Ufraw2Uc tif  vxlImage.mhd  .raw.gz    # read .tif and write .raw files"<< std::endl;
	return 1;
}

int main(int argc, char** argv)  {

	if(argc<3)		return usage();
	imgExt(string(argv[1]));
	std::string headerName(argv[2]);
	if(headerName.size()<4 || headerName.compare(0,headerName.size(),"vxlImage.mhd") != 0)  usage();
	
	string wUmag;
	if(argc>3)	
	{
		wUmag=string(argv[3]);
		if(wUmag[0]=='U') std::cout<<"Writing Umag: "<<wUmag<<endl;
		else ensure(wUmag[0]!='.',"third argument can only be wUmag or UmagOnly or start with '.' (.raw...), ignoring it: "+wUmag);
	}else wUmag="ignor";

	string writeFrmt = (wUmag[0]=='.') ? wUmag :  imgExt();

	std::cout<<" Ufraw2Uc "+imgExt()+"  "+headerName+"  "<<writeFrmt<<endl;

	voxelImage vimage(headerName);
	ensure(vimage.nz(), headerName+" not read",-1);
	int3 n=vimage.size3();
	if(writeFrmt=="dat") vimage.write("vxlImage.dat");
	vimage.data_.resize(0);

	voxelImageT<float> Umg;
	if (wUmag[0]=='U') Umg.reset(n[0],n[1],n[2],0.);

	{
		voxelImageT<float> fField(n[0]+1,n[1],n[2],0.);
		fField.readBin("Ufx"+imgExt());
		if(writeFrmt!=imgExt())  fField.writeBin("Ufx"+writeFrmt);

		for (int k=0; k<fField.nz(); ++k)
		 for (int j=0; j<fField.ny(); ++j)
		  for (int i=0; i<fField.nx()-1; ++i)
				fField(i,j,k)=0.5*(fField(i,j,k)+fField(i+1,j,k));
		if(wUmag!="UmagOnly")
			fField.writeBin("Uccx"+writeFrmt, 0,n[0],0,n[1],0,n[2]);
		if(wUmag[0]=='U')
			forAllkji(Umg) Umg(i,j,k) += sqr(fField(i,j,k));
	}
	{
		voxelImageT<float> fField(n[0],n[1]+1,n[2],0.);
		fField.readBin("Ufy"+imgExt());
		if(writeFrmt!=imgExt())  fField.writeBin("Ufy"+writeFrmt);

		for (int k=0; k<fField.nz(); ++k)
		 for (int j=0; j<fField.ny()-1; ++j)
		  for (int i=0; i<fField.nx(); ++i)
			fField(i,j,k)=0.5*(fField(i,j,k)+fField(i,j+1,k));
		if(wUmag!="UmagOnly")
			fField.writeBin("Uccy"+writeFrmt, 0,n[0],0,n[1],0,n[2]);
		if(wUmag[0]=='U')
			forAllkji(Umg) Umg(i,j,k) += sqr(fField(i,j,k));
	}
	{
		voxelImageT<float> fField(n[0],n[1],n[2]+1,0.);
		fField.readBin("Ufz"+imgExt());
		if(writeFrmt!=imgExt())  fField.writeBin("Ufz"+writeFrmt);

		for (int k=0; k<fField.nz()-1; ++k)
		 for (int j=0; j<fField.ny(); ++j)
		  for (int i=0; i<fField.nx(); ++i)
			fField(i,j,k)=0.5*(fField(i,j,k)+fField(i,j,k+1));
		if(wUmag!="UmagOnly")
			fField.writeBin("Uccz"+writeFrmt, 0,n[0],0,n[1],0,n[2]);
		if(wUmag[0]=='U')
			forAllkji(Umg) Umg(i,j,k) += sqr(fField(i,j,k));
	}
	if(wUmag[0]=='U')  {
		forAlliii_(Umg) Umg(iii) = sqrt(Umg(iii));
		Umg.writeBin("Umag"+writeFrmt);
	}

	if(writeFrmt!=imgExt()) 
	{
		voxelImageT<float> pField(n[0],n[1],n[2],0.);
		pField.readBin("p"+imgExt());
		pField.writeBin("p"+writeFrmt);

		pField.readBin("psi"+imgExt());
		pField.writeBin("psi"+writeFrmt);
	}

	std::cout<< "end" << std::endl;


	return 0;
}


