/*-------------------------------------------------------------------------*\

You can redistribute this code and/or modify this code under the 
terms of the GNU General Public License (GPL) as published by the  
Free Software Foundation, either version 3 of the License, or (at 
your option) any later version. see <http://www.gnu.org/licenses/>.

Please see our website for relavant literature making use of this code:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini: a.q.raeini@imperial.ac.uk

\*-------------------------------------------------------------------------*/


#include <fstream>
#include <iostream>
#include <vector>

#include <assert.h>
#include "voxelImage.h"

using namespace std;

int usage()
{
	std::cout<<"Ufraw2Uc  "<<std::endl;
		std::cout
		<<" Converts face centred velocities (Uf*s) to cell centred velocities (Uc*s)\n"
		<<" Writes velocity magnitide if third argument (optional) is Umag or UmagOnly\n\n"
		<<"Usage examples, type:\n"
		<<" #cd PATH/TO/Ufx.*etc first"<< std::endl
		<<" Ufraw2Uc raw  vxlImage.mhd          # generate Uc*s"<< std::endl
		<<" Ufraw2Uc raw  vxlImage.mhd Umag     # write mag(U) as well"<< std::endl
		<<" Ufraw2Uc raw  vxlImage.mhd UmagOnly # write mag(U) only"<< std::endl
		<<" Ufraw2Uc dat  vxlImage.mhd          # write Uc*s in ascii format"<< std::endl;
	return 1;
}

int main(int argc, char** argv)
{

	if(argc<3)		return usage();
	std::string outFormat(argv[1]);
	imgExt(outFormat);
	std::string headerName(argv[2]);
	if(headerName.size()<4 || headerName.compare(0,headerName.size(),"vxlImage.mhd") != 0) return usage();
	
	string wUmag;
	if(argc>3)	
	{
		wUmag=string(argv[1]);
		if(wUmag[0]=='U') std::cout<<"Writing Umag"<<wUmag<<endl;
		else	 std::cout<<"Warning third argument can only be wUmag or UmagOnly, ignoring it:"<<wUmag<<endl;
	}else wUmag="ignor";

	voxelImage vimage("vxlImage.mhd");
	if(!vimage.nz()) {std::cout<<"Error: vxlImage.mhd not read"<<std::endl; return 1;}
	int3 n=vimage.size3();
	if(outFormat=="dat") vimage.write("vxlImage.dat");
	vimage.data_.resize(0);

	voxelImageT<float> Umg;
	if (wUmag[0]=='U') Umg.reset(n[0],n[1],n[2],0.0);

	{
		voxelImageT<float> fField(n[0]+1,n[1],n[2],0.0);
		fField.readBin("Ufx"+imgExt());

		for (int k = 0; k<int(fField.nz()) ; k++ )
		 for ( int j = 0; j<int(fField.ny()) ; j++ )
		  for ( int i = 0; i<int(fField.nx())-1 ; i++ )
				fField(i,j,k)=0.5*(fField(i,j,k)+fField(i+1,j,k));
		if(wUmag!="UmagOnly")
			fField.writeBin("Uccx"+imgExt(), 0,n[0],0,n[1],0,n[2]);
		if(wUmag[0]=='U')
			forAllkji(Umg) Umg(i,j,k) += fField(i,j,k)*fField(i,j,k);
	}
	{
		voxelImageT<float> fField(n[0],n[1]+1,n[2],0.0);
		fField.readBin("Ufy"+imgExt());

		for (int k = 0; k<int(fField.nz()) ; k++ )
		 for ( int j = 0; j<int(fField.ny())-1 ; j++ )
		  for ( int i = 0; i<int(fField.nx()) ; i++ )
			fField(i,j,k)=0.5*(fField(i,j,k)+fField(i,j+1,k));
		if(wUmag!="UmagOnly")
			fField.writeBin("Uccy"+imgExt(), 0,n[0],0,n[1],0,n[2]);
		if(wUmag[0]=='U')
			forAllkji(Umg) Umg(i,j,k) += fField(i,j,k)*fField(i,j,k);
	}
	{
		voxelImageT<float> fField(n[0],n[1],n[2]+1,0.0);
		fField.readBin("Ufz"+imgExt());
		for (int k = 0; k<int(fField.nz())-1 ; k++ )
		 for ( int j = 0; j<int(fField.ny()) ; j++ )
		  for ( int i = 0; i<int(fField.nx()) ; i++ )
			fField(i,j,k)=0.5*(fField(i,j,k)+fField(i,j,k+1));
		if(wUmag!="UmagOnly")
			fField.writeBin("Uccz"+imgExt(), 0,n[0],0,n[1],0,n[2]);
		if(wUmag[0]=='U')
			forAllkji(Umg) Umg(i,j,k) += fField(i,j,k)*fField(i,j,k);
	}
	if(wUmag[0]=='U')
	{
		forAlliii_(Umg) Umg(iii) = sqrt(Umg(iii));
		Umg.writeBin("Umag"+imgExt());
	}

	if(outFormat=="dat") 
	{
		//voxelImageT<float> pField(n[0],n[1],n[2],0.0);
		
		//else pField.readBin("p.tif");
		//pField.writeAscii("p.dat");
	}

	std::cout<< "end" << std::endl;


	return 0;
}


