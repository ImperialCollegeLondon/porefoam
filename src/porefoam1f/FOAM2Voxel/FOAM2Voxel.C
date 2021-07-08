/*-------------------------------------------------------------------------*\
 Conversion utility, from Openfoam results to 3D images

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

#include <stdio.h>
#include <stdlib.h>     /* malloc, free, rand */
#include <new>
#include <fstream>
#include <iostream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <assert.h>

#include "voxelImage.h"


#include "FOAMProc2Voxel.H"

//template<typename Type>
//std::string  _s(Type str)
//{
	//std::stringstream ss;
	//ss<<str;
	//return ss.str();
//}

int main(int argc, char *argv[])
{

	std::string outputFormat="binary";
	std::string headerName = (argc>1) ? std::string (argv[1]) : "";
	int nProcs=1;
	char _case[] = "-case\0";
	bool initialiseOF=true; ///. only one initialization is needed
	bool writeVoltage=false; ///. only one initialization is needed


		if(headerName.size()<2 || headerName.substr(0,2)=="-h")
			Info<<"\nUsage:"<<endl
				<< "  FOAM2Voxel  headerFileName  N_processors  output-format " << endl<< endl
				<< "  output-format can be \"ascii\", \"binary\"(default), \"all\" or \"oldAscii\" or \"oldBinary\"." << endl
				<< "  Works for both serial (N_processors = 1, default) and parallel runs." << endl
				<< "  \"oldAscii\" output-format is for compatiblity with old IC dispression codes. " << endl<< endl;
		if(headerName.empty() || headerName.substr(0,2)=="-h")
		 { if(headerName.empty()) Info<<"Error: please try again  providing headerFileName"<<endl; return 1;	}	


	if(argc>2) nProcs = atoi(argv[2]);
	if(argc>3) outputFormat = std::string (argv[3]); 
	if(argc>4) writeVoltage = argv[4][0]=='t' || argv[4][0]=='T'; 
	//if(argc>4) skipOutlet = std::string(argv[4]); 
	imgExt(outputFormat);
	Info<<"FOAM2Voxel "<<headerName<<" "<<nProcs<<" "<<outputFormat<<" "<<writeVoltage<<endl;

	voxelImage vximage(headerName,2);
	int3 n = vximage.size3();
	dbl3 xmin=vximage.X0();
	dbl3 dx=vximage.dx();
	vximage.reset(n,1);

	Info<<"\nn: "<<n.x<<" "<<n.y<<" "<<n.z<<" \n"<<endl;

	//dx*=1e-6; 
	//xmin*=1e-6; 


  std::ofstream("OpenMeInParaview.xmf")<<
		"<?xml version=\"1.0\" ?>\n"
		"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n"
		"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">\n"
		" <Domain>\n"
		"  <Grid GridType=\"Uniform\">\n"
		" 	<Topology TopologyType=\"3DCORECTMesh\" Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<"\"/>\n"
		" 	<Geometry GeometryType=\"ORIGIN_DXDYDZ\">\n"
		" 		<DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\"> 0 0 0 </DataItem>\n"
		" 		<DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\"> "<<dx[2]<<" "<<dx[1]<<" "<<dx[0]<<" </DataItem>\n"
		" 	</Geometry>\n"
		" 	<Attribute Name=\"rock\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
		" 		<DataItem Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<"\" DataType=\"UChar\" Precision=\"4\" Format=\"Binary\">	vxlImage"+imgExt()+ " </DataItem>\n"
		" 	</Attribute>\n"
		" 	<Attribute Name=\"pressure\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">\n"
		" 	  <DataItem Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	p"+imgExt()+" </DataItem>\n"
		" 	</Attribute>\n"
		" 	<Attribute Name=\"Uf\" Active=\"0\" AttributeType=\"Vector\" Center=\"Node\">\n"
		" 		<DataItem Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<" 3\" ItemType=\"Function\"  Function=\"JOIN($0 , $1, $2)\" >\n"
		" 		 <DataItem ItemType=\"HyperSlab\"  Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<"\"  Type=\"HyperSlab\">\n"
		" 			 <DataItem Dimensions=\"3 3\" Format=\"XML\">   0 0 0    1 1 1   "<<n.z<<" "<<n.y<<" "<<n.x<<"   </DataItem>\n"
		" 			 <DataItem Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x+1<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	Ufx"+imgExt() +" </DataItem>\n"
		" 		 </DataItem>\n"
		" 		 <DataItem ItemType=\"HyperSlab\"  Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<"\"  Type=\"HyperSlab\">\n"
		" 			 <DataItem Dimensions=\"3 3\" Format=\"XML\">   0 0 0    1 1 1   "<<n.z<<" "<<n.y<<" "<<n.x<<"   </DataItem>\n"
		" 			 <DataItem Dimensions=\""<<n.z<<" "<<n.y+1<<" "<<n.x<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	Ufy"+imgExt() +" </DataItem>\n"
		" 		 </DataItem>\n"
		" 		 <DataItem ItemType=\"HyperSlab\"  Dimensions=\""<<n.z<<" "<<n.y<<" "<<n.x<<"\"  Type=\"HyperSlab\">\n"
		" 			 <DataItem Dimensions=\"3 3\" Format=\"XML\">   0 0 0    1 1 1   "<<n.z<<" "<<n.y<<" "<<n.x<<"   </DataItem>\n"
		" 			 <DataItem Dimensions=\""<<n.z+1<<" "<<n.y<<" "<<n.x<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	Ufz"+imgExt() +" </DataItem>\n"
		" 	    </DataItem>\n"
		" 	  </DataItem>\n"
		" 	</Attribute>\n"
		"  </Grid>\n"
		" </Domain>\n"
		"</Xdmf>\n"
		<<std::endl;




	bool saveMemory=true;
	try  {
		size_t testSize = size_t(n.x)*n.y*n.z*4l*4l*(1.+min(4./nProcs,1.));
		cout <<"testing memory, testSize (GB): 4 x "<< testSize/1000000000 <<std::endl ;
		char* testmemory1 = static_cast<char*>(malloc(testSize));
		char* testmemory2 = static_cast<char*>(malloc(testSize));
		char* testmemory3 = static_cast<char*>(malloc(testSize));
		char* testmemory4 = static_cast<char*>(malloc(testSize));
		if(testmemory4) {saveMemory = false; free(testmemory4); }
		if(testmemory3) { free(testmemory3); }
		if(testmemory2) { free(testmemory2); }
		if(testmemory1) { free(testmemory1); }
		Info<<"\nsave memory: "<<int(saveMemory)<<endl;
	}
	catch (std::bad_alloc &ba)  {
	 cout<<"allocation failure!  switching to save-memory version" ;
	}

	if (!saveMemory)  {
		voxelField<float> pVoxel,peVoxel,vface0,vface1,vface2;
		pVoxel.reset(n,0.);
		if(writeVoltage)  peVoxel.reset(n,0.);
		vface0.reset(n.x+1,n.y,  n.z,  0.);
		vface1.reset(n.x,  n.y+1,n.z,  0.);
		vface2.reset(n.x,  n.y,  n.z+1,0.);


		for (int p=0;p<nProcs;p++)  {
			int argcProc = 3;
			string caseName("./processor"+_s(p)+"\0");
			if (nProcs==1) caseName=".\0";
			char *argvProc[3]={argv[0], _case,&(caseName[0u])};

			
			Info<< "\nprocessor: "<<p<<"========================================" << endl;

			procMain( argcProc, argvProc, vximage, pVoxel, vface0, vface1, vface2, peVoxel, initialiseOF);
			initialiseOF=false;
		}


		if (outputFormat[0]=='o')  {
			 vface0.writeAscii("./Ux.dat", 0,n.x, 0,n.y, 0,n.z);
			 vface1.writeAscii("./Uy.dat", 0,n.x, 0,n.y, 0,n.z);	
			 vface2.writeAscii("./Uz.dat", 0,n.x, 0,n.y, 0,n.z);	
			 pVoxel.writeAscii("./p.dat");
			 if(writeVoltage) peVoxel.writeAscii("./psi.dat");
			 vximage.writeAscii("./vxlImage.dat");
			 vximage.writeHeader("./vxlImage.dat");
		}
		else  {
			vface0.writeNoHdr("./Ufx"+imgExt());
			vface1.writeNoHdr("./Ufy"+imgExt());	
			vface2.writeNoHdr("./Ufz"+imgExt());	
			pVoxel.writeNoHdr("./p"+imgExt());
			if(writeVoltage) peVoxel.writeNoHdr("./psi"+imgExt());
			vximage.writeNoHdr("./vxlImage"+imgExt());
			vximage.writeHeader("./vxlImage-"+imgExt());
		}


		Info<< "end" << endl;

		vximage.printInfo();
	}
	else  {
		{///. vximage

			for (int p=0;p<nProcs;p++)	{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainV( 3, argvProc, vximage, n, xmin, dx, initialiseOF);
				initialiseOF=false;
			}


			if (outputFormat[0]=='o')	{//old
				 vximage.writeAscii("./vxlImage.dat");
				 vximage.writeHeader("./vxlImage.dat_header");
			}
			else  {
				 vximage.write("./vxlImage"+imgExt());
				 vximage.writeHeader("./vxlImage"+imgExt());
			}
			
			vximage.printInfo();
			vximage.reset(0,0,0,0);
		}

		{	voxelField<float> pVoxel(n,0.);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainP( 3, argvProc, pVoxel, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				 pVoxel.writeAscii("./p.dat");
			else
				 pVoxel.writeNoHdr("./p"+imgExt());
		}

	///. voltage / concentrat
		if(writeVoltage)  {
			voxelField<float> pVoxel(n,0.);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainPE( 3, argvProc, pVoxel, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				 pVoxel.writeAscii("./psi.dat");
			else
				 pVoxel.writeNoHdr("./psi"+imgExt());
		}

		{	voxelField<float> vface0(n.x+1,n.y,n.z,0.);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainUx( 3, argvProc, vface0, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				vface0.writeAscii("./Ux.dat", 0,n.x, 0,n.y, 0,n.z);
			else
				vface0.writeNoHdr("./Ufx"+imgExt());
		}

		{	voxelField<float> vface1(n.x,n.y+1,n.z,0.);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainUy( 3, argvProc, vface1, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				vface1.writeAscii("./Uy.dat", 0,n.x, 0,n.y, 0,n.z);	
			else
				vface1.writeNoHdr("./Ufy"+imgExt());	
		}

		{	voxelField<float> vface2(n.x,n.y,n.z+1,0.);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainUz( 3, argvProc, vface2, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				vface2.writeAscii("./Uz.dat", 0,n.x, 0,n.y, 0,n.z);	
			else
				vface2.writeNoHdr("./Ufz"+imgExt());	
		}



	}



	return 0;
}


// ************************************************************************* //
