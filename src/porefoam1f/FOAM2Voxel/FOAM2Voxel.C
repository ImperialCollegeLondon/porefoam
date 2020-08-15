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
* 
Please see our website for relavant literature:
http://www3.imperial.ac.uk/earthscienceandengineering/research/perm/porescalemodelling

For further information please contact us by email:
Ali Q Raeini:    a.q.raeini@imperial.ac.uk
Branko Bijeljic: b.bijeljic@imperial.ac.uk
Martin J Blunt:  m.blunt@imperial.ac.uk
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
	vximage.reset(n[0],n[1],n[2],1);

	Info<<"\nn: "<<n[0]<<" "<<n[1]<<" "<<n[2]<<endl<<endl;

	//dx*=1.0e-6; 
	//xmin*=1.0e-6; 


  {std::ofstream  vxlResults("OpenMeInParaview.xmf");
	vxlResults<<"<?xml version=\"1.0\" ?>"<<std::endl;
	vxlResults<<"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>"<<std::endl;
	vxlResults<<"<Xdmf xmlns:xi=\"http://www.w3.org/2003/XInclude\" Version=\"2.2\">"<<std::endl;
	vxlResults<<"  <Domain>"<<std::endl;
	vxlResults<<"	 <Grid GridType=\"Uniform\">"<<std::endl;
	vxlResults<<"		<Topology TopologyType=\"3DCORECTMesh\" Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<"\"/>"<<std::endl;
	vxlResults<<"		<Geometry GeometryType=\"ORIGIN_DXDYDZ\">"<<std::endl;
	vxlResults<<"		  <DataItem Name=\"Origin\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">"<<std::endl;
	vxlResults<<"				0 0 0"<<std::endl;
	vxlResults<<"			</DataItem>"<<std::endl;
	vxlResults<<"			<DataItem Name=\"Spacing\" Dimensions=\"3\" NumberType=\"Float\" Precision=\"4\" Format=\"XML\">"<<std::endl;
	vxlResults<<"				"<<dx[2]<<" "<<dx[1]<<" "<<dx[0]<<" "<<std::endl;
	vxlResults<<"			</DataItem>"<<std::endl;
	vxlResults<<"		</Geometry>"<<std::endl;
	vxlResults<<"		<Attribute Name=\"rock\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">"<<std::endl;
	vxlResults<<"		  <DataItem Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<"\" DataType=\"UChar\" Precision=\"4\" Format=\"Binary\">	vxlImage"+imgExt()<<std::endl;
	vxlResults<<"		  </DataItem>"<<std::endl;
	vxlResults<<"		</Attribute>"<<std::endl;
	vxlResults<<"		<Attribute Name=\"pressure\" Active=\"1\" AttributeType=\"Scalar\" Center=\"Cell\">"<<std::endl;
	vxlResults<<"		  <DataItem Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	p"+imgExt()<<std::endl;
	vxlResults<<"		  </DataItem>"<<std::endl;
	vxlResults<<"		</Attribute>"<<std::endl;
	vxlResults<<"		<Attribute Name=\"Uf\" Active=\"0\" AttributeType=\"Vector\" Center=\"Node\">"<<std::endl;
	vxlResults<<"			<DataItem Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<" 3\" ItemType=\"Function\"  Function=\"JOIN($0 , $1, $2)\" >"<<std::endl;
	vxlResults<<"			 <DataItem ItemType=\"HyperSlab\"  Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<"\"  Type=\"HyperSlab\">"<<std::endl;
	vxlResults<<"				 <DataItem Dimensions=\"3 3\" Format=\"XML\">   0 0 0    1 1 1   "<<n[2]<<" "<<n[1]<<" "<<n[0]<<"   </DataItem>"<<std::endl;
	vxlResults<<"				 <DataItem Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]+1<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	Ufx"+imgExt() +" </DataItem>"<<std::endl;
	vxlResults<<"			 </DataItem>"<<std::endl;
	vxlResults<<"			 <DataItem ItemType=\"HyperSlab\"  Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<"\"  Type=\"HyperSlab\">"<<std::endl;
	vxlResults<<"				 <DataItem Dimensions=\"3 3\" Format=\"XML\">   0 0 0    1 1 1   "<<n[2]<<" "<<n[1]<<" "<<n[0]<<"   </DataItem>		  "<<std::endl;
	vxlResults<<"				 <DataItem Dimensions=\""<<n[2]<<" "<<n[1]+1<<" "<<n[0]<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	Ufy"+imgExt() +" </DataItem>"<<std::endl;
	vxlResults<<"		    </DataItem>"<<std::endl;
	vxlResults<<"			 <DataItem ItemType=\"HyperSlab\"  Dimensions=\""<<n[2]<<" "<<n[1]<<" "<<n[0]<<"\"  Type=\"HyperSlab\">"<<std::endl;
	vxlResults<<"				 <DataItem Dimensions=\"3 3\" Format=\"XML\">   0 0 0    1 1 1   "<<n[2]<<" "<<n[1]<<" "<<n[0]<<"   </DataItem>		  "<<std::endl;
	vxlResults<<"				 <DataItem Dimensions=\""<<n[2]+1<<" "<<n[1]<<" "<<n[0]<<"\" DataType=\"Float\" Precision=\"4\" Format=\"Binary\">	Ufz"+imgExt() +" </DataItem>"<<std::endl;
	vxlResults<<"		    </DataItem>"<<std::endl;
	vxlResults<<"		  </DataItem>"<<std::endl;
	vxlResults<<"		</Attribute>"<<std::endl;
	vxlResults<<"	 </Grid>"<<std::endl;
	vxlResults<<"  </Domain>"<<std::endl;
	vxlResults<<"</Xdmf>"<<std::endl;
	vxlResults<<""<<std::endl;
  }









	bool saveMemory=true;
	try
	{
		size_t testSize = size_t(n[0])*n[1]*n[2]*4l*4l*(1.0+min(4.0/nProcs,1.0));
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
	catch (std::bad_alloc &ba) 
	{
	 cout<<"allocation failure!  switching to save-memory version" ;
	}

	if (!saveMemory)
	{
		voxelField<float> pVoxel,peVoxel,vface0,vface1,vface2;
		pVoxel.reset(n[0],n[1],n[2],0.0);
		if(writeVoltage)  peVoxel.reset(n[0],n[1],n[2],0.0);
		vface0.reset(n[0]+1,n[1],n[2],0.0);
		vface1.reset(n[0],n[1]+1,n[2],0.0);
		vface2.reset(n[0],n[1],n[2]+1,0.0);



		for (int p=0;p<nProcs;p++)
		{
			int argcProc = 3;
			string caseName("./processor"+_s(p)+"\0");
			if (nProcs==1) caseName=".\0";
			char *argvProc[3]={argv[0], _case,&(caseName[0u])};

			
			Info<< "\nprocessor: "<<p<<"========================================" << endl;

			procMain( argcProc, argvProc, vximage, pVoxel, vface0, vface1, vface2, peVoxel, initialiseOF);
			initialiseOF=false;
		}


		if (outputFormat[0]=='o')
		{
			 vface0.writeAscii("./Ux.dat", 0,n[0], 0,n[1], 0,n[2]);
			 vface1.writeAscii("./Uy.dat", 0,n[0], 0,n[1], 0,n[2]);	
			 vface2.writeAscii("./Uz.dat", 0,n[0], 0,n[1], 0,n[2]);	
			 pVoxel.writeAscii("./p.dat");
			 if(writeVoltage) peVoxel.writeAscii("./psi.dat");
			 vximage.writeAscii("./vxlImage.dat");
			 vximage.writeHeader("./vxlImage.dat");
		}
		else
		{
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
	else
	{
		{///. vximage

			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainV( 3, argvProc, vximage, n, xmin, dx, initialiseOF);
				initialiseOF=false;
			}


			if (outputFormat[0]=='o')//old
			{
				 vximage.writeAscii("./vxlImage.dat");
				 vximage.writeHeader("./vxlImage.dat_header");
			}
			else
			{
				 vximage.write("./vxlImage"+imgExt());
				 vximage.writeHeader("./vxlImage"+imgExt());
			}
			
			vximage.printInfo();
			vximage.reset(0,0,0,0);
		}

		{	voxelField<float> pVoxel(n[0],n[1],n[2],0.0);
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
		if(writeVoltage)
		{
			voxelField<float> pVoxel(n[0],n[1],n[2],0.0);
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

		{	voxelField<float> vface0(n[0]+1,n[1],n[2],0.0);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainUx( 3, argvProc, vface0, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				vface0.writeAscii("./Ux.dat", 0,n[0], 0,n[1], 0,n[2]);
			else
				vface0.writeNoHdr("./Ufx"+imgExt());
		}

		{	voxelField<float> vface1(n[0],n[1]+1,n[2],0.0);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainUy( 3, argvProc, vface1, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				vface1.writeAscii("./Uy.dat", 0,n[0], 0,n[1], 0,n[2]);	
			else
				vface1.writeNoHdr("./Ufy"+imgExt());	
		}

		{	voxelField<float> vface2(n[0],n[1],n[2]+1,0.0);
			for (int p=0;p<nProcs;p++)
			{
				string caseName((nProcs==1) ? ".\0" : "./processor"+_s(p)+"\0");
				char *argvProc[3]={argv[0], _case,&(caseName[0u])};
				Info<< "\nprocessor: "<<p<<"========================================" << endl;
				procMainUz( 3, argvProc, vface2, n, xmin, dx, initialiseOF);
			}
			if (outputFormat[0]=='o')//old
				vface2.writeAscii("./Uz.dat", 0,n[0], 0,n[1], 0,n[2]);	
			else
				vface2.writeNoHdr("./Ufz"+imgExt());	
		}



	}



	return 0;
}


// ************************************************************************* //
