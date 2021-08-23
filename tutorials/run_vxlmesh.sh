#!/bin/bash

# check
[[ `pwd` != *"src/"* ]] || ! echo "Please run me from outside src tree" && exit 1
which AllRunImageTwoPhaseCFMesh || source porefoam/src/script/bashrc
which AllRunImageTwoPhaseCFMesh || ! echo -e 'Please first `source .../src/script/bashrc`' && exit 1

# download the input image
if ! [ -e ketton.raw ]; then echo "Downloading ketton image"
 wget -O ketton.raw https://www.digitalrocksportal.org/projects/125/images/101255/download/ 
fi



#keep the unit as micrometers otherwise the AllRunImagePar needs to be edited
# Create ketton.mhd header 
cat  <<MHD > ketton.mhd
ObjectType =  Image
NDims =       3
ElementType = MET_UCHAR
ElementByteOrderMSB = False
ElementNumberOfChannels = 1
CompressedData = False
HeaderSize = 0
DimSize     = 365  255  225
ElementSize = 1e-6 1e-6 1e-6 
Offset      = 0    0    0
ElementDataFile  = ketton.raw

pore  1 2   // replace 0-1 with 0 and the rest with 1

resample 2 //< make image coarser for quicker testing of mesh generation
resample 2 //< make image coarser for quicker testing of mesh generation

FaceMedian06   10
PointMedian032 15  20
FaceMedian06   10
PointMedian032 15  20

keepLargest0


//Notes:
//	image name:    Ketton_rock_trapped_oil_Segmented_SSb
//	download link: https://doi.org/10.17612/P7D95F
MHD


export meshTag="Vxl" # set tag to avoid conflict with AllRunImageTwoPhase
echo "Echo mesh smoothing may degrade the mesh quality"
export smoothMeshNIter=6
export smoothMeshBoundary=yes
export smoothMeshInternal=yes
export smoothMeshRelax=0.05

# Clean up in case
rm -rf ketton$meshTag*/ AllRunImageTwoPhaseCFMesh

# AllRunImageTwoPhase copies itself to the current directory for edit and run
AllRunImageTwoPhase 

# generate mesh, using voxelToFoam
./AllRunImageTwoPhase


./AllRunImageTwoPhase
./AllRunImageTwoPhase
