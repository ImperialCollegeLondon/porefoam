#!/bin/bash


printf "DimSize = 22 32 32 \n\
	Offset =      0    0    0 \n\
	replaceRange 0 255 1 \n\
	reset  dx 1 1 1 \n\
	Paint cylinder 0 16 16   22 16 16  8 \n\
	" > voxcyl32c2f.mhd
	#ElementDataFile = NO_READreset  dx 1e-6 1e-6 1e-6
	#nIterationsGaussVP   8 \n\
	#kernelRadiusGaussVP  4 \n\

export RefineLevel=0.90  #< this should be less than one for cfMesh to work, check!
export endTime=0.001 #< consider increasing this to 0.01
export UD1s=1000
export nProc=4
export oilFilldFracs=0.5

export smoothMeshNIter=6
export smoothMeshBoundary=yes
export smoothMeshInternal=yes
export smoothMeshRelax=0.1

rm -rf voxcyl32c2f*/ AllRunImageTwoPhase
AllRunImageTwoPhase voxcyl32c2f.mhd  # this just copies AllRunImageTwoPhase
AllRunImageTwoPhase voxcyl32c2f.mhd 
AllRunImageTwoPhase voxcyl32c2f.mhd 
AllRunImageTwoPhase voxcyl32c2f.mhd 



#exit `python3 -c  "from msrc import *; print(fileFloatDiffersFrom('voxcyl32c2f/voxcyl32c2f-1-X/summary_voxcyl32c2f-1-X.txt','K_x= ',2.44436e-12))"`
