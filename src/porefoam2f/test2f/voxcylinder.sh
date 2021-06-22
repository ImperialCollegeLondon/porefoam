#!/bin/bash


printf "DimSize = 32 32 32 \n\
	Offset =      0    0    0 \n\
	replaceRange 0 255 1 \n\
	reset  dx 1 1 1 \n\
	Paint cylinder 0 16 16   32 16 16  8 \n\
	" > voxcylinder.mhd
	#ElementDataFile = NO_READreset  dx 1e-6 1e-6 1e-6
	#nIterationsGaussVP   8 \n\
	#kernelRadiusGaussVP  4 \n\

export RefineLevel=0.90  #< this should be less than one for cfMesh to work, check!
export endTime=0.001 #< consider increasing this to 0.01
export UD1s=1000
export nProc=4
rm -rf voxcylinder*/ AllRunImageTwoPhase 
AllRunImageTwoPhase voxcylinder.mhd  # this just copies AllRunImageTwoPhase
AllRunImageTwoPhase voxcylinder.mhd 
AllRunImageTwoPhase voxcylinder.mhd 
AllRunImageTwoPhase voxcylinder.mhd 



#exit `python3 -c  "from msrc import *; print(fileFloatDiffersFrom('voxcylinder/voxcylinder-1-X/summary_voxcylinder-1-X.txt','K_x= ',2.44436e-12))"`
