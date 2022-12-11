#!/bin/bash


source $(dirname "${BASH_SOURCE[0]}")/../../script/bashrc

printf "DimSize = 22 32 32 \n\
	Offset =      0    0    0 \n\
	replaceRange 0 255 1 \n\
	reset  dx 1 1 1 \n\
	Paint cylinder 0 16 16   22 16 16  8 \n\
	" > voxcyl32cfm.mhd
	#ElementDataFile = NO_READreset  dx 1e-6 1e-6 1e-6
	#nIterationsGaussVP   8 \n\
	#kernelRadiusGaussVP  4 \n\

export RefineLevel=0.90  #< this should be less than one for cfMesh to work, check!
export endTime=0.001 #< consider increasing this to 0.01
export UD1s=1000
export nProc=4
export oilFilldFracs=0.5
rm -rf voxcyl32cfm*/ AllRunImageTwoPhaseCFMesh
AllRunImageTwoPhaseCFMesh voxcyl32cfm.mhd  # this just copies AllRunImageTwoPhase
AllRunImageTwoPhaseCFMesh voxcyl32cfm.mhd
AllRunImageTwoPhaseCFMesh voxcyl32cfm.mhd
AllRunImageTwoPhaseCFMesh voxcyl32cfm.mhd



#exit `python3 -c  "from msrc import *; print(fileFloatDiffersFrom('voxcyl32cfm/voxcyl32cfm-1-X/summary_voxcyl32cfm-1-X.txt','K_x= ',2.44436e-12))"`
