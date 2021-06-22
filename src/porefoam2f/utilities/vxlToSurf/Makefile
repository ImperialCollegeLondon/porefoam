# wmake  

all: msApp


# Makefile.in order is important
USE_VXIMG=1
USE_OMP=1
#include  ../Makefile.in

msSrc ?= $(abspath ../../..)

msCFLAGS := ${msCFLAGS}   -O3   -DGCC_COMPILER 




srcs=$(shell find | grep -F .cpp | sed s:./::)
srcs+=${msSrc}/libvoxel/voxelImage.cpp

USE_CPP17=1
USE_ZLIB=1
USE_TIFF=1
USE_msMAKE=1
include  ${msSrc}/script/Makefile.in

