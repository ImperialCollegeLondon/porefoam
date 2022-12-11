

all: msApp

msSrc ?= $(abspath ../../..)

USE_MINIFOAMX3=1
OFCs=$(shell find | grep -F .C | sed s:./::)
OFincs =  meshTools  surfMesh
OFlibs =  meshTools  surfMesh

USE_ZLIB=1
USE_TIFF=1
USE_OMP=1
USE_msMAKE=1
include  ${msSrc}/script/Makefile.in
