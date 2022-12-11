

all: msLib

msSrc ?= $(abspath ../..)

USE_MINIFOAMX3=1
OFCs=$(shell find | grep -F .C | sed s:./::)
OFincs =   meshTools finiteVolume transportModels/incompressible transportModels
OFlibs =   foam  finiteVolume meshTools incompressibleTransportModels

msCFLAGS += -IalphaContactAngle -fPIC  -DNoRepository
msLFLAGS +=  -fPIC -pthread -Wl,--enable-new-dtags
USE_ZLIB=1
USE_TIFF=1
USE_msMAKE=1
include  ${msSrc}/script/Makefile.in
