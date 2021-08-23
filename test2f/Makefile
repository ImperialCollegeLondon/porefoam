
#all: ;
all:  ;
clean:  ;
#all: test ;


tsts = $(wildcard *.mhd) $(wildcard *.sh)

$(warning tsts: ${tsts})

msSrc ?= $(abspath ../..)
#msTestDir ?= /data3/Ali/bench2/$(shell basename ${CURDIR})${OPT}

USE_msTEST=1
include  ${msSrc}/script/Makefile.in

