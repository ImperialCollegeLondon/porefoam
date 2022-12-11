
all:  ;
clean:;


tsts = $(wildcard *.mhd) $(wildcard *.sh)

$(warning tsts: ${tsts})

msSrc ?= $(abspath ../..)

USE_msTEST=1
include  ${msSrc}/script/Makefile.in
