
msSrc ?= ..

all:
	$(warning "running ${msSrc}/script/AllMake on $(shell dirname $(firstword $(MAKEFILE_LIST)))")
	${msSrc}/script/AllMake


clean:
	$(warning "running ${msSrc}/script/AllClean on $(shell dirname $(realpath $(firstword $(MAKEFILE_LIST))))")
	${msSrc}/script/AllClean



testDir = ${msSrc}/../test/Porefoam1f
test:
	mkdir -p ${testDir}
	$(msSrc)/script/testApp "${testDir}"  python3 "${CURDIR}/test.py"
