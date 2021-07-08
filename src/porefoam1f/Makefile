
msSrc ?= ..

all:  ; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseMake USE_msRecurse=1
clean:; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseClean USE_msRecurse=1

testDir = ${msSrc}/../test/Porefoam1f
test:
	mkdir -p ${testDir}
	$(msSrc)/script/testApp "${testDir}"  python3 "${CURDIR}/test.py"
