
msSrc ?= ..

ifneq "${OPT}" ".exe"
all:  ; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseMake USE_msRecurse=1
else
all: ; @echo skipping ${CURDIR}
endif

clean:; $(MAKE) -f  ${msSrc}/script/Makefile.in  recurseClean USE_msRecurse=1


testDir = ${msSrc}/../test/Porefoam2f
test: ;	cd test2f && make test
#	mkdir -p ${testDir}
#	$(msSrc)/script/testApp "${testDir}"  python3 "${CURDIR}/test.py"
