BINDIR = ./Bin
FIGDIR = ./Fig
SRCDIR = ./Src

reproduce:
	! pushd ${FIGDIR} ; make ; popd
 
clean:
	! pushd ${SRCDIR} ; make clean; popd
	! pushd ${BINDIR} ; make clean; popd
	! pushd ${FIGDIR} ; make clean; popd
