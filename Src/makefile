#
# SEPLIB directories
#
SEPLIBDIR = /usr/local/SEP/lib
SEPINCDIR = /usr/local/SEP/include
#
# Source files and Binary files directories
#
BINDIR    = ./../Bin
SRCDIR    = .
#
# Unix commands
#
INSTALL   = /usr/bin/install -c -m 775
RM        = /bin/rm
FPP       = /lib/cpp -P -traditional
#--------------#
# F90 DEFAULTS #
#--------------#
F90C              = ifort -I.
F90LD             = ifort -L${SEPLIBDIR}
F90EXT            = f90
F90OLEVEL         =
F90OPTIMFLAGS     = -O${F90OLEVEL}
OF90DEFINES  	  = -I. -I${SEPINCDIR}
OF90FLAGS         = -g -FR  -I${SEPINCDIR} -c
OF90LIBS          = -Vaxlib -lm
F90LDFLAGS        = $(filter %.o,$^)
SEPLIB3DF90       = -lsep3df90 -lsep3d -lsepf90 -lsep
F90LIBS           = ${SEPLIB3DF90} ${OF90LIBS}

define F90_COMPILE
	${F90C} ${OF90FLAGS} ${F90OPTIMFLAGS} -o $@
endef

define F90_LN
	${F90LD} ${F90LDFLAGS} ${F90LIBS}
	${INSTALL} a.out $@
	${RM} a.out
endef

%.o: %.${F90EXT}
	${F90_COMPILE}  $*.${F90EXT}

${BINDIR}/%.x %.x: %.o
	${F90_LN}

#
#------------------------------------------------------------
#

default: ${BINDIR}/test_dot1d.x ${BINDIR}/test_dot2d.x ${BINDIR}/test_rt1d.x ${BINDIR}/test_rt2d.x ${BINDIR}/Rt1dinv.x ${BINDIR}/Rt2dinv.x ${BINDIR}/test_rt1d_conv.x ${BINDIR}/test_rt2d_conv.x

#
#	dot1d : Dot-product test for 1D RTM algorithm (rt1d.rt)
#	dot2d : Dot-product test for 2D RTM algorithm (rt2d.rt)
#
dot1d:	${BINDIR}/test_dot1d.x
	${BINDIR}/test_dot1d.x nodoc

dot2d:	${BINDIR}/test_dot2d.x
	${BINDIR}/test_dot2d.x nodoc

${BINDIR}/test_dot1d.x: test_dot1d.o rt1d.o subs.o

${BINDIR}/test_dot2d.x: test_dot2d.o rt2d.o rt2d_conv.o subs.o

${BINDIR}/test_rt1d.x: test_rt1d.o rt1d.o subs.o

${BINDIR}/test_rt1d_conv.x: test_rt1d_conv.o rt1d_conv.o subs.o 

${BINDIR}/Rt1dinv.x: Rt1dinv.o rt1dinv.o rt1dinv2.o rt1d_conv.o rt1d.o subs.o

${BINDIR}/test_rt2d.x: test_rt2d.o rt2d.o subs.o

${BINDIR}/test_rt2d_conv.x: test_rt2d_conv.o rt2d_conv.o subs.o 

${BINDIR}/Rt2dinv.x: Rt2dinv.o rt2dinv.o rt2dinv2.o rt2d_conv.o rt2d.o subs.o

clean:
	/bin/rm -f *.o a.out core
