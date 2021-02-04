# Makefile for F90 version of LSMR and zLSQR.
# Maintained by Michael Saunders <saunders@stanford.edu>
#
# 16 Jul 2010: LSMR version derived from LSQR equivalent.
# 29 Jun 2013: Support for complex LSQR (zLSQR)

  FC      =  gfortran
  FFLAGS  = -g -O
# FFLAGS  = -g -O0 -pedantic -Wall -Wextra -fbounds-check -ftrace=full

  CC      =  gcc
  CFLAGS  = -g -O

# Clear suffix list, then define the ones we want
  .SUFFIXES:
  .SUFFIXES: .c .f .f90 .o

  .f90.o:; ${FC} ${FFLAGS} -c -o $@ $<
  .f.o:;   ${FC} ${FFLAGS} -c -o $@ $<
  .c.o:;   $(CC) $(CFLAGS) -c -o $@ $<

  files = lsmrDataModule.o lsmrblas.o        lsmrblasInterface.o \
          lsmrModule.o     lsmrCheckModule.o lsmrTestModule.o    \
          lsmrTestProgram.o

  zfiles = lsmrDataModule.o zlsmrDataModule.o  zlsmrblas.o  \
           mm_ioModule.o    zlsmrblasInterface.o lsmrReadMtxModule.o \
           zlsmrModule.o    zlsmrCheckModule.o zlsmrTestModule.o    \
           zlsmrTestProgram.o

all: lsmr zlsmr

lsmr: ${files}
	${FC} ${FFLAGS} -o TestProgram ${files}

zlsmr: ${zfiles}
	${FC} ${FFLAGS} -o zTestProgram ${zfiles}

clean:
	\rm -f *.o *.mod
