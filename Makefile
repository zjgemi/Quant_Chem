#####################################################################
#
# Typing the command below => results in the following actions
#    make             =>      Executable a.out
#    make clean       =>      Remove object files during compilation
#    make realclean   =>      Remove object files, executables
#
######################################################################

FF = ifort
CFLAG = -O3
#CFLAG = -fast
#CFLAG = -mcmodel=medium
#PAL = -D PARALLEL
LIB = -mkl

FILES = constant.o gaussian.o gaussian-stype.o SCF.F90

all: scf.out
scf.out: constant.o gaussian.o gaussian-stype.o
	$(FF) $(LIB) $(FILES) -o $@

clean:
	rm -f *.o *.out
realclean:
	rm -f *.o *.out *.e* *.o*

.SUFFIXES:

constant.o: constant.f90
	$(FF) -c constant.f90

gaussian.o: gaussian.f90
	$(FF) -c gaussian.f90

gaussian-stype.o: gaussian-stype.f90
	$(FF) -c gaussian-stype.f90
