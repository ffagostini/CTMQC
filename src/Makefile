FC = gfortran
CFLAG =  -Wall -pedantic -std=f95 -O3 -fopenmp -Wmaybe-uninitialized -fcheck=all
#CFLAG =  -Wall -pedantic -std=f95 -O3 -fopenmp -Wmaybe-uninitialized 
LFLAG = -llapack -lblas

.SUFFIXES : .o .f90 .f
.f90.o:
	$(FC) $(CFLAG) -c $<
.f.o:
	$(FC) $(CFLAG) -c $<

%.o: %.mod

MAIN = main.o

MODS = \
linear_algebra.o \
read_diabatic_properties.o \
read_BOfiles.o \
classical_evolution.o \
coefficients_evolution.o \
electronic_problem.o \
coherence_corrections.o \
kinds.o \
output.o \
time_evolution.o \
tools.o \
variables.o \
wigner_distribution.o

all:  main.x


main.x: $(MAIN) $(MODS) $(SUBS)
	$(FC) $(CFLAG) -o $@ $(MAIN) $(MODS) $(SUBS) \
	$(JDLIB) $(WSMPLIB) $(LFLAG)

.PHONY: clean
clean:
	rm -f *.o *.mod *.x

include make.depends
