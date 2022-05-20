# Intel Fortran
#FORT = ifort
#EXE = dft
#FLAGS = -O3 -no-wrap-margin

# GNU Fortran
FORT = gfortran
EXE = dft
FLAGS = -Ofast -ffree-line-length-none

all:
	$(FORT) $(FLAGS) $(EXE).f90 -o $(EXE).exe

clean: 
	rm $(EXE).exe
	rm *.mod

roda: all
	./$(EXE).exe
	python3 plots.py
