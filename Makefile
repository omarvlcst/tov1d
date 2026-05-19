FC = gfortran
OPTIONS = -O3 -fPIC -finit-local-zero

Solver.out: main.f90 consts.o grid.o rk4.o funciones.o
	$(FC) $(OPTIONS) -o Solver.out main.f90 consts.o grid.o rk4.o funciones.o

consts.o: consts.f90
	$(FC) $(OPTIONS) -c consts.f90

grid.o: grid.f90
	$(FC) $(OPTIONS) -c grid.f90

rk4.o: rk4.f90
	$(FC) $(OPTIONS) -c rk4.f90

funciones.o: funciones.f90
	$(FC) $(OPTIONS) -c funciones.f90

clean:
	rm -rf *.o *.out *.mod *.dat