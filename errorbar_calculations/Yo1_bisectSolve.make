#makefile to compile fortran code used to calculate corrections
#for Yo. A shared object is created which can be called from within R
OBJS=	Yo1_bisectSolve.o

FC= gfortran 

Yo1_bisectSolve.so: $(OBJS)
	$(FC) -shared -o Yo1_bisectSolve.so $(OBJS)
	
Yo1_bisectSolve.o:Yo1_bisectSolve.f90
	$(FC) -c Yo1_bisectSolve.f90

clean:
	rm $(OBJS)

#End of makefile