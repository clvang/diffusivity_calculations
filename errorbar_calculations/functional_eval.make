#makefile for asymptotic theory calculations. this creates
#shared objects between R and Fortran
OBJS=	functional_eval.o \
		Fx_Eval.o 	

FC= gfortran 

functional_eval.so: $(OBJS)
	$(FC) -shared -o functional_eval.so $(OBJS)

functional_eval.o: functional_eval.f90
	$(FC) -c functional_eval.f90
	
Fx_Eval.o:Fx_Eval.f90
	$(FC) -c Fx_Eval.f90

clean:
	rm $(OBJS) 

#End of makefile