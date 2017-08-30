#makefile for asymptotic theory calculations
OBJS= molecularDiff.Rev02.o 	\
	tynCalus.o 	\
	AndradeVSL.o 	\
	rhoCalc.o 	\
	TsatAntoinne.o 	\
	linspace.o 	\
	DABuncertainty.o

FC= gfortran 

molecularDiff: $(OBJS)
	$(FC) -o molecularDiff $(OBJS)

molecularDiff.Rev02.o:molecularDiff.Rev02.f90
	$(FC) -c molecularDiff.Rev02.f90

tynCalus.o:tynCalus.f90
	$(FC) -c tynCalus.f90

AndradeVSL.o:AndradeVSL.f90
	$(FC) -c AndradeVSL.f90

rhoCalc.o:rhoCalc.f90
	$(FC) -c rhoCalc.f90	

TsatAntoinne.o:TsatAntoinne.f90
	$(FC) -c TsatAntoinne.f90	

linspace.o:linspace.f90
	$(FC) -c linspace.f90

DABuncertainty.o:DABuncertainty.f90
	$(FC) -c DABuncertainty.f90

clean:
	rm $(OBJS)

#End of makefile