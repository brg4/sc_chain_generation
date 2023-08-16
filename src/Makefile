# Makefile for program gen_sc_chain
F90	= ${FC}

OPTIONS	= -fdump-core -fcheck=all -Wunused -Wuninitialized -fbacktrace -fopenmp -lgomp -lgfortran -O3

OBJECTS = params.o quaternions.o sc_chain_object.o sc_chain_procedures.o read_write.o gen_sc_chain.o
PROGRAM	= gen_sc_chain 	#target- executable that is generated
ARCHIVE = gen_sc_chain

help:
	@echo Try:
	@echo make $(PROGRAM) .... to build the program named $(PROGRAM)
	@echo make clean .... to clean up, removing object files and program $(PROGRAM)
	@echo make listing .... to make a printable listing ... you may want to adjust the page length
	@echo make archive .... to make an archive tar file you can transfer or submit


%.o:	%.f90
	$(F90) $(OPTIONS) -c $<

$(PROGRAM):	$(OBJECTS)	
	$(F90) $(OBJECTS) $(OPTIONS) -o $(PROGRAM)	

clean:
	rm $(OBJECTS) $(PROGRAM)
