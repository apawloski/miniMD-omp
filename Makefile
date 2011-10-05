# Multiple-machine Makefile

SHELL = /bin/sh

# Files

SRC =	ljs.cpp input.cpp integrate.cpp atom.cpp force.cpp neighbor.cpp \
	thermo.cpp comm.cpp timer.cpp output.cpp setup.cpp
INC =	ljs.h atom.h force.h neighbor.h thermo.h timer.h comm.h integrate.h

# Definitions

ROOT =	miniMD
EXE =	$(ROOT)_$@
OBJ =	$(SRC:.cpp=.o)

# Help

help:
	@echo 'Type "make target" where target is one of:'
	@echo '      mpicxx     (for cygwin box with mpicxx compiler)'
	@echo '      openmpi    (for linux box with OpenMPI and gcc 4.x compiler)'
	@echo '      openmpi-omp    (for linux box with OpenMPI+OpenMP and gcc 4.x compiler)'

# Targets

mpicxx openmpi openmpi-omp:
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE)  "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)

# Clean

clean: clean_mpicxx clean_openmpi clean_openmpi-omp

clean_mpicxx:
	rm -rf Obj_mpicxx

clean_openmpi:
	rm -rf Obj_openmpi

clean_openmpi-omp:
	rm -rf Obj_openmpi-omp

