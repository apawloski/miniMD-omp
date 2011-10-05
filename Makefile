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

# Targets

mpicxx:
	@if [ ! -d Obj_$@ ]; then mkdir Obj_$@; fi
	@cp -p $(SRC) $(INC) Obj_$@
	@cp Makefile.$@ Obj_$@/Makefile
	@cd Obj_$@; \
	$(MAKE)  "OBJ = $(OBJ)" "INC = $(INC)" "EXE = ../$(EXE)" ../$(EXE)
#	@if [ -d Obj_$@ ]; then cd Obj_$@; rm $(SRC) $(INC) Makefile*; fi

# Clean

clean_mpicxx:
	rm -r Obj_mpicxx

