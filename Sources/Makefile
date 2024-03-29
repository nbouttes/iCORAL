.SUFFIXES: .f .F .F90 .f90 .o .mod
.SHELL: /bin/sh

.PHONY : usage
usage:
	@echo ""
	@echo "    * USAGE * "
	@echo ""
	@echo " make all        : compiles the test program main"
	@echo ""


objdir = .

# Command-line options at make call
debug ?= 1

## GFORTRAN OPTIONS ##
FC = ifort
# modify here

NETCDFINCLUDE	= -I/usr/local/install/netcdf-4.3.2/include
RPATHNETCDF     = /usr/lib64
NETCDFLIB	    = -L/usr/lib64 -lnetcdf -lnetcdff -L/lib64
NCIO_PATH       = /home/acclimate/nbouttes/iloveclim_degla3
NCIOLIB         = -L$(NCIO_PATH)/ncio/library -lncio
NCIOINCLUDE     = -I$(NCIO_PATH)/ncio/.obj

LIB = -Wl,-rpath,$(RPATHNETCDF) $(NETCDFLIB) $(NCIOLIB)
INC = $(NETCDFINCLUDE) $(NCIOINCLUDE)

FLAGS  = -I$(objdir) $(INC)
LFLAGS = $(LIB)
MYFLAGS = -mcmodel large

DFLAGS = -O3
ifeq ($(debug), 1)
    DFLAGS   = -C -fpe0 -g -p
endif

## Individual libraries or modules ##
$(objdir)/coral_mod.o: coral_mod.f90
	$(FC) $(DFLAGS) $(FLAGS) -c -o $@ $<

## Complete programs

all: $(objdir)/coral_mod.o
	$(FC) $(DFLAGS) $(FLAGS) $(MYFLAGS) -o main $^ main.f90 $(LFLAGS)
	@echo " "
	@echo "    main @ ok_coral is ready."
	@echo " "

clean:
	rm -f main $(objdir)/*.o $(objdir)/*.mod
