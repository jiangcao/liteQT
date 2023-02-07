# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2019.0.117/linux/mkl
FC90 = /usr/sepp/bin/ifort-2020-af
#F90_FLAGS =  -r8 -check bounds -traceback -fpp 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
F90_FLAGS =  -r8 -O2 -fpp -mkl -traceback -qopenmp # -check bounds 
LIBS = -I$(MKLROOT)/include


# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: liteQT.x 

liteQT.x : main.f90 mod_string.o wannierHam.o green.o bse_mod.o phononDyn.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) $(LIBARPACK) -module $(MODDIR)

.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -module $(MODDIR)
