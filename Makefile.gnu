# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2019.0.117/linux/mkl/
FC90 = /usr/sepp/bin/ifort-2020-af
FC90 = /usr/sepp/bin/gfortran-9.2.0-af
#F90_FLAGS =  -r8 -check bounds -traceback -fpp 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
F90_FLAGS =  -ffree-line-length-none -fopenmp
LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_lp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl 


# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: liteQT_g.x 

liteQT_g.x : main.f90 mkl_dfti.o wannierHam.o green.o green_rgf.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -I$(MODDIR) -J$(MODDIR)


.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -I$(MODDIR) -J$(MODDIR) 
