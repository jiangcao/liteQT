# Fortran 90 compiler
MKLROOT = /usr/pack/intel_compiler-2020-af/x64/compilers_and_libraries_2019.0.117/linux/mkl/
FC90 = /usr/sepp/bin/ifort-2020-af
#F90_FLAGS =  -r8 -check bounds -traceback -fpp 
MACFLAGS = -L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib
#F90_FLAGS =  -I"${MKLROOT}/include" -r8 -O2 -fpp -mkl -traceback  -qopenmp -qopt-matmul # -check bounds  

F90_FLAGS =  -r8 -O2 -fpp -mkl -traceback  -qopenmp -qopt-matmul -heap-arrays # -check bounds  

#LIBS = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_sequential.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm -ldl

# Modules directory
MODDIR = compiled

# Source directory
SRCDIR = src/

# Search directories
vpath %.f90 $(SRCDIR)
vpath %.o $(MODDIR)

# Targets.

all: liteQT_3d.x 

liteQT_3d.x : main3d.f90 mkl_dfti.o wannierHam3d.o green.o green_rgf.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -module $(MODDIR)


liteQT_n.x : main.f90 mkl_dfti.o wannierHam.o green.o green_rgf.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -module $(MODDIR)

bse.x : bse_abs.f90 wannierHam.o bse_mod.o

	$(FC90) -o $@ $< $(MODDIR)/*.o $(F90_FLAGS) $(LIBS) -module $(MODDIR)


.PHONY : clean;

clean :
	rm -f $(MODDIR)/*.o $(MODDIR)/*.mod

# implicit rules

%.o : %.f90
	$(FC90) -o $(MODDIR)/$@ $< -c $(F90_FLAGS) -module $(MODDIR)
