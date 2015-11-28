
#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
F90 = gfortran
#---------------------

#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
# GNU
F90FLAGS = -O 

PROGS = sim_pixcmd fit_pixcmd

COMMON = nrtype.o nr.o pixcmd_vars.o nrutil.o pixcmd_utils.o ran_state.o \
	ran1.o gammln.o poidev.o locate.o init_random_seed.o convolve.o hist_2d.o \
	func.o powell.o linmin.o mnbrak.o brent.o get_model.o myran.o \
	emcee_advance.o add_obs_err.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

sim_pixcmd : sim_pixcmd.o $(COMMON)
	$(F90) -o sim_pixcmd.exe sim_pixcmd.o $(COMMON)

fit_pixcmd : fit_pixcmd.o $(COMMON)
	$(F90) -o fit_pixcmd.exe fit_pixcmd.o $(COMMON)

#------------------------------------------------------------------------#

sim_pixcmd.o : nrtype.o pixcmd_vars.o nr.o pixcmd_utils.o 

fit_pixcmd.o : nrtype.o pixcmd_vars.o nr.o pixcmd_utils.o 


%.o : %.f90
	$(F90) $(F90FLAGS) -o  $@ -c $<

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


