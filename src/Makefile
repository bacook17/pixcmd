
#======================
# Fortran 90 compiler
# (Uncomment only one)
#======================
# GNU
#F90 = gfortran
F90 = mpifort
#F90 = mpif90
#---------------------

#========================
# Compiler Optimizations
#  (Uncomment only one)
#========================
# GNU
F90FLAGS = -O

PROGS = fit_pixcmd dohess write_a_model

COMMON = nrtype.o nr.o pixcmd_vars.o nrutil.o pixcmd_utils.o ran_state.o \
	ran1.o gammln.o poidev.o locate.o init_random_seed.o convolve.o \
	hist_2d.o func.o powell.o linmin.o mnbrak.o brent.o getmodel.o \
	myran.o emcee_advance_mpi.o add_obs_err.o setup_models.o tsum.o \
	four2.o rlft2.o fourrow.o gasdev.o function_parallel_map.o \
	free_workers.o fit_oneatatime.o drawn.o mypoidev.o

all : $(PROGS)

clean :
	rm -rf *.o *.mod *.MOD *~

fit_pixcmd : fit_pixcmd.o $(COMMON)
	$(F90) -o fit_pixcmd.exe fit_pixcmd.o $(COMMON)

dohess : dohess.o $(COMMON)
	$(F90) -o dohess.exe dohess.o $(COMMON)

write_a_model : write_a_model.o $(COMMON)
	$(F90) -o write_a_model.exe write_a_model.o $(COMMON)

#------------------------------------------------------------------------#

fit_pixcmd.o : nrtype.o pixcmd_vars.o nr.o pixcmd_utils.o ran_state.o

dohess.o : nrtype.o pixcmd_vars.o nr.o pixcmd_utils.o 

write_a_model.o : nrtype.o pixcmd_vars.o nr.o pixcmd_utils.o 

ran_state.o: nrutil.o

%.o : %.f90
	$(F90) $(F90FLAGS) -o  $@ -c $<

%.o : nr/%.f90
	$(F90) $(F90FLAGS) -o $@ -c $<


