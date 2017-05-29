#FC=gfortran
#FFLAGS=-fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none -fopenmp
#FWARN= #-Wall -fbounds-check
#FOPT=-march=native -flto -funroll-loops -O3 #-fblas  (check this blas one to replace matmul calls)
#FPROF= #-g -pg

#Intel Compiler
FC=mpif90
FWARN= #-CB
FFLAGS=-r8 -fpp
FOPT=-O3 -xHost -ipo -openmp #-parallel
FPROF= #-pg -g

#FFTW_INC=-I/usr/local/fftw-3.3.3/include/
#FFTW_LIB=-L/usr/local/fftw-3.3.3/lib

#FFTW = -lfftw3_mpi -lfftw3
FFTW = -lfftw3_mpi -lfftw3_omp -lfftw3

#ifdef OMP
#FFLAGS+= -fopenmp
#FFTW+= -lfftw3_threads -lfftw3
#endif

#ifdef THREADS
#THREAD_LIB = -lpthread
#FFTW = -lfftw3_threads -lfftw3

#HAM = hamiltonian_conformal.o 
#HAM = hamiltonian_cosmic.o
#HAM = hamiltonian_minkowski.o
#HAM = hamiltonian_conformal_noncanonical_rksplit.o
#HAM = hamiltonian_conformal_nonminimal_gravity.o
HAM = hamiltonian_noncanonical_conformal_flatten.o
OBJS = fftw_mod_mpi.o params.o analysis_spec.o $(HAM) initialize.o integrator.o

spectral: %: $(OBJS) evolve_spectral.o
	$(FC) $(FFTW_INC) $(FFLAGS) $(FOPT) $(FPROF) $(FWARN) -o runspec evolve_spectral.f90 $(OBJS) $(FFTW_LIB) $(FFTW) -lm $(THREAD_LIB)

%.o: %.f90
	$(FC) $(FFLAGS) $(FFTW_INC) $(FOPT) $(FPROF) $(FWARN) -c $< -o $@ $(FFTW_LIB) $(FFTW) -lm $(THREAD_LIB)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f runspec