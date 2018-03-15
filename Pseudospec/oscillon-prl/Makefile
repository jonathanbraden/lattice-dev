#FC=gfortran
#FFLAGS=-fdefault-real-8 -fdefault-double-8 -cpp -ffree-line-length-none -fopenmp
#FWARN= #-Wall -fbounds-check
#FOPT=-march=native -flto -funroll-loops -O3 #-fblas  (check this blas one to replace matmul calls)
#FPROF= #-g -pg

#Intel Compiler
FC=ifort
FFLAGS=-r8 -fpp
FOPT= -O3 -ipo -openmp #-xHost -parallel
#FPROF= #-pg -g

#FFTW_INC=-I/usr/local/fftw-3.3.3/include/
#FFTW_LIB=-L/usr/local/fftw-3.3.3/lib
FFTW_INC=-I/share/apps/fftw-3.3.5/include/
FFTW_LIB=-L/share/apps/fftw-3.3.5/lib

FFTW = -lfftw3_omp -lfftw3

#ifdef OMP
#FFLAGS+= -fopenmp
#FFTW+= -lfftw3_threads -lfftw3
#endif

#ifdef THREADS
#THREAD_LIB = -lpthread
#FFTW = -lfftw3_threads -lfftw3

HAM = hamiltonian_conformal.o #hamiltonian_minkowski.o hamiltonian_conformal.o 
OBJS = fftw_mod.o params.o analysis_spec.o $(HAM)  #integrator.o

spectral: %: $(OBJS) evolve_spectral.o
	$(FC) $(FFTW_INC) $(FFLAGS) $(FOPT) $(FPROF) $(FWARN) -o runspec evolve_spectral.f90 $(OBJS) $(FFTW_LIB) $(FFTW) -lm $(THREAD_LIB)

%.o: %.f90
	$(FC) $(FFLAGS) $(FFTW_INC) $(FOPT) $(FPROF) $(FWARN) -c $< -o $@ $(FFTW_LIB) $(FFTW) -lm $(THREAD_LIB)

clean:
	rm -f *.o
	rm -f *.mod
	rm -f runspec