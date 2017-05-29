################################################################
# $Id: Makefile,v 1.2 2007/09/19 23:14:02 frolov Exp frolov $
# DEFROST Makefile (Intel Fortran compiler)
################################################################

.SUFFIXES: .f .F .f90

ifndef MPIFC
MPIFC = mpif90
endif

PG_INFO=$(shell $(MPIFC) -V 2>&1 | grep pgf90)
ifneq ($(PG_INFO),)
USE_PG=1
endif

GCC_INFO=$(shell $(MPIFC) -v 2>&1 | grep "gcc version")
ifneq ($(GCC_INFO),)
USE_GNU=1
endif

FFLAGS_DEFS =
BASEINC = -I$(shell cd $(dir $(shell which fftw-wisdom))../include && pwd)

BASELIBDIR = $(shell cd $(dir $(shell which fftw-wisdom))../lib && pwd)
BASELIB = -L$(BASELIBDIR)
LDFLAGS = $(BASELIB)
export LD_RUN_PATH = $(BASELIBDIR)

ifdef USE_PG
FFLAGS = -O3 -fast -r8 # -Mconcur -mp
FFLAGS_COMMON = -Mpreprocess
FFLAGS_FREE = -Mfreeform -Mextend
FFLAGS_FIXED = -Mfixed -Mextend
LDFLAGS += 
else
ifdef USE_GNU
FFLAGS = -O3 -march=pentium4 -fdefault-real-8 -fdollar-ok # -fopenmp
FFLAGS_COMMON = -x f95-cpp-input 
FFLAGS_FREE = -ffree-line-length-none
FFLAGS_FIXED = -ffixed-form -ffixed-line-length-none
LDFLAGS += 
else
# Fortran compiler (adjust for your machine, -r8 and -fpp are mandatory)
#FFLAGS = -O3 -ipo -xP -r8 -pc80 # -parallel -openmp
FFLAGS = -r8
FFLAGS_COMMON = -fpp
FFLAGS_FREE = 
FFLAGS_FIXED = -132
#FFLAGS_FIXED = 
LDFLAGS += -static-intel
endif
endif

FFTWINC = 
FFTWLIB = -lfftw3

NO_FFTW_THREADS=1
ifneq ($(shell ls $(BASELIBDIR) | grep libfftw3_threads),)
ifndef NO_FFTW_THREADS
FFTWLIB += -lfftw3_threads
FFLAGS_DEFS += -DUSE_FFTW_THREADS
endif
endif

# Don't need these with precompiled p3dfft
FFT3DLIB = p3dfft/build/module.o p3dfft/build/fft_spec.o p3dfft/build/fft_init.o p3dfft/build/fft_exec.o p3dfft/build/wrap.o
FFLAGS_DEFS += -DFFTW -DDOUBLE_PREC

#FFT3DLIB = -L/home/jbraden/etc/ -lp3dfft
SPHERELIB = -L/home/jbraden/etc/SHTOOLS/lib -lSHTOOLS2.6
SPHEREINC = -I/home/jbraden/etc/SHTOOLS/modules/

# SILO library (native VisIt data format)
ifdef USE_SILO
SILOINC = 
SILOLIB = -L/usr/local/lib64 -lsilo -lhdf5

FFLAGS_DEFS += -DSILO
endif

FIXKIC = 
FIXKIC = -DDONT_FIX_ICS

KFAC = 1.0  # 2.0 is defrost standard

################################################################

L_VALUES = 10 # 100 200 400
N_VALUES = 256 # 64 128 256 512 1024
SEED_VALUES = 2 # 1 2 3 4 5
#EXES = $(foreach L,$(L_VALUES),$(foreach N,$(N_VALUES),$(foreach SEED,$(SEED_VALUES),defrost_$L_$N_$(SEED))))
#OBJS = $(foreach E,$(EXES),$E.o)

OBJ = modparam.o params.o latpar.o model.o init_cond.o analysis.o scafld.o homlat.o outputfile.o iofile.o evolve.o

all: defrost

clean:
	rm -f *.bak gmon.out core *.o *.mod $(EXES)
	rm -f p3dfft/build/*.o
	rm -f defrost

cleandat:
	rm -f LOG*
	rm -f PSD*
	rm -f CDF*
	rm -f PDF*
	rm -f RG_*
	rm -f T_SERIES

defrost: %: %.o $(FFT3DLIB) $(OBJ)
	$(MPIFC) $(FFLAGS) $(FFLAGS_DEFS) $^ -o defrost $(LDFLAGS) $(FFTWLIB) $(SILOLIB) $(SPHERELIB)

################################################################

# Move this into a library archive that I store somewhere
p3dfft/build/module.o: p3dfft/build/module.F p3dfft/build/fft_spec.o
p3dfft/build/fft_init.o: p3dfft/build/fft_init.F p3dfft/build/module.o p3dfft/build/fft_spec.o
p3dfft/build/fft_exec.o: p3dfft/build/fft_exec.F p3dfft/build/module.o p3dfft/build/fft_spec.o
p3dfft/build/wrap.o: p3dfft/build/wrap.f p3dfft/build/module.o
defrost.o: p3dfft/build/module.o $(OBJ)

modparam.o: modparam.f90

latpar.o: p3dfft/build/module.o params.f90 latpar.f90    # check if I actually need this

params.o: modparam.f90 params.f90 

model.o: modparam.f90 params.f90 latpar.f90 model.f90

init_cond.o: modparam.f90 params.f90 latpar.f90 model.f90 init_cond.f90

analysis.o: modparam.f90 params.f90 model.f90 analysis.f90

scafld.o: p3dfft/build/module.o modparam.f90 params.f90 model.f90 latpar.f90 scafld.f90

homlat.o: modparam.f90 params.f90 model.f90 latpar.f90 scafld.f90 homlat.f90

outputfile.o: modparam.f90 params.f90 model.f90 analysis.f90 outputfile.f90

iofile.o: modparam.f90 params.f90 latpar.f90 model.f90 scafld.f90 homlat.f90 analysis.f90 iofile.f90

evolve.o: modparam.f90 params.f90 latpar.f90 scafld.f90 homlat.f90 evolve.f90

# need to put this in with defrost
#driver.o: modparam.f90 params.f90 model.f90 latpar.f90 fields.f90 homlat.f90 evolve.f90 outputfile.f90 iofile.f90 driver.f90

defrost.o: driver.f90
	$(MPIFC) $(BASEINC) $(SPHEREINC) -c $(FFLAGS_COMMON) $(FFLAGS_FREE) $(FFLAGS) $(FFLAGS_DEFS) -DSEED_VALUES=2 -DKCUT_FAC=$(KFAC) $(FFTWINC) $(SILOINC) $(SPHERELIB) $< -o $@
#defrost.o: driver.f90
#	$(MPIFC) $(BASEINC) -c $(FFLAGS_COMMON) $(FFLAGS_FREE) $(FFLAGS) $(FFLAGS_DEFS) -DL_VALUE=$(word 2,$(subst _, ,$@)).0 -DN_VALUE=$(word 3,$(subst _, ,$@)) -DSEED=$(word 4,$(subst _, ,$(subst ., ,$@))) -DKCUT_FAC=$(KFAC) $(FIXKIC) $(FFTWINC) $(SILOINC) $< -o $@

%.o: %.f90
	$(MPIFC) $(SPHEREINC) -fpp2 -r8 -c $(FFTWINC) $(FFLAGS_DEFS) $< -o $@
#	$(MPIFC) -vec- -fpp2 -r8 -ipo -pc80 -xP -O3 -c $(FFTWINC) $< -o $@

# Not sure why I bothered to separate these, again only need for the p3dfft lirary, so I should just make a library
%.o: %.F
	$(MPIFC) $(BASEINC) -c $(FFLAGS_COMMON) $(FFLAGS_FIXED) $(FFLAGS) $(FFLAGS_DEFS) $(FFTWINC) $< -o $@

%.o: %.f
	$(MPIFC) $(BASEINC) -c $(FFLAGS_COMMON) $(FFLAGS_FIXED) $(FFLAGS) $(FFLAGS_DEFS) $(FFTWINC) $< -o $@

