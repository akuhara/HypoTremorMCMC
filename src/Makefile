#-----------------------------------------------------------------------
# Complier  
#-----------------------------------------------------------------------

# Example 1. GNU fortran for debug
FC     = mpifort
FFLAGS = -g -Wall -pedantic -std=f2008 -fbounds-check -O0 \
         -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	 -fconvert=big-endian

# Example 2. GNU fortran for fast computation
#FC     = mpifort
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check

# Example 3. Intel compiler for debug
#FC      = mpif90
#FFLAGS  = -assume byterecl -traceback -g -check bounds -fpe0 -lmpi 
#MKLROOT = /home/app/intel/compilers_and_libraries_2019.5.281/linux/mkl


#-----------------------------------------------------------------------
# Libararies
#-----------------------------------------------------------------------

# Example A. General case
#FFTW   = -I/usr/local/include -lfftw3
#LAPACK = -llapack -lblas

# Example B. EIC
#FFTW    = -I$(MKLROOT)/include/fftw -mkl=parallel
#LAPACK  =  # leave this empty


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------
COV =
BINDIR = ../bin

HYPO_TREMOR_MCMC = $(BINDIR)/hypo_tremor_mcmc
OBJS_HYPO_TREMOR_MCMC = \
		hypo_tremor_mcmc.o cls_param.o cls_line_text.o \
	 	cls_c3_data.o cls_convertor.o


TARGET = $(HYPO_TREMOR_MCMC)

.PHONY: all
all: $(TARGET)

hypo_tremor_mcmc.o: cls_param.mod cls_line_text.mod cls_convertor.mod
cls_param.o: cls_line_text.mod
cls_convertor.o: cls_c3_data.mod cls_param.mod cls_line_text.mod


$(HYPO_TREMOR_MCMC): $(OBJS_HYPO_TREMOR_MCMC)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $^ -o $@ 


.PHONY: clean
clean:
	rm -f *.mod $(BINDIR)/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) $(COV) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
