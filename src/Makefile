#-----------------------------------------------------------------------
# Complier  
#-----------------------------------------------------------------------

# Example 1. GNU fortran for debug
FC     = mpif90
FFLAGS = -pg -Wall -pedantic -fbounds-check -O0 \
         -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	 -fconvert=big-endian

# Example 2. GNU fortran for fast computation
#FC     = mpifort
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check \
#         -fconvert=big-endian

# Example 3. Intel compiler for debug
#FC      = mpif90
#FFLAGS  = -assume byterecl -traceback -g -check bounds -fpe0 -lmpi 
#MKLROOT = /home/app/intel/compilers_and_libraries_2019.5.281/linux/mkl


#-----------------------------------------------------------------------
# Libararies
#-----------------------------------------------------------------------

# Example A. General case
FFTW   = -I/usr/local/include -lfftw3
LAPACK = -llapack -lblas

# Example B. EIC
#FFTW    = -I$(MKLROOT)/include/fftw -mkl=parallel
#LAPACK  =  # leave this empty


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------
COV =
BINDIR = ../bin
HYPO_TREMOR_DETECT = $(BINDIR)/hypo_tremor_detect
OBJS_HYPO_TREMOR_DETECT = mod_mpi.o \
		hypo_tremor_detect.o cls_param.o cls_line_text.o \
	 	cls_c3_data.o cls_convertor.o cls_detector.o \
                mod_signal_process.o mod_sort.o


TARGET = $(HYPO_TREMOR_DETECT)

.PHONY: all
all: $(TARGET)

hypo_tremor_detect.o: mod_mpi.mod cls_param.mod cls_line_text.mod cls_convertor.mod \
                    cls_detector.mod
cls_param.o: cls_line_text.mod
cls_convertor.o: cls_c3_data.mod cls_line_text.mod mod_signal_process.mod
cls_detector.o: cls_c3_data.mod cls_line_text.mod mod_signal_process.mod \
                mod_sort.mod


$(HYPO_TREMOR_DETECT): $(OBJS_HYPO_TREMOR_DETECT)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  


.PHONY: clean
clean:
	rm -f *.mod $(BINDIR)/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) $(COV) $(FFTW) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
