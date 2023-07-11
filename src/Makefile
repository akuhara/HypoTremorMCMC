#-----------------------------------------------------------------------
# Complier  
#-----------------------------------------------------------------------

# Example 1. GNU fortran for debug
FC     = mpif90
FFLAGS = -pg -Wall -pedantic -fbounds-check -O0 -std=f2008 \
         -Wuninitialized -ffpe-trap=invalid,zero,overflow -fbacktrace \
	 -fconvert=big-endian -g -fcheck=all

# Example 2. GNU fortran for fast computation
#FC     = mpif90
#FFLAGS = -ffast-math -march=native -mtune=native -O3 -fno-range-check \
         -fconvert=big-endian

# Example 3. Intel compiler for fast computation
#FC      = mpiifort
#FFLAGS  = -assume byterecl -lmpi -xAVX -O3 -convert big_endian
#MKLROOT = /home/app/intel/compilers_and_libraries_2019.5.281/linux/mkl


#-----------------------------------------------------------------------
# Libararies
#-----------------------------------------------------------------------

# Example A. General case
FFTW   = -I/usr/include -lfftw3 -L/usr/lib64

# Example B. EIC
#FFTW    = -I$(MKLROOT)/include/fftw -mkl=parallel


#-----------------------------------------------------------------------
# DONOT CHANGE BELOW
#-----------------------------------------------------------------------
COV =
BINDIR = ../bin
HYPO_TREMOR_CONVERT  = $(BINDIR)/hypo_tremor_convert
HYPO_TREMOR_CORRELATE   = $(BINDIR)/hypo_tremor_correlate
HYPO_TREMOR_MEASURE = $(BINDIR)/hypo_tremor_measure
HYPO_TREMOR_SELECT   = $(BINDIR)/hypo_tremor_select
HYPO_TREMOR_MCMC     = $(BINDIR)/hypo_tremor_mcmc
HYPO_TREMOR_STATISTICS     = $(BINDIR)/hypo_tremor_statistics

OBJS_HYPO_TREMOR_CONVERT = mod_mpi.o hypo_tremor_convert.o cls_param.o \
                           cls_c3_data.o cls_convertor.o mod_signal_process.o \
                           cls_line_text.o
OBJS_HYPO_TREMOR_CORRELATE = mod_mpi.o hypo_tremor_correlate.o cls_param.o \
                           cls_c3_data.o cls_correlator.o mod_signal_process.o \
                           cls_line_text.o mod_sort.o
OBJS_HYPO_TREMOR_MEASURE = mod_mpi.o \
		hypo_tremor_measure.o cls_param.o cls_line_text.o \
                cls_measurer.o mod_sort.o mod_signal_process.o
OBJS_HYPO_TREMOR_SELECT = mod_mpi.o \
                hypo_tremor_select.o cls_param.o cls_line_text.o cls_selector.o \
                mod_regress.o 
OBJS_HYPO_TREMOR_MCMC = mod_mpi.o \
                hypo_tremor_mcmc.o cls_param.o cls_line_text.o \
                cls_model.o mod_random.o cls_obs_data.o \
                cls_mcmc.o cls_parallel.o cls_forward.o

OBJS_HYPO_TREMOR_STATISTICS = mod_mpi.o \
                hypo_tremor_statistics.o cls_param.o cls_line_text.o \
                cls_statistics.o mod_sort.o

TARGET = $(HYPO_TREMOR_CONVERT) $(HYPO_TREMOR_MEASURE) $(HYPO_TREMOR_CORRELATE)\
         $(HYPO_TREMOR_SELECT) $(HYPO_TREMOR_MCMC) $(HYPO_TREMOR_STATISTICS)
.PHONY: all
all: $(TARGET)

hypo_tremor_convert.o: mod_mpi.mod cls_param.mod cls_line_text.mod cls_convertor.mod
hypo_tremor_correlate.o: mod_mpi.mod cls_param.mod cls_line_text.mod cls_correlator.mod \
                       cls_c3_data.mod
hypo_tremor_measure.o: mod_mpi.mod cls_param.mod cls_line_text.mod cls_measurer.mod
hypo_tremor_select.o: mod_mpi.mod cls_param.mod cls_line_text.mod  cls_selector.mod 
hypo_tremor_mcmc.o: mod_mpi.mod cls_model.mod cls_param.mod cls_line_text.mod \
                  cls_model.mod cls_obs_data.mod mod_random.mod cls_parallel.mod   \
	           cls_forward.mod
hypo_tremor_statistics.o: mod_mpi.mod cls_param.mod cls_line_text.mod cls_statistics.mod
cls_mcmc.o: mod_random.mod cls_model.mod
cls_parallel.o: mod_mpi.mod mod_random.mod cls_mcmc.mod 
cls_obs_data.o: cls_line_text.mod mod_random.mod
cls_model.o: mod_random.mod
cls_param.o: cls_line_text.mod
cls_convertor.o: cls_c3_data.mod cls_line_text.mod mod_signal_process.mod
cls_correlator.o: cls_c3_data.mod cls_line_text.mod mod_signal_process.mod \
                mod_sort.mod
cls_measurer.o: cls_line_text.mod mod_mpi.mod mod_signal_process.mod mod_sort.mod
cls_selector.o: cls_line_text.mod mod_regress.mod
cls_statistics.o: cls_line_text.mod mod_mpi.mod mod_sort.mod
cls_forward.o: cls_model.mod cls_obs_data.mod

$(HYPO_TREMOR_CONVERT): $(OBJS_HYPO_TREMOR_CONVERT)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  

$(HYPO_TREMOR_CORRELATE): $(OBJS_HYPO_TREMOR_CORRELATE)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  

$(HYPO_TREMOR_MEASURE): $(OBJS_HYPO_TREMOR_MEASURE)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  

$(HYPO_TREMOR_SELECT): $(OBJS_HYPO_TREMOR_SELECT)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  

$(HYPO_TREMOR_MCMC): $(OBJS_HYPO_TREMOR_MCMC)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  

$(HYPO_TREMOR_STATISTICS): $(OBJS_HYPO_TREMOR_STATISTICS)
	@if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR); fi
	$(FC) $(FFLAGS) $(COV) $(FFTW) $^ -o $@  


.PHONY: clean
clean:
	rm -f *.mod $(BINDIR)/* *.o

%.o: %.f90
	$(FC) $(FFLAGS) $(COV) $(FFTW) -c $< -o $*.o 
%.mod: %.f90 %.o
	@:
