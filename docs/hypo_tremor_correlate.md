# Step 2. Calculate cross-correlations

The `hypo_tremor_correlate` program reads the envelope file `<station name>.merged.env` generated `hypo_tremor_convert` and calculates cross-correlation functions. of all station pirs using the specified time window length (_t_win_corr_) and time interval (_t_step_corr_). The values for _t_win_corr_ and _t_step_corr_ must be specified in the parameter file. 

Two types of output files are generated after the successful running of this program: `<station name 1>.<station name 2>.corr`, and `<station name 1>.<station name 2>.max_corr`. The former contains "time, lag-time, & correlation value" triplets in a binary format. The latter contains "time & maximum correlation value" doublets in a binary format.

## Required Parameters

The following parameters must appear in the [parameter file](parameter_file.md). 

__n_procs__

Number of processes for MPI parallel computation. This value must be equal to the one you specify in the command line. The purpose of this parameter is to prevent unintended number specified in the command line, which for some programs results in unintended results.

__station_file__

Absolute or relative path to a [station file](station_file.md) that lists station names, their position (X-Y-Z), and sensitivities of two horizontal components.

__t_win_corr__

Length of time window for cross-correlation [sec].

__t_step_corr__

Interval of time windows [sec].
