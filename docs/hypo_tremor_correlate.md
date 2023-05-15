# Step 2. Calculate cross-correlations

The `hypo_tremor_correlate` program reads the envelope file `<station name>.merged.env` generated `hypo_tremor_convert` and calculates cross-correlation functions. of all station pirs using the specified time window length (_t_win_corr_) and time interval (_t_step_corr_). The values for _t_win_corr_ and _t_step_corr_ must be specified in the parameter file. 

Two types of output files are generated after the successful running of this program: `<station name 1>.<station name 2>.corr`, and `<station name 1>.<station name 2>.max_corr`. 