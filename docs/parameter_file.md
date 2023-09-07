# Parameter file

All programs in this package take a parameter file as its first command line argument, which allows users to set tuning parameters or specify paths to SAC data files and station lists according to their preferences. Each line in the parameter file consists of a parameter name and its value separated by the "=" symbol. Users can add comments after the "#" symbol, just like in many programming languages. Additionally, spaces are automatically ignored, so "station_list = hoge" and "station_list = h o ge" are equivalent, although we recognize that this may seem unusual to most users. 

The parameter file only accepts pre-defined parameter names that are listed in the templete below, and the program will terminate with an error message if any unexpected parameter names are found in the file. If the required parameter names are missing from the file, the program will also terminate with an error message.


__Templete of the parameter file__
```
#------------------------------------------#
# Sample Parameter File for HypoTremorMCMC #
#------------------------------------------#

#------------------#
# For all programs #
#------------------#
#
# n_procs      : Number of processes for MPI parallel computation.
#                This value must be equal to the one you specify in
#                the command line. The purpose of this parameter is to
#                prevent unintended number specified in the command line,
#                which for some programs results in unintended results. 
#
# station_file : Input file that lists station names, their position (X-Y-Z), 
#                and sensitivity
#

n_procs = 20

station_file = station_xy.list

#-------------------------#
# For hypo_tremor_convert #
#-------------------------#
#
# data_dir        : Directory containing SAC data files
#
#
# time_id_file    : Name of the input file containing a list of time identifiers 
#                   in chronological order, which appear in the file names of SAC
#                   data files
#
# cmp1            : Name of the first horizontal component specified in the
#                   filename of SAC data files
#
# cmp2            : Name of the second horizontal component specified in the
#                   filename of SAC data files
#
# filename_format : Filename convention for SAC data files constituted of
#                   the three variables ($STA, $ID, and $CMP, representing, 
#                   station name, time ID, and component, 
#                   respectively), the concatenation operator (+), and 
#                   arbitral character strings
#
# t_win_conv      : Length of a time window to which the envelope 
#                   transformation is applied [sec]
#

data_dir        = ./data

time_id_file    = dataID_2d.list

cmp1            = EH1

cmp2            = EH2

filename_format = $STA + / + $ID + . + $CMP 

t_win_conv      = 3000.0                    

#---------------------------#
# For hypo_tremor_correlate #
#---------------------------#
#
# t_win_corr  : Length of a time window for which cross-correlation is  
#               is calculated [sec]
#
#
# t_step_corr : Interval of time windows [sec]
#

t_win_corr      = 300.0

t_step_corr     = 150.0 

#-------------------------#
# For hypo_tremor_measure #
#-------------------------#
#
# alpha        : Percentile on the cross-correlation histogram, specifying
#                the cross-correlation value required for tremor detection
#                (0.0 <= alpha < 1.0)
#
# n_pair_thred : The number of station pairs that meet the cross-correlation 
#                requirement for tremor detection
#

alpha           = 0.98

n_pair_thred    = 300 

#------------------------#
# For hypo_tremor_select #
#------------------------#
#
# z_guess : Assumed depth for tremor sources [km]
# 
# vs_min  : Minimum S-wave velocity to be accepted [km/s]
#
# vs_max  : Maximum S-wave velocity to be accepted [km/s]
#
# b_min   : Minimum attenuation strength to be accepted
#
# b_max   : Maximum attenuation strength to be accepted
# 

z_guess = 7.0

vs_min = 2.0

vs_max = 4.0

b_min  = 0.015

b_max  = 0.03

#----------------------#
# For hypo_tremor_mcmc #
#----------------------#
#
# < Group 1 - Iterations >
#
# n_iter             : Number of iterations, including a burn-in period
# 		          
# n_burn             : Number of iterations within the burn-in period
# 		          
# n_interval         : Interval of iterations at which model parameters are saved
# 		          
# n_chains           : Number of MCMC chains per process
# 		          
# n_cool             : Number of non-tempered MCMC chains per process
# 		          
# temp_high          : Heighest temperature for parallel tempering scheme
#

n_iter = 4000000

n_burn = 2000000

n_interval = 1000

n_chains = 5

n_cool   = 1

temp_high = 200.d0

#
# < Group 2 - Prior probability > 
#        
# prior_z            : Prior information on the shallowest limit of event depth
#                      [km]
#			     
# prior_width_z      : Parameter controlling prior width of event depth [km]
#                      The prior probability has the maximum amplitude at 
#                      <prior_width_z> + <prior_z> km
#		       		            
# prior_width_xy     : Parameter controlling prior width of event horizontal 
#                      locations (i.e., standard deviation of Gaussian 
#                      distribution) [km]
#                      
# prior_vs           : Parameter controlling prior mean of S-wave velocity [km/s] 
# 		          
# prior_width_vs     : Parameter controlling prior width of S-wave velocity 
#                      (i.e., standard deviation of Gaussian distribution)
#                      [km/s]
#
# prior_qs           : Parameter controlling prior mean of Qs
#
# prior_width_qs     : Parameter controlling prior width of Qs 
#                      i.e., standard deviation of Gaussian distribution)
#
# prior_width_t_corr : Parameter controlling prior width of arrival time 
#                      correction (i.e., standard deviation of Guassian
#                      distribution) [sec]
#
# prior_width_a_corr : Parameter controlling prior width of amplitude correction
#                      (i.e., standard deviation of Guassian distribution)
#


prior_z            = 0.0

prior_width_z      = 10.0

prior_width_xy     = 30.0

prior_vs           = 3.0

prior_width_vs     = 1.0

prior_qs           = 250

prior_width_qs     = 100

prior_width_t_corr = 0.5

prior_width_a_corr = 0.02

#
# < Group 3 - Random walk >
#
# step_size_z      : Parameter controlling random walk's step size for event depth
#                    (i.e., standard deviation of Gaussian distribution) [km]
# 
# step_size_xy     : Parameter controlling random walk's step size for event 
#                    horizontal locations (i.e., standard deviation of Gaussian 
#                    distribution) [km]
# 
# step_size_vs      : Parameter controlling random walk's step size for S-wave
#                     velocity (i.e., standard deviation of Gaussian 
#                     distribution) [km/s]
# 
# step_size_qs      : Parameter controlling random walk's step size for Qs
#                    (i.e., standard deviation of Gaussian distribution) [km]
# 
# 
# step_size_t_corr  : Parameter controlling random walk's step size for arrival
#                     time correction (i.e., standard deviation of Gaussian 
#                     distribution) [sec]
# 
# step_size_a_corr  : Parameter controlling random walk's step size for amplitude
#                     correction(i.e., standard deviation of Gaussian 
#                     distribution)
# 

step_size_z = 0.4

step_size_xy = 2.0

step_size_vs = 0.2

step_size_qs = 5.0

step_size_t_corr = 0.03

step_size_a_corr = 0.005

#
# < Group 4 - Others >
#
# solve_vs       : Whether to solve S-wave velocity [T/F]
#
# solve_t_corr   : Whether to solve arrival time corrections [T/F]
#
# solve_qs       : Whether to solve Qs [T/F]
#
# solve_a_corr   : Whether to solve amplitude corrections [T/F]
#
# use_time        : Whether to use arrival time data [T/F]
#
# use_amp         : Whether to use amplitude data [T/F]
#

solve_vs = T    

solve_t_corr = T

solve_qs = T

solve_a_corr = T

use_time = T

use_amp = T


```

