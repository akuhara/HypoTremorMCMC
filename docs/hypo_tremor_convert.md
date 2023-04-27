# Step 1. Envelope conversion

To locate tremor sources, the first step is to generate envelope waveforms. Tremor signals are predominantly composed of horizontally polarized S-wave energy, which requires two horizontal components to undergo envelope conversion. The `hypo_tremor_convert` program performs this process by reading continuous waveform data from SAC files and applying several signal processing techniques such as detrending, tapering, band-pass filtering, envelope conversion, smoothing, merging the two components, and decimation.

The `hypo_tremor_convert` program  utilizes the Open MPI framework for parallelization, enabling efficient processing of multiple stations' data by assigning each process with multiple stations to handle. To execute the program, a parameter file must be specified as the first command line argument, containing crucial information about SAC data file configurations, station information, and other relevant details. The parameter file format can be found in this resource (parameter file) for further reference.

The conversion to envelopes typically involves processing a vast amount of time series data spanning months to years, making it impractical to handle all at once. Usually, the time series is divided into smaller fragments capable of signal processing with fast Fourier transform. However, this segmentation poses a challenge as the endpoints of each fragment are subject to tapering, leading to information loss. The `hypo_tremor_convert` program addresses this issue by processing segments with half-overlapping adjacent segments. From the resultant fragmented envelopes, the program extracts the portions not affected by the tapering and connects them to create a smooth envelope across the entire observation period. 

## Required Parameters

The following parameters must appear in the parameter file.

### station_file 
Absolute or relative path to a file that lists station names, their position (X-Y-Z), and sensitivities of two horizontal components.

### data_dir
Absolute or relative path to a directory that contains SAC data files.

### time_id_file
Absolute or relative path to a file containing a list of time identifiers in chronological order, which appear in the filenames of SAC data files. 

### cmp1
Name of the first horizontal component specified in the filename of SAC data files.

### cmp2
Name of the second horizontal component specified in the filename of SAC data files.

### filename_format 
Filename convention for SAC data files constituted of the three reserved variables ($STA, $ID, and $CMP, representing station name, time ID, and component, respectively), the concatenation operator (+), and arbitral character strings.

### t_win_conv 
Length of a time window to which the envelope transformation is applied [sec]

## Tips

![files](../img/files_convert.png)


