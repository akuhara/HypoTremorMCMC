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

## Input file configuration

To structure your input SAC data files as shown in the figure below, you can update the [parameter file](./parameter_file.md) using the following settings:
* Set the _filename_format_ parameter to "$STA + / + $ID + . + $CMP"
* Set the _data_dir_ parameter to "SAC_data"
* Specify _cmp1_ as "BHE" and _cmp2_ as "BHN"
The path to the parameter file must be provided as the first command line argument. In the case of the figure below, where "hypo_tremor.in" is the parameter file, the program should be executed as:
 
`mpirun -np 20 hypo_tremor_convert hypo_tremor.in`.

![files](./img/files_convert.png)

In addition, you will need to provide a station file ("station.txt" in the figure above) and specify its path in the parameter file as _station_file_. The station file contains information on the location of the stations in the X-Y-Z coordinate system and the sensitivity of the sensors. For detailed instructions on how to create this file, please refer to [this guide (station file)](statiln_file.md).

You will also need a file ("time_id.txt" in the figure above) that lists the time identifiers (ID) for each of your SAC data files. The time ID typically follows the format of yymmdd.HHMMSS, where yy represents the year, mm represents the month, dd represents the day, HH represents the hour, MM represents the minute, and SS represents the second. However, it is not necessary for the IDs to explicitly indicate time information. The program will read the SAC files that have the ID specified in the list, in the order specified by the list, and concatenate them to create a single time series. Specifically, it uses the _filename_format_ parameter specified in the parameter file to find the corresponding SAC files.

## Output envelope

After running the program successfully, you will have `<station name>.merged.env` files in the current directory. These files have a binary format and contain a series of time-amplitude doublets representing a smoothed envelope. To read and visualize these files, the simplest way is to use `gnuplot` as shown below:

```
gnuplot
plot "<station name>.merged.env" binary form=%double%double" rec=(-1) endian=little with lines
```
Note that the value for `endian` can be `big`, depending on your environment. 




