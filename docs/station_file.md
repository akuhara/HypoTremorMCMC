# Station file

The station file is a critical component of the package, as it contains information about the stations used in the analysis, 
including their names, locations in the X-Y-Z coordinate system, and sensor sensitivities. 
The path to the station file must be specified in the parameter file using the _station_file_ parameter.

Each row of the station file contains information about a particular station, and is formatted with six columns. 
The first column contains the station name, while the second, third, and fourth columns contain the X, Y, and Z coordinates of the station, respectively. 
The Z coordinate is positive downward, in keeping with convention in seismology.

The final two columns of the station file contain the sensitivities of the first and second components of the station's horizontal sensors, respectively. 
These values represent the conversion factor from counts to meters per second for each component of the sensor. 
It is important to note that these values must be carefully determined for each station, as they can have a significant impact on the results of subsequent analyses.


| Column No.| 1           | 2       | 3       | 4                           | 5                                        | 6                                        |
| --------- | ----------- | ------- | ------- | --------------------------- | ---------------------------------------- | ---------------------------------------- | 
|           | station name| X in km | Y in km | Z in km (positive downward) | Sensitivity of 1st component (m/s/count) | Sensitivity of 2nd component (m/s/count) | 

