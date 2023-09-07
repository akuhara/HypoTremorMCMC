# HypoTremorMCMC (Beta version)

__Locationg Tectonic Tremors with Uncertainty Estimates__

Copyright (C) 2023 __Takeshi Akuhara__[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-6129-8459)

---

This repository contains software written in modern fortran (F2008 style) for 
determining the source locations tremor signals without clear phase onsets. The details of the method is provided in a preprint 
[(Akuhara et al. 2023, EarthArXiv)](https://doi.org/10.31223/X59S9J).


 
## __NOTE: This software is currently under development and the documentation may not be complete at this time. Any changes to the sorce code or documentation may occur without prior notification.__

## Installation 

To install the software, you need to have a Fortran compiler, OpenMPI, and FFTW installed on your system.  To get started, clone the repository to your local machine by running:

```
git clone https://github.com/akuhara/HypoTremorMCMC.git 
```

Then, compile the program by running the following command in the src directory of the repository:

```
cd dir
make
``` 

It is important to note that you may need to edit the Makefile in accordance with your environment.

## Usage

To estimate tremor locations, follow these five steps:

1. __Envelope Conversion:__ Read continuous seismic records of two horizontal components from SAC format files and convert them into a smoothed envelope.

```
 mpirun -np [process number] hypo_tremor_convert [parameter file]
```

2. __Calculate Cross-correlation:__ Calculate cross-correlation functions of the envelope between all station pairs.

```
 mpirun -np [process number] hypo_tremor_correlate [parameter file]
```

3. __Measure and Optimize Time- and Amplitude-difference:__ Measure the arrival time- and amplitude-difference between station pairs and then optimize them to obtain station-specific relative measurements.

```
 mpirun -np [process number] hypo_tremor_measure [parameter file]
```

4. __Select Good-quality Events:__ Select events with good-quality based on rough estimates of propagation speed and attenuation strengths of a seismic wave.

```
 mpirun -np [process number] hypo_tremor_select [parameter file]
```

5. __Perform MCMC:__ Perform Bayesian inversion using Markov-chain Monte Carlo (MCMC) method to obtain the posterior probability of source locations.

```
 mpirun -np [process number] hypo_tremor_mcmc [parameter file]
```

6. __Do statistics:__ Extract statistical information from MCMC samples.

```
 mpirun -np [process number] hypo_tremor_statistics [parameter file]
```

You can find detailed instructions on how to use this software by visiting the [online documentation](https://hypotremormcmc.readthedocs.io/en/latest/).

## Contributing

Contributions to the software are welcome. Please [open an issue](https://github.com/akuhara/HypoTremorMCMC/issues) or submit a pull request on the [GitHub repository](https://github.com/akuhara/HypoTremorMCMC).
