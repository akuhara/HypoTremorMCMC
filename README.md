# HypoTremorMCMC (Beta version)

__Locationg Tectonic Tremors with Uncertainty Estimates__

Copyright (C) 2023 __Takeshi Akuhara__[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-6129-8459)

---

This repository contains software written in modern fortran (F2008 style) for 
determining the source locations tremor signals without clear phase onsets. The details of the method is provided in a preprint 
[(Akuhara et al. 2023, EarthArXiv)](https://doi.org/10.31223/X59S9J).


 
## __NOTE: This software is currently under development and the documentation may not be complete at this time. Any changes to the sorce code or documentation may occur without prior notification.__

## Installation 

To install the software, you need to have a Fortran compiler, OpenMPI, and FFTW installed on your system. 


1. Clone the repository to your local machine:
```
git clone https://github.com/username/HypoTremorMCMC.git` 
```

2. Comple the program by running the following command in the `src` directory of the repository: 

```
cd dir
make
``` 
## How It Works

Follow the five steps below to estimate tremor locations.

### Step 1: Envelope conversion

This step reads continuous seismic records of two horizontal components from SAC format files and convert them into a smoothed envelope. 

#### Usage
```
 mpirun -np [process number] hypo_tremor_convert [parameter file]
```

### Step 2: Calculate cross-correlation

This step calculates cross-correlation functions of the envelope between all station pairs. 

#### Usage
```
 mpirun -np [process number] hypo_tremor_correlate [parameter file]
```

### Step 3: Measure & optimize time- and amplitude- difference

This step first measure the arrival time- and amplitude-difference between station pairs and then optimize them to obtain staiton-specific relative measurements.

```
 mpirun -np [process number] hypo_tremor_measure [parameter file]
```

### Step 4: Select good-quality events

This step selects events with good-quality on the basis of the rough estimates of propagation speed and attenuation strengths of a seismic wave.

#### Usage
```
 mpirun -np [process number] hypo_tremor_select [parameter file]
```

### Step 5: Perform MCMC 

This step peforms Bayesian inversion using Markov-chain Monte Carlo (MCMC) method to obtain the posterior probability of source locations.

#### Usage
```
 mpirun -np [process number] hypo_tremor_mcmc [parameter file]
```


## Contributing

Contributions to the software are welcome. Please open an issue or submit a pull request on the GitHub repository.
