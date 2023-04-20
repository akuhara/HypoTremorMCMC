# HypoTremorMCMC

__Locationg Tectonic Tremors with Uncertainty Estimates__

Copyright (C) 2023 __Takeshi Akuhara__[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-6129-8459)

---

This repository contains software written in modern fortran (F2008 style) for 
determining the source locations tremor signals without clear phase onsets. The details of the method is provided in a preprint 
[(Akuhara et al. 2023, EarthArXiv)](https://doi.org/10.31223/X59S9J).
 


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

### Usage

To use the program, follow these steps: 

1. Create a parameter file with the necessary input parameters, and other files containing required information, such as input filename convention and station arrangement. You also need to prepare data file in SAC format. See the example file `hypo_tremor_mcmc.in` for guidance on the required format. 


2. Run the five programs one by one:

```
mpirun -np [process number] hypo_tremor_convert   [parameter file]
mpirun -np [process number] hypo_tremor_correlate [parameter file]
mpirun -np [process number] hypo_tremor_measure   [parameter file]
mpirun -np [process number] hypo_tremor_select    [parameter file]
mpirun -np [process number] hypo_tremor_mcmc      [parameter file]
```

The above five programs (`hypo_tremor_convert`, `hypo_tremor_correlate`, `hypo_tremor_measure`, `hypo_tremor_select`, and `hypo_tremor_mcmc`) must be executed in this order because some necessary input files are created by the previous program. 


## Contributing

Contributions to the software are welcome. Please open an issue or submit a pull request on the GitHub repository.
