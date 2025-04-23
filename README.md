# genesis-nonbonded-kernels

## Description and programming model

This is kernel program of real-part nonbonding interaction forces. The program works with a single MPI with OpenMP parallerization.

The code is given from GENESIS v2.1.0. 

The original GENESIS is provided in the websit;

https://github.com/genesis-release-r-ccs/genesis

The modules for reading inputs are given from the kernel codes from Priority Issue Target Applications:

https://github.com/RIKEN-RCCS/fs2020-tapp-kernels

License : GNU Lessar General Public Licese version 3 or later.

## Kernel codes

- Generic: kernel givn from 'Generic' kernel of GENESIS. https://github.com/sugitalab/genesis-mkl-shared/blob/v2.1.0/src/spdyn/sp_energy_nonbond_intel.fpp

- Intel: kernel givn from 'Intel' kernel of GENESIS. https://github.com/sugitalab/genesis-mkl-shared/blob/v2.1.0/src/spdyn/sp_energy_nonbond_intel.fpp
- Intel_mod1: Modified 'Intel' kernel.

## Input data

The input files for the kernel program are located in the data/ directory. 
These files contain dumped memory data required to execute the kernel. 
Since the data size is very large, the files are compressed using bzip2. 
Please extract them before running the kernel.
The original system is composed of 92,224 atoms, taken from the GENESIS 
benchmark set (NPT simulation). 
The data was dumped from the 0-th MPI process out of a total of 16 processes.

 - data/data_kernel_generic.bz2 : input for Generic
 - data/data_kernel_Oct.bz2 : input for Intel/Intel_mod
