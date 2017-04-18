# Problem description

The acoustic Full Waveform Inversion (FWI) aims to generate high resolution subsoil velocity models from acquired seismic data through an iterative process. Time-dependent seismic wave equation is solved forward and backward in time in order to estimate the velocity model. From the differences between acquired data and the computed velocity model, a gradient volume is calculated and used to update the velocity model on the next iteration.

The inverse problem is not linear and ill-conditioned. This makes solving the problem at high frequencies difficult. Instead, the initial stimulus is decomposed into a spectrum of frequencies. Then, low frequencies are solved first on a coarse grid providing a good guess for higher frequencies.
Conceptually, FWI can be divided into three main steps. A pre-processing step estimates the computational resources needed to solve the problem according to the number of shots, wavelet frequency and domain dimensions. Then, the wave propagator solves the timedomain formulation of the wave equation forward and backward in time. Finally, a postprocessing step gathers the information from the computation of all different frequencies into a single final velocity model.

# Compile and run on different architectures/machines

This application uses CMake to discover all dependences and build the application. Each directory contains its own `CMakeLists.txt` script. The build process follows the typical build process of every cmake application.

```bash
cmake -DCMAKE_C_COMPILER=<foo-compiler> [ -D<OPTION_1>=<yes|no> -D<OPTION_2>=<YES|NO> ... ]  <path-to-project-base-dir>
```

__WARNING:__ *Always* make *out-of-source* builds (don't execute cmake from the project root directory):
```bash
cd FWIDIR/<repo-name>/
mkdir build && cd build
cmake <options> ..
make
```
We provide different scripts to set up a specific environment on different machines. These scripts are localted under `FWIDIR/Scripts` directory. To set up an environment, simply source it:

```bash
source FWIDIR/Scripts/<environment_file.sh>
```
# Compilers and versions
    - Intel Compiler (>= 15)
    - GCC Compiler (6.2.0)

NOTE: 
1. To execute Offload in BSC machines:
   - module load ompss/stable
   - module load intel/16.0.0
   - module load gcc/6.2.0
2. To run on BSC's KNLs
   - Use the slurm.sh file to reserve a timeslot (e.g 1 hour)
   - Check if the job has been assigned to a node (use squeue command)
   - Ssh to the assigned node (ssh knl0x, where x can be a number between 1 and 4)
   - Then it is enable to run the mockup
