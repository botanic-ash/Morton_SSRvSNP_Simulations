# Overview

This repository contains the code required for the Simulation component of the IMLS_GCCO marker comparison project, which is quantifying how measures of *ex situ* conservation differ between RADseq and microsatellite markers.
Simulations using different marker types (microsatellites, or "MSAT", and SNP markers, or "DNA") 
are performed using the [fastSimcoal2 software](http://cmpg.unibe.ch/software/fastsimcoal27/) and the [strataG](https://github.com/EricArcher/strataG) R package.

### RScripts
This folder contains the R scripts use to specify simulation parameters, run simulations, convert Arlequin outputs to genind objects,
and process those genind objects.

Functions used repeatedly throughout the scripts are declared in [SSRvSNP_Sim_functions.R](https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/SSRvSNP_Sim_functions.R).

The script [GenerateFSCparams.R](https://github.com/akoontz11/Morton_SSRvSNP_Simulations/blob/main/RScripts/GenerateFSCparams.R) specifies the parameters used for
MSAT and DNA fastSimcoal2 simulations, and runs those simulations, storing the results to .params R objects. Note that every time this script is 
run or sourced, new fastSimcoal2 simulation output files will populate the specified [SimulationOutputs](https://github.com/akoontz11/Morton_SSRvSNP_Simulations/tree/main/SimulationOutputs) directory (see below).

### SimulationOutputs
This folder contains 2 subfolders, one for simulation outputs using microsatellite markers ([MSAT_marker](https://github.com/akoontz11/Morton_SSRvSNP_Simulations/tree/main/SimulationOutputs/MSAT_marker)) 
and one for simulation outputs using SNP markers ([DNA_marker](https://github.com/akoontz11/Morton_SSRvSNP_Simulations/tree/main/SimulationOutputs/DNA_marker)). 

Within each of these are folders for each simulation scenario containing Arlequin output files. 
Additionally, simulation parameter files, and the log files for simulations, are contained in the top directory of each marker type subfolder

## Contact
For questions about this project, open an Issue or contact Austin Koontz.
