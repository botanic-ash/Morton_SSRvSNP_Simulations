# Overview

This repository is forked from Austin Koontz's Morton_SSRvSNP_Simulations which quantifies how measures of *ex situ* conservation differ between RADseq and microsatellite markers.

Simulations using different marker types (microsatellites, or "MSAT", and SNP markers, or "DNA")  are performed using the [fastSimcoal2 software](http://cmpg.unibe.ch/software/fastsimcoal27/) and the [strataG](https://github.com/EricArcher/strataG) R package.

**NOTE: Currently this project does not incorporate SNP markers since strataG would not load onto Ash's computer. However, it could (and ideally would) be incorporated, potentially from a different simulation software like SLiM or plink.**

### SimulationOutputs (untouched from the forked repository)
This folder contains 2 subfolders, one for simulation outputs using microsatellite markers ([MSAT_marker](https://github.com/botanic-ash/Morton_SSRvSNP_Simulations/tree/main/SimulationOutputs/MSAT_marker))
and one for simulation outputs using SNP markers ([DNA_marker](https://github.com/botanic-ash/Morton_SSRvSNP_Simulations/tree/main/SimulationOutputs/DNA_marker)). 

Within each of these are folders for each simulation scenario containing Arlequin output files.
Additionally, simulation parameter files, and the log files for simulations, are contained in the top directory of each marker type subfolder.

All scripts which actually run the simulations have been deleted from the forked repository. Can be found in the RScripts folder of the initial [Morton_SSRvSNP_Simulations](https://github.com/HobanLab/Morton_SSRvSNP_Simulations) repository. 

### RScripts
This folder contains the R scripts used to specify simulation parameters, run simulations, convert Arlequin outputs to genind objects, and process those genind objects.

Functions used repeatedly throughout the scripts are declared in [neccessary_functions.R](https://github.com/botanic-ash/Morton_SSRvSNP_Simulations/blob/main/RScripts/necessary_functions.R).

The script [subset_resample_MSAT_AMH_diff_error_rates.R](https://github.com/botanic-ash/Morton_SSRvSNP_Simulations/blob/main/RScripts/subset_resample_MSAT_AMH_diff_error_rates.R) is used to add error to the simulated data, obtain allelic diversity and heterozygosity estimates, and then randomly resample each dataset to model *ex situ* conservation sampling.

Figures made and analyses performed in [figs_and_analyses.R](https://github.com/botanic-ash/Morton_SSRvSNP_Simulations/blob/main/RScripts/figs_and_analyses.R).
