#!/bin/bash

# Script for starting fastSimcoal2 simulations, based on draft parameter files

# Specify input directory, which contains relevant fasSimcoal2 parameter files
dir=/home/akoontz/Documents/SSRvSNP/Simulations/Code/TEST_ChrStr2x2_fscParams
# In order to avoid space issues, this script needs to be run from the RAID1 directory
# -g indicates genomic data (# of individuals = (# of samples * pop. effective size/2))  
# -p indicates that the gametic phase is known in the Arlequin format
# -S indicates that monomorphic DNA sequence sites should be included in the output files

# Create a variable to capture the number of replications, for each simulation instance
reps=5

# ---RADseq simulations---
# 1 population
fsc2702 -i $dir/params_DNA/DNA_01pops_migHigh.par -n $reps -g -p -S
fsc2702 -i $dir/params_DNA/DNA_01pops_migLow.par -n $reps -g -p -S
# 4 populations
fsc2702 -i $dir/params_DNA/DNA_04pops_migHigh.par -n $reps -g -p -S
fsc2702 -i $dir/params_DNA/DNA_04pops_migLow.par -n $reps -g -p -S
