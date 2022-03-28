#!/bin/bash

# Script for starting fastSimcoal2 simulations, based on draft parameter files

# Specify input directory, which contains relevant fasSimcoal2 parameter files
dir=/home/akoontz/Documents/SSRvSNP/Simulations/Code/draft_20220325_fscParams
# In order to avoid space issues, this script needs to be run from the RAID1 directory

# Create a variable to capture the number of replications, for each simulation instance
reps=5

# ---Microsatellite simulations---
# 1 population
fsc2702 -i $dir/params_MSAT/MSAT_01pops_migHigh.par -n $reps
fsc2702 -i $dir/params_MSAT/MSAT_01pops_migLow.par -n $reps
# 4 populations
fsc2702 -i $dir/params_MSAT/MSAT_04pops_migHigh.par -n $reps
fsc2702 -i $dir/params_MSAT/MSAT_04pops_migLow.par -n $reps
# 16 populations
fsc2702 -i $dir/params_MSAT/MSAT_16pops_migHigh.par -n $reps
fsc2702 -i $dir/params_MSAT/MSAT_16pops_migLow.par -n $reps

# ---RADseq simulations---
# 1 population
fsc2702 -i $dir/params_DNA/DNA_01pops_migHigh.par -n $reps
fsc2702 -i $dir/params_DNA/DNA_01pops_migLow.par -n $reps
# 4 populations
fsc2702 -i $dir/params_DNA/DNA_04pops_migHigh.par -n $reps
fsc2702 -i $dir/params_DNA/DNA_04pops_migLow.par -n $reps
# 16 populations
fsc2702 -i $dir/params_DNA/DNA_16pops_migHigh.par -n $reps
fsc2702 -i $dir/params_DNA/DNA_16pops_migLow.par -n $reps
