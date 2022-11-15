#!/bin/bash

# Script for running fastSimcoal2 "standalone" (i.e. outside of strataG),
# to demontrate issue of concatenation of DNA marker types

# Create a variable to capture the number of replications, for each simulation instance
reps=5

# Run fastSimcoal2, using the parameter file (build by strataG, in the demo_dnaMarkerIssue.R script)
# -g: simulation genomic data (# of individuals = (# of samples * pop. effective size/2))
# -p: the gametic phase is known in the Arlequin format
# -S: monomorphic DNA sequence sites should be included in the output files
# NOTE: in strataG, it is impossible to run simulations with the -g or -p flags on
fsc2709 -i ./demo_dnaMarkerIssue.par -n $reps -g -p -S
