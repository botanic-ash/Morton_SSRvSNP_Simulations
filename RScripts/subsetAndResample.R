# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files, then subsets each file to specify an "ex situ" population
# It then resamples this subset iteratively, and compares the sample to the total genetic
# diversity of the population, and calculates representation metrics

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

sim.wd <- "~/Documents/SSRvSNP/Simulations/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Run simulations
# source("RScripts/GenerateFSCparams.R")
# Alternatively, source the params objects from previously run simulations, using readParams functions
# Microsatellites
readParams_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_marker"))
# DNA
readParams_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker"))

# ---- WORKING: SUMMARIZING REPRESENTATION RATES AND ALLELE FREQUENCY PROPORTIONS ----
# Declare vector of different proportions of individuals that are randomly categorized as "garden"
gardenRates <- c(0.2,0.1,0.05,0.01)
# Matrices to capture allele frequency proportions and ex situ representation rates
MSAT_allFreqProp_mat <- matrix(nrow=6, ncol=3)
MSAT_exSituRep_arr <- array(dim=c(5,5,4))

# CAPTURE ex situ rates for each garden proportion. Report those
# SUMMARIZE allele frequency proportions across all scenarios for each marker

# MSAT ----

# MSAT.scenarios <- c("MSAT_01pop_migLow", "MSAT_01pop_migHigh",
#                     "MSAT_04pop_migLow", "MSAT_04pop_migHigh",
#                     "MSAT_16pop_migLow", "MSAT_16pop_migHigh")
# 
# for(i in 1:length(MSAT.scenarios)){
#   arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/", MSAT.scenarios[i],"/")
#   genindList <- convertAllArp(arp.path = arpPath, params = MSAT.scenarios[i])
# }
# 
# c(MSAT.scenarios[1], ".params")

# 1 pop, mig Low ----
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_01pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_01pop_migLow/")
MSAT_01pop_migLow_genind <- convertAllArp(arp.path = MSAT_01pop_migLow_arpPath, 
                                          params = MSAT_01pop_migLow.params)
# Allele frequency proportions
MSAT_01pop_migLow_freqProportions <- sapply(MSAT_01pop_migLow_genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_01pop_migLow_freqProportions, 1, mean)
# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  browser()
  # Assign individuals to garden population
  MSAT_01pop_migLow_genind <- lapply(MSAT_01pop_migLow_genind, assignGardenSamples, proportion=gardenRates[i])
  # Ex situ representation rates
  MSAT_01pop_migLow_repRates <- sapply(MSAT_01pop_migLow_genind, exSituRepresentation)
  apply(MSAT_01pop_migLow_repRates, 1, mean)
  MSAT_exSituRep_arr[,,i] <- cbind(MSAT_exSituRep_arr[,,i], apply(MSAT_01pop_migLow_repRates, 1, mean))
  
}
MSAT_exSituRep_arr


apply(MSAT_exSituRep_arr, 2, mean, simplify = FALSE)

rowMeans(MSAT_exSituRep_arr, dims = 3)

# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
allFreqProp_mat

apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]

exSituRep_mat

apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd) 

# 1 pop, mig High ----
# Matrices to capture allele frequency proportions and ex situ representation rates
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_01pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_01pop_migHigh/")
MSAT_01pop_migHigh_genind <- convertAllArp(arp.path = MSAT_01pop_migHigh_arpPath, 
                                          params = MSAT_01pop_migHigh.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_01pop_migHigh_genind <- lapply(MSAT_01pop_migHigh_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_01pop_migHigh_freqProportions <- sapply(MSAT_01pop_migHigh_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_01pop_migHigh_freqProportions)
  # Ex situ representation rates
  MSAT_01pop_migHigh_repRates <- sapply(MSAT_01pop_migHigh_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_01pop_migHigh_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)

# 4 pop, mig Low ----
# Matrices to capture allele frequency proportions and ex situ representation rates
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_04pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_04pop_migLow/")
MSAT_04pop_migLow_genind <- convertAllArp(arp.path = MSAT_04pop_migLow_arpPath, 
                                          params = MSAT_04pop_migLow.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_04pop_migLow_genind <- lapply(MSAT_04pop_migLow_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_04pop_migLow_freqProportions <- sapply(MSAT_04pop_migLow_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_04pop_migLow_freqProportions)
  # Ex situ representation rates
  MSAT_04pop_migLow_repRates <- sapply(MSAT_04pop_migLow_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_04pop_migLow_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd) 

# 4 pop, mig High ----
# Matrices to capture allele frequency proportions and ex situ representation rates
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_04pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_04pop_migHigh/")
MSAT_04pop_migHigh_genind <- convertAllArp(arp.path = MSAT_04pop_migHigh_arpPath, 
                                           params = MSAT_04pop_migHigh.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_04pop_migHigh_genind <- lapply(MSAT_04pop_migHigh_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_04pop_migHigh_freqProportions <- sapply(MSAT_04pop_migHigh_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_04pop_migHigh_freqProportions)
  # Ex situ representation rates
  MSAT_04pop_migHigh_repRates <- sapply(MSAT_04pop_migHigh_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_04pop_migHigh_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)

# 16 pop, mig Low ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_16pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_16pop_migLow/")
MSAT_16pop_migLow_genind <- convertAllArp(arp.path = MSAT_16pop_migLow_arpPath, 
                                          params = MSAT_16pop_migLow.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_16pop_migLow_genind <- lapply(MSAT_16pop_migLow_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_16pop_migLow_freqProportions <- sapply(MSAT_16pop_migLow_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_16pop_migLow_freqProportions)
  # Ex situ representation rates
  MSAT_16pop_migLow_repRates <- sapply(MSAT_16pop_migLow_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_16pop_migLow_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd) 

# 16 pop, mig High ----
# Matrices to capture allele frequency proportions and ex situ representation rates
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_16pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_16pop_migHigh/")
MSAT_16pop_migHigh_genind <- convertAllArp(arp.path = MSAT_16pop_migHigh_arpPath, 
                                           params = MSAT_16pop_migHigh.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_16pop_migHigh_genind <- lapply(MSAT_16pop_migHigh_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_16pop_migHigh_freqProportions <- sapply(MSAT_16pop_migHigh_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_16pop_migHigh_freqProportions)
  # Ex situ representation rates
  MSAT_16pop_migHigh_repRates <- sapply(MSAT_16pop_migHigh_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_16pop_migHigh_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)

# DNA ----
# 1 pop, mig Low ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
DNA_01pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/DNA_marker/DNA_01pop_migLow/")
DNA_01pop_migLow_genind <- convertAllArp(arp.path = DNA_01pop_migLow_arpPath, 
                                          params = DNA_01pop_migLow.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_01pop_migLow_genind <- lapply(DNA_01pop_migLow_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_01pop_migLow_freqProportions <- sapply(DNA_01pop_migLow_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_01pop_migLow_freqProportions)
  # Ex situ representation rates
  DNA_01pop_migLow_repRates <- sapply(DNA_01pop_migLow_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_01pop_migLow_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd) 

# 1 pop, mig High ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
DNA_01pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/DNA_marker/DNA_01pop_migHigh/")
DNA_01pop_migHigh_genind <- convertAllArp(arp.path = DNA_01pop_migHigh_arpPath, 
                                           params = DNA_01pop_migHigh.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_01pop_migHigh_genind <- lapply(DNA_01pop_migHigh_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_01pop_migHigh_freqProportions <- sapply(DNA_01pop_migHigh_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_01pop_migHigh_freqProportions)
  # Ex situ representation rates
  DNA_01pop_migHigh_repRates <- sapply(DNA_01pop_migHigh_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_01pop_migHigh_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)

# 4 pop, mig Low ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
DNA_04pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/DNA_marker/DNA_04pop_migLow/")
DNA_04pop_migLow_genind <- convertAllArp(arp.path = DNA_04pop_migLow_arpPath, 
                                          params = DNA_04pop_migLow.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_04pop_migLow_genind <- lapply(DNA_04pop_migLow_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_04pop_migLow_freqProportions <- sapply(DNA_04pop_migLow_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_04pop_migLow_freqProportions)
  # Ex situ representation rates
  DNA_04pop_migLow_repRates <- sapply(DNA_04pop_migLow_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_04pop_migLow_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd) 

# 4 pop, mig High ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
DNA_04pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/DNA_marker/DNA_04pop_migHigh/")
DNA_04pop_migHigh_genind <- convertAllArp(arp.path = DNA_04pop_migHigh_arpPath, 
                                           params = DNA_04pop_migHigh.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_04pop_migHigh_genind <- lapply(DNA_04pop_migHigh_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_04pop_migHigh_freqProportions <- sapply(DNA_04pop_migHigh_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_04pop_migHigh_freqProportions)
  # Ex situ representation rates
  DNA_04pop_migHigh_repRates <- sapply(DNA_04pop_migHigh_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_04pop_migHigh_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)

# 16 pop, mig Low ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
DNA_16pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/DNA_marker/DNA_16pop_migLow/")
DNA_16pop_migLow_genind <- convertAllArp(arp.path = DNA_16pop_migLow_arpPath, 
                                          params = DNA_16pop_migLow.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_16pop_migLow_genind <- lapply(DNA_16pop_migLow_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_16pop_migLow_freqProportions <- sapply(DNA_16pop_migLow_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_16pop_migLow_freqProportions)
  # Ex situ representation rates
  DNA_16pop_migLow_repRates <- sapply(DNA_16pop_migLow_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_16pop_migLow_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd) 

# 16 pop, mig High ----
# Matrices to capture allele frequency proportions and ex situ representation rates
allFreqProp_mat <- matrix(nrow=3)
exSituRep_mat <- matrix(nrow=5)
# Generate genind files, from Arlequin outputs and strataG params objects
DNA_16pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/DNA_marker/DNA_16pop_migHigh/")
DNA_16pop_migHigh_genind <- convertAllArp(arp.path = DNA_16pop_migHigh_arpPath, 
                                           params = DNA_16pop_migHigh.params)

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_16pop_migHigh_genind <- lapply(DNA_16pop_migHigh_genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_16pop_migHigh_freqProportions <- sapply(DNA_16pop_migHigh_genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_16pop_migHigh_freqProportions)
  # Ex situ representation rates
  DNA_16pop_migHigh_repRates <- sapply(DNA_16pop_migHigh_genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_16pop_migHigh_repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)
