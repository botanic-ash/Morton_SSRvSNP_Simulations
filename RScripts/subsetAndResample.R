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
# Alternatively, source the genind objects from previously run simulations, using readGeninds functions
# Microsatellites
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_marker/data.MSAT/"))
# DNA, low mutation
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker/data.DNA/"))
# DNA, high mutation
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_highMut_marker/data.DNA/"))

# ---- WORKING: EXAMPLES USING A GARDEN RATE OF 5% ----
# 1. Summarize results for MSAT
#     a. Figure out the setup for the matrix that will take in ex Situ representation rates
#     b. Same for the matrix of allele frequency proportions
# 2. Summarize results for DNA, low mutation (1e-8)
# 3. Summarize results for DNA, high mutation (1e-5)
gardenRate <- 0.05
# Make list of genind objects
# (This could be a part of the readGeninds_MSAT/DNA functions...consider updating later)
MSAT_geninds <- list(MSAT_01pop_migLow.genind, MSAT_01pop_migHigh.genind, MSAT_04pop_migLow.genind,
                     MSAT_04pop_migHigh.genind, MSAT_16pop_migLow.genind, MSAT_16pop_migHigh.genind)
DNA_geninds <- list(DNA_01pop_migLow.genind, DNA_01pop_migHigh.genind, DNA_04pop_migLow.genind,
                    DNA_04pop_migHigh.genind, DNA_16pop_migLow.genind, DNA_16pop_migHigh.genind)
DNA_highMut_geninds <- list(DNA_highMut_01pop_migLow.genind, DNA_highMut_01pop_migHigh.genind, DNA_highMut_04pop_migLow.genind,
                    DNA_highMut_04pop_migHigh.genind, DNA_highMut_16pop_migLow.genind, DNA_highMut_16pop_migHigh.genind)
# Make arrays for allele frequency proportions and ex situ representation rates (both markers)
# MSATs
MSAT_alleleFreqSummaries <- array(dim = c(3, 2, length(MSAT_geninds)))
MSAT_repRateSummaries <- array(dim = c(5, 2, length(MSAT_geninds)))
# DNA: low mutation (1e-8)
DNA_alleleFreqSummaries <- array(dim = c(3, 2, length(DNA_geninds)))
DNA_repRateSummaries <- array(dim = c(5, 2, length(DNA_geninds)))
# DNA: high mutation (1e-8)
DNA_highMut_alleleFreqSummaries <- array(dim = c(3, 2, length(DNA_highMut_geninds)))
DNA_highMut_repRateSummaries <- array(dim = c(5, 2, length(DNA_highMut_geninds)))

# %%% MSAT ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
for (i in 1:length(MSAT_geninds)){
  # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
  alleleFrequencies <- sapply(MSAT_geninds[[i]], getAlleleFreqProportions)
  MSAT_alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
  # Assign individuals to garden population
  MSAT_geninds[[i]] <- lapply(MSAT_geninds[[i]], assignGardenSamples, proportion=gardenRate)
  # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
  representationRates <- sapply(MSAT_geninds[[i]], exSituRepresentation)
  MSAT_repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
}

# %%% DNA, low mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
for (i in 1:length(DNA_geninds)){
  # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
  alleleFrequencies <- sapply(DNA_geninds[[i]], getAlleleFreqProportions)
  DNA_alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
  # Assign individuals to garden population
  DNA_geninds[[i]] <- lapply(DNA_geninds[[i]], assignGardenSamples, proportion=gardenRate)
  # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
  representationRates <- sapply(DNA_geninds[[i]], exSituRepresentation)
  DNA_repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
}

# %%% DNA, high mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
for (i in 1:length(DNA_highMut_geninds)){
  # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
  alleleFrequencies <- sapply(DNA_highMut_geninds[[i]], getAlleleFreqProportions)
  DNA_highMut_alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
  # Assign individuals to garden population
  DNA_highMut_geninds[[i]] <- lapply(DNA_highMut_geninds[[i]], assignGardenSamples, proportion=gardenRate)
  # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
  representationRates <- sapply(DNA_highMut_geninds[[i]], exSituRepresentation)
  DNA_highMut_repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
}

# 1 pop, mig Low ----
# Assign individuals to garden population
MSAT_01pop_migLow.genind <- lapply(MSAT_01pop_migLow.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
MSAT_01pop_migLow.repRates <- sapply(MSAT_01pop_migLow.genind, exSituRepresentation); print(MSAT_01pop_migLow.repRates)

# 1 pop, mig High ----
# Assign individuals to garden population
MSAT_01pop_migHigh.genind <- lapply(MSAT_01pop_migHigh.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
MSAT_01pop_migHigh.repRates <- sapply(MSAT_01pop_migHigh.genind, exSituRepresentation); print(MSAT_01pop_migHigh.repRates)

# 4 pops, mig Low ----
# Assign individuals to garden population
MSAT_04pop_migLow.genind <- lapply(MSAT_04pop_migLow.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
MSAT_04pop_migLow.repRates <- sapply(MSAT_04pop_migLow.genind, exSituRepresentation); print(MSAT_04pop_migLow.repRates)

# 4 pops, mig High ----
# Assign individuals to garden population
MSAT_04pop_migHigh.genind <- lapply(MSAT_04pop_migHigh.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
MSAT_04pop_migHigh.repRates <- sapply(MSAT_04pop_migHigh.genind, exSituRepresentation); print(MSAT_04pop_migHigh.repRates)

# 16 pops, mig Low ----
# Assign individuals to garden population
MSAT_16pop_migLow.genind <- lapply(MSAT_16pop_migLow.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
MSAT_16pop_migLow.repRates <- sapply(MSAT_16pop_migLow.genind, exSituRepresentation); print(MSAT_16pop_migLow.repRates)

# 16 pops, mig High ----
# Assign individuals to garden population
MSAT_16pop_migHigh.genind <- lapply(MSAT_16pop_migHigh.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
MSAT_16pop_migHigh.repRates <- sapply(MSAT_16pop_migHigh.genind, exSituRepresentation); print(MSAT_16pop_migHigh.repRates)
summarize_exSituRepresentation(MSAT_01pop_migHigh.repRates)

# %%% DNA ----
# 1 pop, mig Low 
DNA_01pop_migLow.freqProportions <- sapply(DNA_01pop_migLow.genind, getAlleleFreqProportions)
summarize_alleleFreqProportions(DNA_01pop_migLow.freqProportions)

# 1 pop, mig High 
DNA_01pop_migHigh.freqProportions <- sapply(DNA_01pop_migHigh.genind, getAlleleFreqProportions)
summarize_alleleFreqProportions(DNA_01pop_migHigh.freqProportions)

# 4 pops, mig Low 
DNA_04pop_migLow.freqProportions <- sapply(DNA_04pop_migLow.genind, getAlleleFreqProportions)
summarize_alleleFreqProportions(DNA_04pop_migLow.freqProportions)

# 4 pops, mig High 
DNA_04pop_migHigh.freqProportions <- sapply(DNA_04pop_migHigh.genind, getAlleleFreqProportions)
summarize_alleleFreqProportions(DNA_04pop_migHigh.freqProportions)

# 16 pops, mig Low 
DNA_16pop_migLow.freqProportions <- sapply(DNA_16pop_migLow.genind, getAlleleFreqProportions)
summarize_alleleFreqProportions(DNA_16pop_migLow.freqProportions)

# 16 pops, mig High 
DNA_16pop_migHigh.freqProportions <- sapply(DNA_16pop_migHigh.genind, getAlleleFreqProportions)
summarize_alleleFreqProportions(DNA_16pop_migHigh.freqProportions)


# 1 pop, mig Low ----
# Assign individuals to garden population
DNA_01pop_migLow.genind <- lapply(DNA_01pop_migLow.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
DNA_01pop_migLow.repRates <- sapply(DNA_01pop_migLow.genind, exSituRepresentation); print(DNA_01pop_migLow.repRates)

# 1 pop, mig High ----
# Assign individuals to garden population
DNA_01pop_migHigh.genind <- lapply(DNA_01pop_migHigh.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
DNA_01pop_migHigh.repRates <- sapply(DNA_01pop_migHigh.genind, exSituRepresentation); print(DNA_01pop_migHigh.repRates)

# 4 pops, mig Low ----
# Assign individuals to garden population
DNA_04pop_migLow.genind <- lapply(DNA_04pop_migLow.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
DNA_04pop_migLow.repRates <- sapply(DNA_04pop_migLow.genind, exSituRepresentation); print(DNA_04pop_migLow.repRates)

# 4 pops, mig High ----
# Assign individuals to garden population
DNA_04pop_migHigh.genind <- lapply(DNA_04pop_migHigh.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
DNA_04pop_migHigh.repRates <- sapply(DNA_04pop_migHigh.genind, exSituRepresentation); print(DNA_04pop_migHigh.repRates)

# 16 pops, mig Low ----
# Assign individuals to garden population
DNA_16pop_migLow.genind <- lapply(DNA_16pop_migLow.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
DNA_16pop_migLow.repRates <- sapply(DNA_16pop_migLow.genind, exSituRepresentation); print(DNA_16pop_migLow.repRates)

# 16 pops, mig High ----
# Assign individuals to garden population
DNA_16pop_migHigh.genind <- lapply(DNA_16pop_migHigh.genind, assignGardenSamples, proportion=gardenRate)
# Ex situ representation rates
DNA_16pop_migHigh.repRates <- sapply(DNA_16pop_migHigh.genind, exSituRepresentation); print(DNA_16pop_migHigh.repRates)

# ---- WORKING: SUMMARIZING REPRESENTATION RATES AND ALLELE FREQUENCY PROPORTIONS ----
# Declare vector of different proportions of individuals that are randomly categorized as "garden"
gardenRates <- c(0.2,0.1,0.05,0.01)
# Matrices to capture allele frequency proportions and ex situ representation rates
MSAT_allFreqProp_mat <- matrix(nrow=6, ncol=3)
MSAT_exSituRep_arr <- array(dim=c(5,5,4))

# CAPTURE ex situ rates for each garden proportion. Report those
# SUMMARIZE allele frequency proportions across all scenarios for each marker

MSAT_alleleProps_Mat <- matrix(nrow=6,ncol=3)
DNA_alleleProps_Mat <- matrix(nrow=6,ncol=3)

# 1 pop, mig Low ----
# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_01pop_migLow_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_01pop_migLow/")
MSAT_01pop_migLow.genind <- convertAllArp(arp.path = MSAT_01pop_migLow_arpPath, 
                                          params = MSAT_01pop_migLow.params)
# Allele frequency proportions
MSAT_01pop_migLow_freqProportions <- sapply(MSAT_01pop_migLow.genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_01pop_migLow_freqProportions, 1, mean)
MSAT_alleleProps_Mat[1,] <- MSAT_allFreqProp_mat
# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  browser()
  # Assign individuals to garden population
  MSAT_01pop_migLow.genind <- lapply(MSAT_01pop_migLow.genind, assignGardenSamples, proportion=gardenRates[i])
  # Ex situ representation rates
  MSAT_01pop_migLow.repRates <- sapply(MSAT_01pop_migLow.genind, exSituRepresentation)
  apply(MSAT_01pop_migLow.repRates, 1, mean)
  MSAT_exSituRep_arr[,,i] <- cbind(MSAT_exSituRep_arr[,,i], apply(MSAT_01pop_migLow.repRates, 1, mean))
  
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

# Generate genind files, from Arlequin outputs and strataG params objects
MSAT_01pop_migHigh_arpPath <- paste0(sim.wd,"SimulationOutputs/MSAT_marker/MSAT_01pop_migHigh/")
MSAT_01pop_migHigh.genind <- convertAllArp(arp.path = MSAT_01pop_migHigh_arpPath, 
                                          params = MSAT_01pop_migHigh.params)
# Allele frequency proportions
MSAT_01pop_migHigh_freqProportions <- sapply(MSAT_01pop_migHigh.genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_01pop_migHigh_freqProportions, 1, mean)
MSAT_alleleProps_Mat[2,] <- MSAT_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_01pop_migHigh.genind <- lapply(MSAT_01pop_migHigh.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_01pop_migHigh_freqProportions <- sapply(MSAT_01pop_migHigh.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_01pop_migHigh_freqProportions)
  # Ex situ representation rates
  MSAT_01pop_migHigh.repRates <- sapply(MSAT_01pop_migHigh.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_01pop_migHigh.repRates)
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
MSAT_04pop_migLow.genind <- convertAllArp(arp.path = MSAT_04pop_migLow_arpPath, 
                                          params = MSAT_04pop_migLow.params)
# Allele frequency proportions
MSAT_04pop_migLow_freqProportions <- sapply(MSAT_04pop_migLow.genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_04pop_migLow_freqProportions, 1, mean)
MSAT_alleleProps_Mat[3,] <- MSAT_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_04pop_migLow.genind <- lapply(MSAT_04pop_migLow.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_04pop_migLow_freqProportions <- sapply(MSAT_04pop_migLow.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_04pop_migLow_freqProportions)
  # Ex situ representation rates
  MSAT_04pop_migLow.repRates <- sapply(MSAT_04pop_migLow.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_04pop_migLow.repRates)
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
MSAT_04pop_migHigh.genind <- convertAllArp(arp.path = MSAT_04pop_migHigh_arpPath, 
                                           params = MSAT_04pop_migHigh.params)
# Allele frequency proportions
MSAT_04pop_migHigh_freqProportions <- sapply(MSAT_04pop_migHigh.genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_04pop_migHigh_freqProportions, 1, mean)
MSAT_alleleProps_Mat[4,] <- MSAT_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_04pop_migHigh.genind <- lapply(MSAT_04pop_migHigh.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_04pop_migHigh_freqProportions <- sapply(MSAT_04pop_migHigh.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_04pop_migHigh_freqProportions)
  # Ex situ representation rates
  MSAT_04pop_migHigh.repRates <- sapply(MSAT_04pop_migHigh.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_04pop_migHigh.repRates)
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
MSAT_16pop_migLow.genind <- convertAllArp(arp.path = MSAT_16pop_migLow_arpPath, 
                                          params = MSAT_16pop_migLow.params)
# Allele frequency proportions
MSAT_16pop_migLow_freqProportions <- sapply(MSAT_16pop_migLow.genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_16pop_migLow_freqProportions, 1, mean)
MSAT_alleleProps_Mat[5,] <- MSAT_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_16pop_migLow.genind <- lapply(MSAT_16pop_migLow.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_16pop_migLow_freqProportions <- sapply(MSAT_16pop_migLow.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_16pop_migLow_freqProportions)
  # Ex situ representation rates
  MSAT_16pop_migLow.repRates <- sapply(MSAT_16pop_migLow.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_16pop_migLow.repRates)
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
MSAT_16pop_migHigh.genind <- convertAllArp(arp.path = MSAT_16pop_migHigh_arpPath, 
                                           params = MSAT_16pop_migHigh.params)
# Allele frequency proportions
MSAT_16pop_migHigh_freqProportions <- sapply(MSAT_16pop_migHigh.genind, getAlleleFreqProportions)
MSAT_allFreqProp_mat <- apply(MSAT_16pop_migHigh_freqProportions, 1, mean)
MSAT_alleleProps_Mat[6,] <- MSAT_allFreqProp_mat


# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  MSAT_16pop_migHigh.genind <- lapply(MSAT_16pop_migHigh.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  MSAT_16pop_migHigh_freqProportions <- sapply(MSAT_16pop_migHigh.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, MSAT_16pop_migHigh_freqProportions)
  # Ex situ representation rates
  MSAT_16pop_migHigh.repRates <- sapply(MSAT_16pop_migHigh.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, MSAT_16pop_migHigh.repRates)
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
DNA_01pop_migLow.genind <- convertAllArp(arp.path = DNA_01pop_migLow_arpPath, 
                                          params = DNA_01pop_migLow.params)
# Allele frequency proportions
DNA_01pop_migLow_freqProportions <- sapply(DNA_01pop_migLow.genind, getAlleleFreqProportions)
DNA_allFreqProp_mat <- apply(DNA_01pop_migLow_freqProportions, 1, mean)
DNA_alleleProps_Mat[1,] <- DNA_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_01pop_migLow.genind <- lapply(DNA_01pop_migLow.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_01pop_migLow_freqProportions <- sapply(DNA_01pop_migLow.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_01pop_migLow_freqProportions)
  # Ex situ representation rates
  DNA_01pop_migLow.repRates <- sapply(DNA_01pop_migLow.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_01pop_migLow.repRates)
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
DNA_01pop_migHigh.genind <- convertAllArp(arp.path = DNA_01pop_migHigh_arpPath, 
                                           params = DNA_01pop_migHigh.params)
# Allele frequency proportions
DNA_01pop_migHigh_freqProportions <- sapply(DNA_01pop_migHigh.genind, getAlleleFreqProportions)
DNA_allFreqProp_mat <- apply(DNA_01pop_migHigh_freqProportions, 1, mean)
DNA_alleleProps_Mat[2,] <- DNA_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_01pop_migHigh.genind <- lapply(DNA_01pop_migHigh.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_01pop_migHigh_freqProportions <- sapply(DNA_01pop_migHigh.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_01pop_migHigh_freqProportions)
  # Ex situ representation rates
  DNA_01pop_migHigh.repRates <- sapply(DNA_01pop_migHigh.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_01pop_migHigh.repRates)
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
DNA_04pop_migLow.genind <- convertAllArp(arp.path = DNA_04pop_migLow_arpPath, 
                                          params = DNA_04pop_migLow.params)
# Allele frequency proportions
DNA_04pop_migLow_freqProportions <- sapply(DNA_04pop_migLow.genind, getAlleleFreqProportions)
DNA_allFreqProp_mat <- apply(DNA_04pop_migLow_freqProportions, 1, mean)
DNA_alleleProps_Mat[3,] <- DNA_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_04pop_migLow.genind <- lapply(DNA_04pop_migLow.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_04pop_migLow_freqProportions <- sapply(DNA_04pop_migLow.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_04pop_migLow_freqProportions)
  # Ex situ representation rates
  DNA_04pop_migLow.repRates <- sapply(DNA_04pop_migLow.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_04pop_migLow.repRates)
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
DNA_04pop_migHigh.genind <- convertAllArp(arp.path = DNA_04pop_migHigh_arpPath, 
                                           params = DNA_04pop_migHigh.params)
# Allele frequency proportions
DNA_04pop_migHigh_freqProportions <- sapply(DNA_04pop_migHigh.genind, getAlleleFreqProportions)
DNA_allFreqProp_mat <- apply(DNA_04pop_migHigh_freqProportions, 1, mean)
DNA_alleleProps_Mat[4,] <- DNA_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_04pop_migHigh.genind <- lapply(DNA_04pop_migHigh.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_04pop_migHigh_freqProportions <- sapply(DNA_04pop_migHigh.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_04pop_migHigh_freqProportions)
  # Ex situ representation rates
  DNA_04pop_migHigh.repRates <- sapply(DNA_04pop_migHigh.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_04pop_migHigh.repRates)
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
DNA_16pop_migLow.genind <- convertAllArp(arp.path = DNA_16pop_migLow_arpPath, 
                                          params = DNA_16pop_migLow.params)
# Allele frequency proportions
DNA_16pop_migLow_freqProportions <- sapply(DNA_16pop_migLow.genind, getAlleleFreqProportions)
DNA_allFreqProp_mat <- apply(DNA_16pop_migLow_freqProportions, 1, mean)
DNA_alleleProps_Mat[5,] <- DNA_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_16pop_migLow.genind <- lapply(DNA_16pop_migLow.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_16pop_migLow_freqProportions <- sapply(DNA_16pop_migLow.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_16pop_migLow_freqProportions)
  # Ex situ representation rates
  DNA_16pop_migLow.repRates <- sapply(DNA_16pop_migLow.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_16pop_migLow.repRates)
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
DNA_16pop_migHigh.genind <- convertAllArp(arp.path = DNA_16pop_migHigh_arpPath, 
                                           params = DNA_16pop_migHigh.params)
# Allele frequency proportions
DNA_16pop_migHigh_freqProportions <- sapply(DNA_16pop_migHigh.genind, getAlleleFreqProportions)
DNA_allFreqProp_mat <- apply(DNA_16pop_migHigh_freqProportions, 1, mean)
DNA_alleleProps_Mat[6,] <- DNA_allFreqProp_mat

# Loop through different garden proportions, appending results to matrices
for(i in 1:length(gardenRates)){
  # Assign individuals to garden population
  DNA_16pop_migHigh.genind <- lapply(DNA_16pop_migHigh.genind, assignGardenSamples, proportion=gardenRates[i])
  # Allele frequency proportions
  DNA_16pop_migHigh_freqProportions <- sapply(DNA_16pop_migHigh.genind, getAlleleFreqProportions)
  allFreqProp_mat <- cbind(allFreqProp_mat, DNA_16pop_migHigh_freqProportions)
  # Ex situ representation rates
  DNA_16pop_migHigh.repRates <- sapply(DNA_16pop_migHigh.genind, exSituRepresentation)
  exSituRep_mat <- cbind(exSituRep_mat, DNA_16pop_migHigh.repRates)
}
# Summarize allele frequency proportions
allFreqProp_mat <- allFreqProp_mat[,-1]
apply(allFreqProp_mat, 1, mean); apply(allFreqProp_mat, 1, sd) 
# Summarize ex situ representation rates
exSituRep_mat <- exSituRep_mat[,-1]
apply(exSituRep_mat, 1, mean); apply(exSituRep_mat, 1, sd)
