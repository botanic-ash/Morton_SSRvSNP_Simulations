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
# DNA (1e-8 mutation rate)
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker/data.DNA/"))
# DNA, high mutation (1e-5)
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_highMut_marker/data.DNA/"), prefix="DNA_highMut")

# Specify the proportion of individuals assigned to gardens
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
colnames(MSAT_alleleFreqSummaries) <- colnames(MSAT_repRateSummaries) <-c("mean", "sd")
# DNA: low mutation (1e-8)
DNA_alleleFreqSummaries <- array(dim = c(3, 2, length(DNA_geninds)))
DNA_repRateSummaries <- array(dim = c(5, 2, length(DNA_geninds)))
colnames(DNA_alleleFreqSummaries) <- colnames(DNA_repRateSummaries) <-c("mean", "sd")
# DNA: high mutation (1e-5)
DNA_highMut_alleleFreqSummaries <- array(dim = c(3, 2, length(DNA_highMut_geninds)))
DNA_highMut_repRateSummaries <- array(dim = c(5, 2, length(DNA_highMut_geninds)))
colnames(DNA_highMut_alleleFreqSummaries) <- colnames(DNA_highMut_repRateSummaries) <-c("mean", "sd")

# %%% MSAT ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
for (i in 1:length(MSAT_geninds)){
  # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
  alleleFrequencies <- sapply(MSAT_geninds[[i]], getWildAlleleFreqProportions)
  MSAT_alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
  # Assign individuals to garden population
  MSAT_geninds[[i]] <- lapply(MSAT_geninds[[i]], assignGardenSamples, proportion=gardenRate)
  # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
  representationRates <- sapply(MSAT_geninds[[i]], exSituRepresentation)
  MSAT_repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
}
round(MSAT_alleleFreqSummaries, 2) ; round(MSAT_repRateSummaries, 2)

# %%% DNA, low mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
for (i in 1:length(DNA_geninds)){
  # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
  alleleFrequencies <- sapply(DNA_geninds[[i]], getWildAlleleFreqProportions)
  DNA_alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
  # Assign individuals to garden population
  DNA_geninds[[i]] <- lapply(DNA_geninds[[i]], assignGardenSamples, proportion=gardenRate)
  # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
  representationRates <- sapply(DNA_geninds[[i]], exSituRepresentation)
  DNA_repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
}
round(DNA_alleleFreqSummaries, 2) ; round(DNA_repRateSummaries, 2)

# %%% DNA, high mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
for (i in 1:length(DNA_highMut_geninds)){
  # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
  alleleFrequencies <- sapply(DNA_highMut_geninds[[i]], getWildAlleleFreqProportions)
  DNA_highMut_alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
  # Assign individuals to garden population
  DNA_highMut_geninds[[i]] <- lapply(DNA_highMut_geninds[[i]], assignGardenSamples, proportion=gardenRate)
  # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
  representationRates <- sapply(DNA_highMut_geninds[[i]], exSituRepresentation)
  DNA_highMut_repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
}
round(DNA_highMut_alleleFreqSummaries, 2) ; round(DNA_highMut_repRateSummaries, 2)
