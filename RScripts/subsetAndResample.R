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
# DNA, medium mutation (1e-6)
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_medMut_marker/data.DNA/"), prefix="DNA_medMut")

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
DNA_medMut_geninds <- list(DNA_medMut_01pop_migLow.genind, DNA_medMut_01pop_migHigh.genind, DNA_medMut_04pop_migLow.genind,
                            DNA_medMut_04pop_migHigh.genind, DNA_medMut_16pop_migLow.genind, DNA_medMut_16pop_migHigh.genind)

# %%% MSAT ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
summarize_simulations(MSAT_geninds, gardenRate = gardenRate)

# %%% DNA, low mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
summarize_simulations(DNA_geninds, gardenRate = gardenRate)

# %%% DNA, high mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
summarize_simulations(DNA_highMut_geninds, gardenRate = gardenRate)

# %%% DNA, medium mutation ----
# Summarize allele frequency proportions and ex situ representation rates across scenarios
summarize_simulations(DNA_medMut_geninds, gardenRate = gardenRate)
