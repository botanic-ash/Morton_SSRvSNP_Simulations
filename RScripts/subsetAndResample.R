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
source("RScripts/GenerateFSCparams.R")
# Alternatively, source the params objects from previously run simulations
# Microsatellites
setwd(paste0(sim.wd,"SimulationOutputs/MSAT_marker"))
MSAT_01pop_migLow.params <- readRDS(dir(pattern = "^params.MSAT_01pop_migLow"))
MSAT_01pop_migHigh.params <- readRDS(dir(pattern = "^params.MSAT_01pop_migHigh"))
MSAT_04pop_migLow.params <- readRDS(dir(pattern = "^params.MSAT_04pop_migLow"))
MSAT_04pop_migHigh.params <- readRDS(dir(pattern = "^params.MSAT_04pop_migHigh"))
MSAT_16pop_migHigh.params <- readRDS(dir(pattern = "^params.MSAT_16pop_migLow"))
MSAT_16pop_migHigh.params <- readRDS(dir(pattern = "^params.MSAT_16pop_migHigh"))
# DNA
setwd(paste0(sim.wd,"SimulationOutputs/DNA_marker"))
DNA_01pop_migLow.params <- readRDS(dir(pattern = "^params.DNA_01pop_migLow"))
DNA_01pop_migHigh.params <- readRDS(dir(pattern = "^params.DNA_01pop_migHigh"))
DNA_04pop_migLow.params <- readRDS(dir(pattern = "^params.DNA_04pop_migLow"))
DNA_04pop_migHigh.params <- readRDS(dir(pattern = "^params.DNA_04pop_migHigh"))
DNA_16pop_migHigh.params <- readRDS(dir(pattern = "^params.DNA_16pop_migLow"))
DNA_16pop_migHigh.params <- readRDS(dir(pattern = "^params.DNA_16pop_migHigh"))

# ---- CONVERT ARLEQUIN FILES TO GENIND ----
# MSAT ----
# Move to MSAT directory
setwd(paste0(sim.wd,"MSAT_marker"))

# Convert files
MSAT_01pop_migLow_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_01pop_migLow/")
MSAT_01pop_migLow_genind <- convertAllArp(arp.path = MSAT_01pop_migLow_arpPath, 
                                          params = MSAT_01pop_migLow.params)

# We need to build a vector where each element is either "garden" or "wild", depending on 
# the results of the sample function

pop(MSAT_01pop_migLow_genind[[1]])

test.1 <- assignGardenSamples(MSAT_01pop_migLow_genind[[1]])
pop(test.1)
exSituRepresentation(test.1)

test.2 <- assignGardenSamples(MSAT_01pop_migLow_genind[[2]])
exSituRepresentation(test.2)

test.3 <- assignGardenSamples(MSAT_01pop_migLow_genind[[3]])
exSituRepresentation(test.3)

test.4 <- assignGardenSamples(MSAT_01pop_migLow_genind[[4]])
exSituRepresentation(test.4)

test.5 <- assignGardenSamples(MSAT_01pop_migLow_genind[[5]])
exSituRepresentation(test.5)
