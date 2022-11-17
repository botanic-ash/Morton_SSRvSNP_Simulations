# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files generated from previously run fastSimcoal simulations,
# then subsets each file to specify a group of "ex situ" (garden) individuals.
# It then resamples this subset iteratively, and compares the samples to the total genetic
# diversity of the population, and calculates representation metrics. Resampling arrays are 
# created, and then plotted.

# In order to function iteratively over large objects (i.e. lists of genind objects), the functions
# in this script use many apply family functions

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)
library(parallel)
library(RColorBrewer)
library(scales)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs;
# this is a file located on the RAID1 drive, for space reasons, and linked in the home directory)
sim.wd <- "~/Shared/SSRvSNP_Sim/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Parallelism: set up relevant cores; load adegenet and parallel libraries onto the cluster
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("parallel"))

# %%% READ IN SIMULATIONS AND PROCESS RESULTS %%% ----
# Run the simulations
# source("RScripts/GenerateFSCparams.R")
# Alternatively, source the genind objects from previously run simulations, using readGeninds functions
# Microsatellites
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_marker/data.MSAT/"))
# readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_medMut_marker/data.DNA/"), prefix="DNA_medMut")
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker/data.DNA/"))

# Combine all scenarios for each marker into a list of genind lists
MSAT_geninds <- list(MSAT_01pop_migLow.genind, MSAT_01pop_migHigh.genind, MSAT_04pop_migLow.genind,
                     MSAT_04pop_migHigh.genind, MSAT_16pop_migLow.genind, MSAT_16pop_migHigh.genind)
DNA_geninds <- list(DNA_01pop_migLow.genind, DNA_01pop_migHigh.genind, DNA_04pop_migLow.genind,
                           DNA_04pop_migHigh.genind, DNA_16pop_migLow.genind, DNA_16pop_migHigh.genind)

# %%% ASSIGN GARDEN SAMPLES %%% ----
# Specify the proportion of total individuals assigned to gardens
gardenRate <- 0.05
# Assign individuals to garden populations. Since we're dealing with a list of lists, we use rapply
# (We could also use a nested lapply: genind_list <- lapply (genind_list, lapply, function, args))
MSAT_geninds <- rapply(MSAT_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
DNA_geninds <- rapply(DNA_geninds, assignGardenSamples, proportion=gardenRate, how = "list")

# %%% SUMMARIZE SIMULATIONS %%% ----
# MSAT
summarize_simulations(MSAT_geninds, gardenRate = gardenRate)
# DNA
summarize_simulations(DNA_geninds, gardenRate = gardenRate)
makeAlleleFreqHist(DNA_geninds[[4]][[3]], title="Simulated DNA: 1,000 loci, mutation.rate=5e-6 (4 pops, high migration)")

# %%% RESAMPLING %%% ----
# Specify number of replicates, to use for both marker types
num_reps <- 5
# Export relevant functions and variables
clusterExport(cl, varlist = c("getWildFreqs", "getAlleleCategories", "exSitu_Sample", 
                              "exSitu_Resample", "parResample_genind", "num_reps",
                              "MSAT_geninds", "DNA_geninds"))

# MSAT 
MSAT_resamplingArrays <- rapply(MSAT_geninds, parResample_genind, cluster=cl, how = "list")
names(MSAT_resamplingArrays) <- c("MSAT_01pop_migLow","MSAT_01pop_migHigh","MSAT_04pop_migLow",
                                  "MSAT_04pop_migHigh","MSAT_16pop_migLow","MSAT_16pop_migHigh")
# DNA
DNA_resamplingArrays <- rapply(DNA_geninds, parResample_genind, cluster=cl, how = "list")
names(DNA_resamplingArrays) <- c("DNA_01pop_migLow","DNA_01pop_migHigh","DNA_04pop_migLow",
                                 "DNA_04pop_migHigh","DNA_16pop_migLow","DNA_16pop_migHigh")
# Close cores
stopCluster(cl)

# %%% SUMMARIZE RESAMPLING RESULTS %%% ----
# Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
# MSAT
MSAT_min95_Means <- rapply(MSAT_resamplingArrays, resample_min95_mean, how = "list")
MSAT_meanValues <- rapply(MSAT_resamplingArrays, resample_meanValues, how = "list")
# DNA
DNA_min95Values <- rapply(DNA_resamplingArrays, resample_min95_mean, how = "list")
DNA_meanValues <- rapply(DNA_resamplingArrays, resample_meanValues, how = "list")

# %%% PLOTTING %%% 
# Pick plot colors (for all plots!), with transparency for values other than Total
plotColors <- c("red","red4","darkorange3","coral","purple")

# MSAT
rapply(MSAT_resamplingArrays, resample_Plot, colors=plotColors, title="")
# DNA
rapply(DNA_resamplingArrays, resample_Plot, colors=plotColors, title="")
