# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files generatedfrom previously run fastSimcoal simulations,
# then subsets each file to specify an a group of "ex situ" (garden) individuals
# It then resamples this subset iteratively, and compares the sample to the total genetic
# diversity of the population, and calculates representation metrics

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)
library(parallel)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs;
# this is a file located on the RAID1 drive, for space reasons, and linked in the home directory)
sim.wd <- "~/Shared/SSRvSNP_Sim/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Parallelism
# Set up relevant cores, and make sure adegenet and parallel libraries are present on the cluster
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)
clusterEvalQ(cl, library("adegenet"))
clusterEvalQ(cl, library("parallel"))

# %%% READ IN SIMULATIONS AND PROCESS RESULTS %%% ----
# source("RScripts/GenerateFSCparams.R")
# Alternatively, source the genind objects from previously run simulations, using readGeninds functions
# Microsatellites
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_marker/data.MSAT/"))
# readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_medMut_marker/data.DNA/"), prefix="DNA_medMut")
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker/data.DNA/"))

# Specify the proportion of individuals assigned to gardens
gardenRate <- 0.05
# Make list of genind objects
MSAT_geninds <- list(MSAT_01pop_migLow.genind, MSAT_01pop_migHigh.genind, MSAT_04pop_migLow.genind,
                     MSAT_04pop_migHigh.genind, MSAT_16pop_migLow.genind, MSAT_16pop_migHigh.genind)
DNA_geninds <- list(DNA_01pop_migLow.genind, DNA_01pop_migHigh.genind, DNA_04pop_migLow.genind,
                           DNA_04pop_migHigh.genind, DNA_16pop_migLow.genind, DNA_16pop_migHigh.genind)
# Simulation summaries
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
                              "exSitu_Resample", "parResample_genind", "num_reps", "gardenRate",
                              "MSAT_geninds", "DNA_geninds"))

# MSAT 
MSAT_01pop_migLow.Resampling <- lapply(MSAT_01pop_migLow.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
MSAT_01pop_migHigh.Resampling <- lapply(MSAT_01pop_migHigh.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
MSAT_04pop_migLow.Resampling <- lapply(MSAT_04pop_migLow.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
MSAT_04pop_migHigh.Resampling <- lapply(MSAT_04pop_migHigh.genind, 
                                        function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
MSAT_16pop_migLow.Resampling <- lapply(MSAT_16pop_migLow.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
MSAT_16pop_migHigh.Resampling <- lapply(MSAT_16pop_migHigh.genind, 
                                        function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))

# DNA
DNA_01pop_migLow.Resampling <- lapply(DNA_01pop_migLow.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
DNA_01pop_migHigh.Resampling <- lapply(DNA_01pop_migHigh.genind, 
                                        function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
DNA_04pop_migLow.Resampling <- lapply(DNA_04pop_migLow.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
DNA_04pop_migHigh.Resampling <- lapply(DNA_04pop_migHigh.genind, 
                                        function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
DNA_16pop_migLow.Resampling <- lapply(DNA_16pop_migLow.genind, 
                                       function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))
DNA_16pop_migHigh.Resampling <- lapply(DNA_16pop_migHigh.genind, 
                                        function(x) parResample_genind(x, gRate=gardenRate, cluster=cl))

# Close cores
stopCluster(cl)
