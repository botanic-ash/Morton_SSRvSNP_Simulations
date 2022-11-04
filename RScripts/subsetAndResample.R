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

# Set working directory to local GitHub repository 
sim.wd <- "/RAID1/Simulations/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Run simulations
# source("RScripts/GenerateFSCparams.R")
# Alternatively, source the genind objects from previously run simulations, using readGeninds functions
# Microsatellites
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_marker/data.MSAT/"))
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_medMut_marker/data.DNA/"), prefix="DNA_medMut")
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker/data.DNA/"))
# DNA, medium-high mutation (5e-5)
# readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_medHighMut_marker/data.DNA/"), prefix="DNA_medHighMut")

# Specify the proportion of individuals assigned to gardens
gardenRate <- 0.05
# Make list of genind objects
MSAT_geninds <- list(MSAT_01pop_migLow.genind, MSAT_01pop_migHigh.genind, MSAT_04pop_migLow.genind,
                     MSAT_04pop_migHigh.genind, MSAT_16pop_migLow.genind, MSAT_16pop_migHigh.genind)
# DNA_geninds <- list(DNA_01pop_migLow.genind, DNA_01pop_migHigh.genind, DNA_04pop_migLow.genind,
#                     DNA_04pop_migHigh.genind, DNA_16pop_migLow.genind, DNA_16pop_migHigh.genind)
# DNA_highMut_geninds <- list(DNA_highMut_01pop_migLow.genind, DNA_highMut_01pop_migHigh.genind, DNA_highMut_04pop_migLow.genind,
#                     DNA_highMut_04pop_migHigh.genind, DNA_highMut_16pop_migLow.genind, DNA_highMut_16pop_migHigh.genind)
DNA_medMut_geninds <- list(DNA_medMut_01pop_migLow.genind, DNA_medMut_01pop_migHigh.genind, DNA_medMut_04pop_migLow.genind,
                            DNA_medMut_04pop_migHigh.genind, DNA_medMut_16pop_migLow.genind, DNA_medMut_16pop_migHigh.genind)
DNA_geninds <- list(DNA_01pop_migLow.genind, DNA_01pop_migHigh.genind, DNA_04pop_migLow.genind,
                           DNA_04pop_migHigh.genind, DNA_16pop_migLow.genind, DNA_16pop_migHigh.genind)
# DNA_medHighMut_geninds <- list(DNA_medHighMut_01pop_migLow.genind, DNA_medHighMut_01pop_migHigh.genind, DNA_medHighMut_04pop_migLow.genind,
#                            DNA_medHighMut_04pop_migHigh.genind, DNA_medHighMut_16pop_migLow.genind, DNA_medHighMut_16pop_migHigh.genind)

# %%% MSAT ----
# Range.constraint=15
summarize_simulations(MSAT_geninds, gardenRate = gardenRate)

# Range.constraint=10
# summarize_simulations(MSAT_rang10_geninds, gardenRate = gardenRate)

# Range.constraint=10, higher mutation rate (1e-3) ----
# summarize_simulations(MSAT_highMut_geninds, gardenRate = gardenRate)

# %%% DNA ----
# Lower mutation (1e-8)
# summarize_simulations(DNA_geninds, gardenRate = gardenRate)

# Higher mutation (1e-5)
# summarize_simulations(DNA_highMut_geninds, gardenRate = gardenRate)

# Medium mutation (1e-6)
summarize_simulations(DNA_medMut_geninds, gardenRate = gardenRate)

# Larger numbers of SNP loci
summarize_simulations(DNA_geninds, gardenRate = gardenRate)
makeAlleleFreqHist(DNA_geninds[[4]][[3]], title="Simulated DNA: 1,000 loci, mutation.rate=5e-6 (4 pops, high migration)")

# Medium-high mutation (5e-6)
# summarize_simulations(DNA_medHighMut_geninds, gardenRate = gardenRate)

# !!! WORKING !!!
# %%% RESAMPLING %%% ----
# Set up relevant cores, and make sure adegenet library is present on cluster
num_cores <- detectCores() - 1 ; cl <- makeCluster(num_cores)
clusterEvalQ(cl, library("adegenet"))

# MSAT ----
# Specify number of replicates
num_reps <- 50
# Export relevant functions and variables
clusterExport(cl, varlist = c("getAlleleCategories", "exSitu_Sample", "exSitu_Resample", 
                              "num_reps", "MSAT_geninds"))
# Run resampling in parallel, to generate an array
samplingResults_MSAT.01pop.migLow.1 <- 
  parSapply(cl, 1:num_reps, 
            function(a) exSitu_Resample(gen.obj = MSAT_01pop_migLow.genind[[1]]), 
            simplify = "array")



# Close cores
stopCluster(cl)
