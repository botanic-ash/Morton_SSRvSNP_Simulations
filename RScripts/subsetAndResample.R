# %%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files generated from previously run fastSimcoal simulations,
# then subsets each file to specify a group of "ex situ" (garden) individuals. The ex situ
# representation (i.e. how well do these garden individuals represent total allelic diversity) is calculated.

# Then, the remaining individuals ("wild") are resampled iteratively, and the allelic diversity
# of sample subsets (in comparison to the whole of wild allelic diversity) is calculated, then plotted

# In order to function iteratively over large objects (i.e. lists of genind objects), the steps
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
readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_marker/data.MSAT/"))
readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_marker/data.DNA/"))

# Combine all scenarios for each marker into a list of genind lists
MSAT_geninds <- list(MSAT_01pop_migLow.genList, MSAT_01pop_migHigh.genList, MSAT_04pop_migLow.genList,
                     MSAT_04pop_migHigh.genList, MSAT_16pop_migLow.genList, MSAT_16pop_migHigh.genList)
DNA_geninds <- list(DNA_01pop_migLow.genList, DNA_01pop_migHigh.genList, DNA_04pop_migLow.genList,
                    DNA_04pop_migHigh.genList, DNA_16pop_migLow.genList, DNA_16pop_migHigh.genList)

# %%% ASSIGN GARDEN SAMPLES %%% ----
# Specify the proportion of total individuals assigned to gardens
gardenRate <- 0.05
# Assign individuals to garden populations. Since we're dealing with a list of lists, we use rapply
# (We could also use a nested lapply: genind_list <- lapply (genind_list, lapply, function, args))
MSAT_geninds <- rapply(MSAT_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
DNA_geninds <- rapply(DNA_geninds, assignGardenSamples, proportion=gardenRate, how = "list")

# %%% SUMMARIZE SIMULATIONS %%% ----
# MSAT
summarize_simulations(MSAT_geninds)
# DNA
summarize_simulations(DNA_geninds)
# Plot a histogram of allele frequencies for one of the DNA smiulations
makeAlleleFreqHist(DNA_geninds[[4]][[3]], 
                   title="Simulated DNA: 1,000 loci, mutation.rate=5e-6 (4 pops, high migration)")

# %%% RESAMPLING %%% ----
# Specify number of resampling replicates, to use for both marker types
num_reps <- 5
# Export relevant functions and variables to the cluster
clusterExport(cl, varlist = c("getWildFreqs", "getAlleleCategories", "exSitu_Sample", 
                              "exSitu_Resample", "parResample_genind", "num_reps",
                              "MSAT_geninds", "DNA_geninds"))

# MSAT 
MSAT_resamplingArrays <- rapply(MSAT_geninds, parResample_genind, reps=num_reps, cluster=cl, how = "list")
# DNA
DNA_resamplingArrays <- rapply(DNA_geninds, parResample_genind, reps=num_reps, cluster=cl, how = "list")

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

# Plotting commands nested in invisible function, to prevent text from being printed
# MSAT
invisible(rapply(MSAT_resamplingArrays, resample_Plot, colors=plotColors))
# DNA
invisible(rapply(DNA_resamplingArrays, resample_Plot, colors=plotColors))
