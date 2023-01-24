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
library(parallel)
library(RColorBrewer)
library(scales)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs;
# this is a file located on the RAID1 drive, for space reasons, and linked in the home directory)
sim.wd <- "~/Shared/SSRvSNP_Sim/Code/"
setwd(sim.wd)
# Read in relevant functions
source("RScripts/functions_SSRvSNP_Sim.R")
# Parallelism: specify number of cores to use
num_cores <- 8
# Flags for processing different datasets (nInd=1200 and nInd=4800)
nInd_1200_Flag <- 1
nInd_4800_Flag <- 1

# %%% ORIGINAL SIMULATIONS: NIND 1200 %%% ----
# Based on flag value, analyze NIND 1200 dataset
if(nInd_1200_Flag==1){
  print("%%% ANALYZING NIND=1200 DATASET %%%")
  # %%% Read in simulations and process results ----
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
  
  # %%% Assign garden samples ----
  # Specify the proportion of total individuals assigned to gardens
  gardenRate <- 0.05
  # Assign individuals to garden populations. Since we're dealing with a list of lists, we use rapply
  # (We could also use a nested lapply: genind_list <- lapply (genind_list, lapply, function, args))
  MSAT_geninds <- rapply(MSAT_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  DNA_geninds <- rapply(DNA_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  
  # %%% Summarize simulations ----
  # MSAT
  # summarize_simulations(MSAT_geninds)
  # # DNA
  # summarize_simulations(DNA_geninds)
  # Plot a histogram of allele frequencies for one of the DNA smiulations
  # makeAlleleFreqHist(DNA_geninds[[4]][[3]], 
  #                    title="Simulated DNA: 1,000 loci, mutation.rate=5e-6 (4 pops, high migration)")
  
  # %%% Resampling ----
  # Specify number of resampling replicates, to use for both marker types
  num_reps <- 5
  
  # Run resampling in parallel, and save the resampling array result to a specified location on the disk
  # MSAT 
  MSAT_resamplingArrays <- mclapply(MSAT_geninds, Resample_genList, mc.cores = num_cores)
  saveRDS(MSAT_resamplingArrays, 
          file=paste0(sim.wd, "SimulationOutputs/MSAT_marker/data.MSAT/MSAT_N1200_resampArr.Rdata"))
  # DNA
  DNA_resamplingArrays <- mclapply(DNA_geninds, Resample_genList, mc.cores = num_cores)
  saveRDS(DNA_resamplingArrays,
          file=paste0(sim.wd, "SimulationOutputs/DNA_marker/data.DNA/DNA_N1200_resampArr.Rdata"))
  
  # %%% Summarize resampling results ----
  # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
  # MSAT
  # MSAT_min95_Means <- rapply(MSAT_resamplingArrays, resample_min95_mean, how = "list")
  # MSAT_meanValues <- rapply(MSAT_resamplingArrays, resample_meanValues, how = "list")
  # # DNA
  # DNA_min95Values <- rapply(DNA_resamplingArrays, resample_min95_mean, how = "list")
  # DNA_meanValues <- rapply(DNA_resamplingArrays, resample_meanValues, how = "list")
  
  # %%% Plotting ----
  # Pick plot colors (for all plots!)
  # plotColors <- c("red","red4","darkorange3","coral","purple")
  # # Specify the directories to save plots to
  # N1200_MSAT_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/ResamplingCurves_20221215/N1200/MSAT/"
  # N1200_DNA_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/ResamplingCurves_20221215/N1200/DNA/"
  # 
  # # Plotting commands nested in invisible function, to prevent text from being printed
  # # MSAT
  # invisible(rapply(MSAT_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_MSAT_plotDir))
  # # DNA
  # invisible(rapply(DNA_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N1200_DNA_plotDir))
  
  # %%% Build and run linear models ----
  # Capture all explanatory variables into a single data.frame
  # params <- data.frame(expand.grid(n.pop=c(1,4,16), mig.Rate=c(0.001,0.01), marker=c("MSAT", "DNA")))
}

# %%% INCREASED SIMULATIONS: NIND 4800 %%% ----
# Based on flag value, analyze NIND 4800 dataset
if(nInd_4800_Flag==1){
  print("%%% ANALYZING NIND=4800 DATASET %%%")
  # %%% Read in simulations and process results ----
  # Run the simulations
  # source("RScripts/GenerateFSCparams_N4800.R")
  # Alternatively, source the genind objects from previously run simulations, using readGeninds functions
  readGeninds_MSAT(paste0(sim.wd,"SimulationOutputs/MSAT_N4800_marker/data.MSAT/"), prefix = "MSAT_N4800")
  readGeninds_DNA(paste0(sim.wd,"SimulationOutputs/DNA_N4800_marker/data.DNA/"), prefix = "DNA_N4800")
  
  # Combine all scenarios for each marker into a list of genind lists
  MSAT_N4800_geninds <- list(MSAT_N4800_01pop_migLow.genList, MSAT_N4800_01pop_migHigh.genList, 
                             MSAT_N4800_04pop_migLow.genList, MSAT_N4800_04pop_migHigh.genList, 
                             MSAT_N4800_16pop_migLow.genList, MSAT_N4800_16pop_migHigh.genList)
  DNA_N4800_geninds <- list(DNA_N4800_01pop_migLow.genList, DNA_N4800_01pop_migHigh.genList, 
                            DNA_N4800_04pop_migLow.genList, DNA_N4800_04pop_migHigh.genList, 
                            DNA_N4800_16pop_migLow.genList, DNA_N4800_16pop_migHigh.genList)
  
  # %%% Assign garden samples ----
  # Specify the proportion of total individuals assigned to gardens
  gardenRate <- 0.05
  # Assign individuals to garden populations. Since we're dealing with a list of lists, we use rapply
  # (We could also use a nested lapply: genind_list <- lapply (genind_list, lapply, function, args))
  MSAT_N4800_geninds <- rapply(MSAT_N4800_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  DNA_N4800_geninds <- rapply(DNA_N4800_geninds, assignGardenSamples, proportion=gardenRate, how = "list")
  
  # %%% Summarize simulations ----
  # MSAT
  # summarize_simulations(MSAT_N4800_geninds)
  # # DNA
  # summarize_simulations(DNA_N4800_geninds)
  
  # %%% Resampling ----
  # Specify number of resampling replicates, to use for both marker types
  num_reps <- 5
  
  # Run resampling in parallel, and save the resampling array result to a specified location on the disk
  # MSAT 
  MSAT_N4800_resamplingArrays <- mclapply(MSAT_N4800_geninds, Resample_genList, mc.cores = num_cores)
  saveRDS(MSAT_N4800_resamplingArrays, 
          file=paste0(sim.wd, "SimulationOutputs/MSAT_N4800_marker/data.MSAT/MSAT_N4800_resampArr.Rdata"))
  # DNA
  DNA_N4800_resamplingArrays <- mclapply(DNA_N4800_geninds, Resample_genList, mc.cores = num_cores)
  saveRDS(DNA_N4800_resamplingArrays,
          file=paste0(sim.wd, "SimulationOutputs/DNA_N4800_marker/data.DNA/DNA_N4800_resampArr.Rdata"))
  
  # %%% Summarize resampling results ----
  # # Calculate 95% minimum sampling size, and the mean values of each category, for all scenarios
  # MSAT
  # MSAT_N4800_min95_Means <- rapply(MSAT_N4800_resamplingArrays, resample_min95_mean, how = "list")
  # MSAT_N4800_meanValues <- rapply(MSAT_N4800_resamplingArrays, resample_meanValues, how = "list")
  # # DNA
  # DNA_N4800_min95Values <- rapply(DNA_N4800_resamplingArrays, resample_min95_mean, how = "list")
  # DNA_N4800_meanValues <- rapply(DNA_N4800_resamplingArrays, resample_meanValues, how = "list")
  # 
  # # %%% Plotting ----
  # # Pick plot colors (for all plots!)
  # plotColors <- c("red","red4","darkorange3","coral","purple")
  # # Specify the directories to save plots to
  # N4800_MSAT_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/ResamplingCurves_20221215/N4800/MSAT/"
  # N4800_DNA_plotDir <- "~/Documents/SSRvSNP/Simulations/Documentation/Images/ResamplingCurves_20221215/N4800/DNA/"
  # 
  # # Plotting commands nested in invisible function, to prevent text from being printed
  # # MSAT
  # invisible(rapply(MSAT_N4800_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_MSAT_plotDir))
  # # DNA
  # invisible(rapply(DNA_N4800_resamplingArrays, resample_Plot_PNG, colors=plotColors, data.dir=N4800_DNA_plotDir))
}
