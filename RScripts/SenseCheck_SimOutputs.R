# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SENSE CHECK FSC SIMULATION OUTPUTS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script analyzes the simulation outputs from fastSimcoal2 
# to see if different expectations of simulation results are met.
# It first converts the Arlequin output files from fastSimcoal2 to genind objects that can be processed by adegenet

# Then, it processes these genind objects to see if the results of simulations meet certain expectations,
# according to population genetics theory:
# 1. Simulations with more populations should have higher numbers of total alleles
# 2. Simulations with higher migration rates should have lower Fst values
# 3. Analysis of the allele frequency spectra. 

# These checks are made for simulations using both microsatellite ("MSAT") and SNP ("DNA") marker types

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

sim.wd <- "~/Documents/SSRvSNP/Simulations/Code/SimulationOutputs/"
setwd(sim.wd)
# Read in relevant functions
source("../RScripts/SSRvSNP_Sim_functions.R")
# Run simulations
source("../RScripts/GenerateFSCparams.R")

# ---- CONVERT ARLEQUIN FILES TO GENIND ----
# MSAT ----
# Move to MSAT directory
setwd(paste0(sim.wd,"MSAT_marker"))

# Convert files
MSAT_01pop_migLow_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_01pop_migLow/")
MSAT_01pop_migLow_genind <- convertAllArp(arp.path = MSAT_01pop_migLow_arpPath, 
                                          params = MSAT_01pop_migLow.params)
MSAT_01pop_migHigh_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_01pop_migHigh/")
MSAT_01pop_migHigh_genind <- convertAllArp(arp.path = MSAT_01pop_migHigh_arpPath, 
                                           params = MSAT_01pop_migHigh.params)
MSAT_04pop_migLow_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_04pop_migLow/")
MSAT_04pop_migLow_genind <- convertAllArp(arp.path = MSAT_04pop_migLow_arpPath, 
                                          params = MSAT_04pop_migLow.params)
MSAT_04pop_migHigh_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_04pop_migHigh/")
MSAT_04pop_migHigh_genind <- convertAllArp(arp.path = MSAT_04pop_migHigh_arpPath, 
                                           params = MSAT_04pop_migHigh.params)
MSAT_16pop_migLow_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_16pop_migLow/")
MSAT_16pop_migLow_genind <- convertAllArp(arp.path = MSAT_16pop_migLow_arpPath, 
                                          params = MSAT_16pop_migLow.params)
MSAT_16pop_migHigh_arpPath <- paste0(sim.wd,"MSAT_marker/MSAT_16pop_migHigh/")
MSAT_16pop_migHigh_genind <- convertAllArp(arp.path = MSAT_16pop_migHigh_arpPath, 
                                           params = MSAT_16pop_migHigh.params)
# DNA ----
# Move to DNA directory
setwd(paste0(sim.wd,"DNA_marker"))

# Convert files
DNA_01pop_migLow_arpPath <- paste0(sim.wd,"DNA_marker/DNA_01pop_migLow/")
DNA_01pop_migLow_genind <- convertAllArp(arp.path = DNA_01pop_migLow_arpPath, 
                                         params = DNA_01pop_migLow.params)
DNA_01pop_migHigh_arpPath <- paste0(sim.wd,"DNA_marker/DNA_01pop_migHigh/")
DNA_01pop_migHigh_genind <- convertAllArp(arp.path = DNA_01pop_migHigh_arpPath, 
                                          params = DNA_01pop_migHigh.params)
DNA_04pop_migLow_arpPath <- paste0(sim.wd,"DNA_marker/DNA_04pop_migLow/")
DNA_04pop_migLow_genind <- convertAllArp(arp.path = DNA_04pop_migLow_arpPath, 
                                         params = DNA_04pop_migLow.params)
DNA_04pop_migHigh_arpPath <- paste0(sim.wd,"DNA_marker/DNA_04pop_migHigh/")
DNA_04pop_migHigh_genind <- convertAllArp(arp.path = DNA_04pop_migHigh_arpPath, 
                                          params = DNA_04pop_migHigh.params)
DNA_16pop_migLow_arpPath <- paste0(sim.wd,"DNA_marker/DNA_16pop_migLow/")
DNA_16pop_migLow_genind <- convertAllArp(arp.path = DNA_16pop_migLow_arpPath, 
                                         params = DNA_16pop_migLow.params)
DNA_16pop_migHigh_arpPath <- paste0(sim.wd,"DNA_marker/DNA_16pop_migHigh/")
DNA_16pop_migHigh_genind <- convertAllArp(arp.path = DNA_16pop_migHigh_arpPath,                                            params = DNA_16pop_migHigh.params)

# ---- SENSE CHECK ----
# 1. MORE ALLELES IN SCENARIOS WITH MORE POPULATIONS ----
# MSAT ----
mean(sapply(MSAT_01pop_migLow_genind, function(x) ncol(x@tab)))
mean(sapply(MSAT_04pop_migLow_genind, function(x) ncol(x@tab)))
mean(sapply(MSAT_16pop_migLow_genind, function(x) ncol(x@tab)))

mean(sapply(MSAT_01pop_migHigh_genind, function(x) ncol(x@tab)))
mean(sapply(MSAT_04pop_migHigh_genind, function(x) ncol(x@tab)))
mean(sapply(MSAT_16pop_migHigh_genind, function(x) ncol(x@tab)))

# DNA ----
mean(sapply(DNA_01pop_migLow_genind, function(x) ncol(x@tab)))
mean(sapply(DNA_04pop_migLow_genind, function(x) ncol(x@tab)))
mean(sapply(DNA_16pop_migLow_genind, function(x) ncol(x@tab)))

mean(sapply(DNA_01pop_migHigh_genind, function(x) ncol(x@tab)))
mean(sapply(DNA_04pop_migHigh_genind, function(x) ncol(x@tab)))
mean(sapply(DNA_16pop_migHigh_genind, function(x) ncol(x@tab)))

# 2. HIGHER FST FOR SCENARIOS WITH LOWER MIGRATION RATES ----
# MSAT ----
sapply(MSAT_04pop_migLow_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_04pop_migHigh_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

sapply(MSAT_16pop_migLow_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_16pop_migHigh_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# DNA ----
sapply(DNA_04pop_migLow_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(DNA_04pop_migHigh_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

sapply(DNA_16pop_migLow_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) 
sapply(DNA_16pop_migHigh_genind, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# 3. ALLELE FREQUENCY SPECTRA ----
# QUESTION: when we calculate allele frequencies, do we divide by the number of individuals in the
# population? Or, in the entire species? Currently, doing the entire species...
# MSAT ----
MSAT_01pop_migLow_Freqs <- lapply(MSAT_01pop_migLow_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
MSAT_01pop_migHigh_Freqs <- lapply(MSAT_01pop_migHigh_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
MSAT_04pop_migLow_Freqs <- lapply(MSAT_04pop_migLow_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
MSAT_04pop_migHigh_Freqs <- lapply(MSAT_04pop_migHigh_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
MSAT_16pop_migLow_Freqs <- lapply(MSAT_16pop_migLow_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
MSAT_16pop_migHigh_Freqs <- lapply(MSAT_16pop_migHigh_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)

hist(MSAT_01pop_migLow_Freqs[[1]])
hist(MSAT_01pop_migLow_Freqs[[2]])
hist(MSAT_01pop_migLow_Freqs[[3]])
hist(MSAT_01pop_migLow_Freqs[[4]])
hist(MSAT_01pop_migLow_Freqs[[5]])

hist(MSAT_01pop_migHigh_Freqs[[1]])
hist(MSAT_01pop_migHigh_Freqs[[2]])
hist(MSAT_01pop_migHigh_Freqs[[3]])
hist(MSAT_01pop_migHigh_Freqs[[4]])
hist(MSAT_01pop_migHigh_Freqs[[5]])

hist(MSAT_04pop_migLow_Freqs[[1]])
hist(MSAT_04pop_migLow_Freqs[[2]])
hist(MSAT_04pop_migLow_Freqs[[3]])
hist(MSAT_04pop_migLow_Freqs[[4]])
hist(MSAT_04pop_migLow_Freqs[[5]])

hist(MSAT_04pop_migHigh_Freqs[[1]])
hist(MSAT_04pop_migHigh_Freqs[[2]])
hist(MSAT_04pop_migHigh_Freqs[[3]])
hist(MSAT_04pop_migHigh_Freqs[[4]])
hist(MSAT_04pop_migHigh_Freqs[[5]])

hist(MSAT_16pop_migLow_Freqs[[1]])
hist(MSAT_16pop_migLow_Freqs[[2]])
hist(MSAT_16pop_migLow_Freqs[[3]])
hist(MSAT_16pop_migLow_Freqs[[4]])
hist(MSAT_16pop_migLow_Freqs[[5]])

hist(MSAT_16pop_migHigh_Freqs[[1]])
hist(MSAT_16pop_migHigh_Freqs[[2]])
hist(MSAT_16pop_migHigh_Freqs[[3]])
hist(MSAT_16pop_migHigh_Freqs[[4]])
hist(MSAT_16pop_migHigh_Freqs[[5]])

# DNA ----
DNA_01pop_migLow_Freqs <- lapply(DNA_01pop_migLow_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
DNA_01pop_migHigh_Freqs <- lapply(DNA_01pop_migHigh_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
DNA_04pop_migLow_Freqs <- lapply(DNA_04pop_migLow_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
DNA_04pop_migHigh_Freqs <- lapply(DNA_04pop_migHigh_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
DNA_16pop_migLow_Freqs <- lapply(DNA_16pop_migLow_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)
DNA_16pop_migHigh_Freqs <- lapply(DNA_16pop_migHigh_genind, function(x) colSums(x@tab, na.rm = TRUE)/(nInd*2)*100)

hist(DNA_01pop_migLow_Freqs[[1]])
hist(DNA_01pop_migLow_Freqs[[2]])
hist(DNA_01pop_migLow_Freqs[[3]])
hist(DNA_01pop_migLow_Freqs[[4]])
hist(DNA_01pop_migLow_Freqs[[5]])

hist(DNA_01pop_migHigh_Freqs[[1]])
hist(DNA_01pop_migHigh_Freqs[[2]])
hist(DNA_01pop_migHigh_Freqs[[3]])
hist(DNA_01pop_migHigh_Freqs[[4]])
hist(DNA_01pop_migHigh_Freqs[[5]])

hist(DNA_04pop_migLow_Freqs[[1]])
hist(DNA_04pop_migLow_Freqs[[2]])
hist(DNA_04pop_migLow_Freqs[[3]])
hist(DNA_04pop_migLow_Freqs[[4]])
hist(DNA_04pop_migLow_Freqs[[5]])

hist(DNA_04pop_migHigh_Freqs[[1]])
hist(DNA_04pop_migHigh_Freqs[[2]])
hist(DNA_04pop_migHigh_Freqs[[3]])
hist(DNA_04pop_migHigh_Freqs[[4]])
hist(DNA_04pop_migHigh_Freqs[[5]])

hist(DNA_16pop_migLow_Freqs[[1]])
hist(DNA_16pop_migLow_Freqs[[2]])
hist(DNA_16pop_migLow_Freqs[[3]])
hist(DNA_16pop_migLow_Freqs[[4]])
hist(DNA_16pop_migLow_Freqs[[5]])

hist(DNA_16pop_migHigh_Freqs[[1]])
hist(DNA_16pop_migHigh_Freqs[[2]])
hist(DNA_16pop_migHigh_Freqs[[3]])
hist(DNA_16pop_migHigh_Freqs[[4]])
hist(DNA_16pop_migHigh_Freqs[[5]])
