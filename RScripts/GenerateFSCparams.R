# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% GENERATE FSC PARAMETERS USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses strataG to create the fastsimcoal2 (fsc) parameter files used for the SSR v. SNP comparison project
# Then, those simulation parameter files are run
# After declaring variables (used throughout the script), the code is broken into 
# sections according to which marker type is used for simulations ("MSAT" or "DNA").

# Note that every time this script is run/sourced, new fsc output files are generated

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

# Set working directory to the folder containing fastSimcoal2 simulation outputs
sim.wd <- "~/Documents/SSRvSNP/Simulations/Code/SimulationOutputs/"
setwd(sim.wd)

# ---- VARIABLES ----
num_reps <- 5
fscVersion <- "fsc2709"
# DEMES
# Specify number of total individuals, for all simulations
# Since there are 4 deme and 16 deme scenarios, this value must be divisible by 4 and 16
nInd <- 1200
# 1 Population
demeA <- fscDeme(deme.size = nInd, sample.size = nInd)
demes1 <- fscSettingsDemes(demeA)
# 4 Populations
demeB <- fscDeme(deme.size = nInd/4, sample.size = nInd/4)
demes4 <- fscSettingsDemes(demeB, demeB, demeB, demeB)
# 16 Populations
demeC <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demes16 <- fscSettingsDemes(demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,demeC,
                            demeC,demeC,demeC,demeC)
# MIGRATION
low_mig <- 0.001
high_mig <- 0.01
# 4 Populations
mig.mat.4.Low <- matrix(low_mig, nrow=4, ncol = 4); diag(mig.mat.4.Low) <- 0
mig.mat.4.High <- matrix(high_mig, nrow=4, ncol = 4); diag(mig.mat.4.High) <- 0
mig.mat.4.Final <- matrix(0, nrow=4, ncol = 4)
mig4Low <- fscSettingsMigration(mig.mat.4.Low, mig.mat.4.Final)
mig4High <- fscSettingsMigration(mig.mat.4.High, mig.mat.4.Final)
# 16 Populations
mig.mat.16.Low <- matrix(low_mig, nrow=16, ncol = 16); diag(mig.mat.16.Low) <- 0
mig.mat.16.High <- matrix(high_mig, nrow=16, ncol = 16); diag(mig.mat.16.High) <- 0
mig.mat.16.Final <- matrix(0, nrow=16, ncol = 16)
mig16Low <- fscSettingsMigration(mig.mat.16.Low, mig.mat.16.Final)
mig16High <- fscSettingsMigration(mig.mat.16.High, mig.mat.16.Final)

# HISTORICAL EVENTS
# 4 Populations
hist.event0 <- fscEvent(event.time = 50000, source = 0, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event1 <- fscEvent(event.time = 50000, source = 1, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event2 <- fscEvent(event.time = 50000, source = 2, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event3 <- fscEvent(event.time = 50000, source = 3, sink = 0, prop.migrants = 1, migr.mat = 1)
histEvent4 <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
# 16 Populations
hist.event4 <- fscEvent(event.time = 50000, source = 4, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event5 <- fscEvent(event.time = 50000, source = 5, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event6 <- fscEvent(event.time = 50000, source = 6, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event7 <- fscEvent(event.time = 50000, source = 7, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event8 <- fscEvent(event.time = 50000, source = 8, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event9 <- fscEvent(event.time = 50000, source = 9, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event10 <- fscEvent(event.time = 50000, source = 10, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event11 <- fscEvent(event.time = 50000, source = 11, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event12 <- fscEvent(event.time = 50000, source = 12, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event13 <- fscEvent(event.time = 50000, source = 13, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event14 <- fscEvent(event.time = 50000, source = 14, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event15 <- fscEvent(event.time = 50000, source = 15, sink = 0, prop.migrants = 1, migr.mat = 1)
histEvent16 <- fscSettingsEvents(hist.event0,hist.event1,hist.event2,hist.event3,hist.event4,hist.event5,hist.event6,
                                 hist.event7,hist.event8,hist.event9,hist.event10,hist.event11,hist.event12,
                                 hist.event13,hist.event14,hist.event15)
# MSAT Genetic parameters
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)
# DNA Genetic parameters
dna <- fscBlock_dna(sequence.length = 25, mut.rate = 1e-8)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, num.chrom = 5)

# ---- MICROSATELLITE SIMULATIONS ----
# Outputs are stored within a folder in the parent directory named "MSAT_marker"
setwd(paste0(sim.wd,"MSAT_marker"))
# 1 POPULATION ----
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
MSAT_01pop_migLow.params <- fscWrite(demes = demes1, genetics = MSATgenetics, 
                                     label = "MSAT_01pop_migLow", use.wd=TRUE)
MSAT_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = MSATgenetics, label = "MSAT_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 1 population, low migration")
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 1 population, high migration")
MSAT_01pop_migHigh.params <- fscRun(MSAT_01pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Save params objects to Rdata files, for long term storage
saveRDS(MSAT_01pop_migLow.params, file = paste0("params.MSAT_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_01pop_migHigh.params, file = paste0("params.MSAT_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS ----
# Write parameter files
MSAT_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
MSAT_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 4 populations, low migration")
MSAT_04pop_migLow.params <- fscRun(MSAT_04pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 4 populations, high migration")
MSAT_04pop_migHigh.params <- fscRun(MSAT_04pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Save params objects to Rdata files, for long term storage
saveRDS(MSAT_04pop_migLow.params, file = paste0("params.MSAT_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_04pop_migHigh.params, file = paste0("params.MSAT_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS ----
# Write parameter files
MSAT_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16, 
                                     genetics = MSATgenetics, label = "MSAT_16pop_migLow", use.wd=TRUE)
MSAT_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16, 
                                      genetics = MSATgenetics, label = "MSAT_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("MICROSATELLITES: 16 populations, low migration")
MSAT_16pop_migLow.params <- fscRun(MSAT_16pop_migLow.params, num.sims = num_reps, exec = fscVersion)
print("MICROSATELLITES: 16 populations, high migration")
MSAT_16pop_migHigh.params <- fscRun(MSAT_16pop_migHigh.params, num.sims = num_reps, exec = fscVersion)
# Save params objects to Rdata files, for long term storage
saveRDS(MSAT_16pop_migLow.params, file = paste0("params.MSAT_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(MSAT_16pop_migHigh.params, file = paste0("params.MSAT_16pop_migHigh.",Sys.Date(),".Rdata"))

# ---- DNA SIMULATIONS ----
# Outputs are stored within a folder in the parent directory named "DNA_marker"
setwd(paste0(sim.wd,"DNA_marker"))
# 1 POPULATION ----
# Write parameter files. Make a mighHigh .par file as well, even though it's identical to migLow (with 1 population)
DNA_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 1 population, low migration")
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 1 population, high migration")
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Save params objects to Rdata files, for long term storage
saveRDS(DNA_01pop_migLow.params, file = paste0("params.DNA_01pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_01pop_migHigh.params, file = paste0("params.DNA_01pop_migHigh.",Sys.Date(),".Rdata"))

# 4 POPULATIONS ----
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4, 
                                     genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4, 
                                      genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 4 populations, low migration")
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 4 populations, high migration")
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Save params objects to Rdata files, for long term storage
saveRDS(DNA_04pop_migLow.params, file = paste0("params.DNA_04pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_04pop_migHigh.params, file = paste0("params.DNA_04pop_migHigh.",Sys.Date(),".Rdata"))

# 16 POPULATIONS ----
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16, 
                                     genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16, 
                                      genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
print("DNA: 16 populations, low migration")
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
print("DNA: 16 populations, high migration")
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
# Save params objects to Rdata files, for long term storage
saveRDS(DNA_16pop_migLow.params, file = paste0("params.DNA_16pop_migLow.",Sys.Date(),".Rdata"))
saveRDS(DNA_16pop_migHigh.params, file = paste0("params.DNA_16pop_migHigh.",Sys.Date(),".Rdata"))
