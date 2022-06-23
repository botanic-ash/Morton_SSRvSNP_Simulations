# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FSC PARAMETERS USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses strataG to create the fastsimcoal2 (fsc) parameter files used for the SSR v. SNP comparison project
# Then, those parameter files are run, and the outputs stored
# After declaring variables (used throughout the script), the code is broken into an MSAT section and a DNA section

# TO DO: to optimize this script, check out the "Predefinted Values" section of the fscTutorials() page

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/")

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
demeB_1 <- fscDeme(deme.size = nInd/4, sample.size = nInd/4)
demeB_2 <- fscDeme(deme.size = nInd/4, sample.size = nInd/4)
demeB_3 <- fscDeme(deme.size = nInd/4, sample.size = nInd/4)
demeB_4 <- fscDeme(deme.size = nInd/4, sample.size = nInd/4)
demes4 <- fscSettingsDemes(demeB_1, demeB_2, demeB_3, demeB_4)
# 16 Populations
demeC_1 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_2 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_3 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_4 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_5 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_6 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_7 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_8 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_9 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_10 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_11 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_12 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_13 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_14 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_15 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demeC_16 <- fscDeme(deme.size = nInd/16, sample.size = nInd/16)
demes16 <- fscSettingsDemes(demeC_1,demeC_2,demeC_3,demeC_4,demeC_5,demeC_6,demeC_7,demeC_8,demeC_9,demeC_10,demeC_11,demeC_12,
                            demeC_13,demeC_14,demeC_15,demeC_16)
# MIGRATION
low_mig <- 0.001
high_mig <- 0.01
mig.mat.4.Low <- matrix(low_mig, nrow=4, ncol = 4); diag(mig.mat.4.Low) <- 0
mig.mat.4.High <- matrix(high_mig, nrow=4, ncol = 4); diag(mig.mat.4.High) <- 0
mig.mat.4.Final <- matrix(0, nrow=4, ncol = 4)
mig4Low <- fscSettingsMigration(mig.mat.4.Low, mig.mat.4.Final)
mig4High <- fscSettingsMigration(mig.mat.4.High, mig.mat.4.Final)
mig.mat.16.Low <- matrix(low_mig, nrow=16, ncol = 16); diag(mig.mat.16.Low) <- 0
mig.mat.16.High <- matrix(high_mig, nrow=16, ncol = 16); diag(mig.mat.16.High) <- 0
mig.mat.16.Final <- matrix(0, nrow=16, ncol = 16)
mig16Low <- fscSettingsMigration(mig.mat.16.Low, mig.mat.16.Final)
mig16High <- fscSettingsMigration(mig.mat.16.High, mig.mat.16.Final)
# Historical events
hist.event0 <- fscEvent(event.time = 50000, source = 0, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event1 <- fscEvent(event.time = 50000, source = 1, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event2 <- fscEvent(event.time = 50000, source = 2, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event3 <- fscEvent(event.time = 50000, source = 3, sink = 0, prop.migrants = 1, migr.mat = 1)
histEvent4 <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
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
histEvent16 <- fscSettingsEvents(hist.event0,hist.event1,hist.event2,hist.event3,hist.event4,hist.event5,hist.event6,hist.event7,
                                 hist.event8,hist.event9,hist.event10,hist.event11,hist.event12,hist.event13,hist.event14,hist.event15)
# MSAT Genetic parameters
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)
# DNA Genetic parameters
dna <- fscBlock_dna(sequence.length = 25, mut.rate = 1e-5)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, num.chrom = 5)

# ---- MICROSATELLITE SIMULATIONS ----
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/MSAT_marker/")
# 1 POPULATION ----
# Write parameter files. We make a mighHigh .par file as well, even though it's identical to migLow (with one population)
MSAT_01pop_migLow.params <- fscWrite(demes = demes1, genetics = MSATgenetics, label = "MSAT_01pop_migLow", use.wd=TRUE)
MSAT_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = MSATgenetics, label = "MSAT_01pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps, exec = fscVersion)
MSAT_01pop_migHigh.params <- fscRun(MSAT_01pop_migHigh.params, num.sims = num_reps, exec = fscVersion)

# 4 POPULATIONS ----
# Write parameter files
MSAT_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
MSAT_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_04pop_migLow.params <- fscRun(MSAT_04pop_migLow.params, num.sims = num_reps, exec = fscVersion)
MSAT_04pop_migHigh.params <- fscRun(MSAT_04pop_migHigh.params, num.sims = num_reps, exec = fscVersion)

# 16 POPULATIONS ----
# Write parameter files
MSAT_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16, 
                                     genetics = MSATgenetics, label = "MSAT_16pop_migLow", use.wd=TRUE)
MSAT_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16, 
                                      genetics = MSATgenetics, label = "MSAT_16pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_16pop_migLow.params <- fscRun(MSAT_16pop_migLow.params, num.sims = num_reps, exec = fscVersion)
MSAT_16pop_migHigh.params <- fscRun(MSAT_16pop_migHigh.params, num.sims = num_reps, exec = fscVersion)

# ---- DNA SIMULATIONS ----
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/DNA_marker/")
# 1 POPULATION ----
# Write parameter files. We make a mighHigh .par file as well, even though it's identical to migLow (with one population)
DNA_01pop_migLow.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes1, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)

# 4 POPULATIONS ----
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes4, migration = mig4Low, events = histEvent4, 
                                     genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes4, migration = mig4High, events = histEvent4, 
                                      genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)

# 16 POPULATIONS ----
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes16, migration = mig16Low, events = histEvent16, 
                                     genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes16, migration = mig16High, events = histEvent16, 
                                      genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)

# ---- CONVERT ARLEQUIN FILES TO GENIND ----
# Function converting Arlequin output to genind (through gtypes format)
strataG_arp2gen <- function(params, repNumber, marker){
  # Read in the Arlequin file, convert it to a gtype object, then to a genind object
  arp <- fscReadArp(params, sim=c(1,repNumber), marker = marker)
  gtype <- df2gtypes(arp, ploidy = 2)
  genind <- gtypes2genind(gtype)
  return(genind)
}

# TO DO: reformat this convertAllArp function to pass a vector of repNumbers to arp2gen, to use with lapply
# Function for converting ALL Arlequin files in a specified directory (using above arp2genind fxn)
# Outputs a list of genind objects, where each list item is a replicate
convertAllArp <- function(arp.path){
  # Retrieve original working directory, to reset to after conversion
  original.WD <- getwd()
  # Navigate to the folder containing simulation outputs
  setwd(arp.path)
  # Create an empty list object to receive list of genind
  genind.list <- list(length=length(dir()[str_detect(dir(), pattern = ".arp")]))
  # Convert all Arlequin files to adegenet, creating a list of genind objects
  for(i in 1:length(genind.list)){
    genind.list[[i]] <- strataG_arp2gen()
  }
  genind.list <- lapply(dir()[str_detect(dir(), pattern = ".arp")], strataG_arp2gen)
  # Reset to original working directory, and return a list of genind objects
  setwd(original.WD)
  return(genind.list)
}

# This is close, but the sim argument...we need to make a list, 
# A list similar to lapply(seq(1:5), function(x) c(x,1)), except items in reverse order...
sapply(MSAT_01pop_migLow.params, fscReadArp, sim=c(1,seq(1:5)))

# Convert microsatellite Arlequin outputs ----
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/MSAT_marker/")

MSAT_01pop_migLow_genind_1 <- arp2genind(MSAT_01pop_migLow.params, 1, marker = "microsat")
MSAT_01pop_migLow_genind_2 <- arp2genind(MSAT_01pop_migLow.params, 2, marker = "microsat")
MSAT_01pop_migLow_genind_3 <- arp2genind(MSAT_01pop_migLow.params, 3, marker = "microsat")
MSAT_01pop_migLow_genind_4 <- arp2genind(MSAT_01pop_migLow.params, 4, marker = "microsat")
MSAT_01pop_migLow_genind_5 <- arp2genind(MSAT_01pop_migLow.params, 5, marker = "microsat")
MSAT_01pop_migLow_genind <- list(MSAT_01pop_migLow_genind_1,MSAT_01pop_migLow_genind_2,
                                 MSAT_01pop_migLow_genind_3,MSAT_01pop_migLow_genind_4,
                                 MSAT_01pop_migLow_genind_5)

MSAT_01pop_migHigh_genind_1 <- arp2genind(MSAT_01pop_migHigh.params, 1, marker = "microsat")
MSAT_01pop_migHigh_genind_2 <- arp2genind(MSAT_01pop_migHigh.params, 2, marker = "microsat")
MSAT_01pop_migHigh_genind_3 <- arp2genind(MSAT_01pop_migHigh.params, 3, marker = "microsat")
MSAT_01pop_migHigh_genind_4 <- arp2genind(MSAT_01pop_migHigh.params, 4, marker = "microsat")
MSAT_01pop_migHigh_genind_5 <- arp2genind(MSAT_01pop_migHigh.params, 5, marker = "microsat")
MSAT_01pop_migHigh_genind <- list(MSAT_01pop_migHigh_genind_1,MSAT_01pop_migHigh_genind_2,
                                   MSAT_01pop_migHigh_genind_3,MSAT_01pop_migHigh_genind_4,
                                   MSAT_01pop_migHigh_genind_5)

MSAT_04pop_migLow_genind_1 <- arp2genind(MSAT_04pop_migLow.params, 1, marker = "microsat")
MSAT_04pop_migLow_genind_2 <- arp2genind(MSAT_04pop_migLow.params, 2, marker = "microsat")
MSAT_04pop_migLow_genind_3 <- arp2genind(MSAT_04pop_migLow.params, 3, marker = "microsat")
MSAT_04pop_migLow_genind_4 <- arp2genind(MSAT_04pop_migLow.params, 4, marker = "microsat")
MSAT_04pop_migLow_genind_5 <- arp2genind(MSAT_04pop_migLow.params, 5, marker = "microsat")
MSAT_04pop_migLow_genind <- list(MSAT_04pop_migLow_genind_1,MSAT_04pop_migLow_genind_2,
                                 MSAT_04pop_migLow_genind_3,MSAT_04pop_migLow_genind_4,
                                 MSAT_04pop_migLow_genind_5)

MSAT_04pop_migHigh_genind_1 <- arp2genind(MSAT_04pop_migHigh.params, 1, marker = "microsat")
MSAT_04pop_migHigh_genind_2 <- arp2genind(MSAT_04pop_migHigh.params, 2, marker = "microsat")
MSAT_04pop_migHigh_genind_3 <- arp2genind(MSAT_04pop_migHigh.params, 3, marker = "microsat")
MSAT_04pop_migHigh_genind_4 <- arp2genind(MSAT_04pop_migHigh.params, 4, marker = "microsat")
MSAT_04pop_migHigh_genind_5 <- arp2genind(MSAT_04pop_migHigh.params, 5, marker = "microsat")
MSAT_04pop_migHigh_genind <- list(MSAT_04pop_migHigh_genind_1,MSAT_04pop_migHigh_genind_2,
                                 MSAT_04pop_migHigh_genind_3,MSAT_04pop_migHigh_genind_4,
                                 MSAT_04pop_migHigh_genind_5)

MSAT_16pop_migLow_genind_1 <- arp2genind(MSAT_16pop_migLow.params, 1, marker = "microsat")
MSAT_16pop_migLow_genind_2 <- arp2genind(MSAT_16pop_migLow.params, 2, marker = "microsat")
MSAT_16pop_migLow_genind_3 <- arp2genind(MSAT_16pop_migLow.params, 3, marker = "microsat")
MSAT_16pop_migLow_genind_4 <- arp2genind(MSAT_16pop_migLow.params, 4, marker = "microsat")
MSAT_16pop_migLow_genind_5 <- arp2genind(MSAT_16pop_migLow.params, 5, marker = "microsat")
MSAT_16pop_migLow_genind <- list(MSAT_16pop_migLow_genind_1,MSAT_16pop_migLow_genind_2,
                                 MSAT_16pop_migLow_genind_3,MSAT_16pop_migLow_genind_4,
                                 MSAT_16pop_migLow_genind_5)

MSAT_16pop_migHigh_genind_1 <- arp2genind(MSAT_16pop_migHigh.params, 1, marker = "microsat")
MSAT_16pop_migHigh_genind_2 <- arp2genind(MSAT_16pop_migHigh.params, 2, marker = "microsat")
MSAT_16pop_migHigh_genind_3 <- arp2genind(MSAT_16pop_migHigh.params, 3, marker = "microsat")
MSAT_16pop_migHigh_genind_4 <- arp2genind(MSAT_16pop_migHigh.params, 4, marker = "microsat")
MSAT_16pop_migHigh_genind_5 <- arp2genind(MSAT_16pop_migHigh.params, 5, marker = "microsat")
MSAT_16pop_migHigh_genind <- list(MSAT_16pop_migHigh_genind_1,MSAT_16pop_migHigh_genind_2,
                                  MSAT_16pop_migHigh_genind_3,MSAT_16pop_migHigh_genind_4,
                                  MSAT_16pop_migHigh_genind_5)

# Convert DNA Arlequin outputs ----
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/DNA_marker/")

DNA_01pop_migLow_genind_1 <- arp2genind(DNA_01pop_migLow.params, 1, marker = "DNA")
DNA_01pop_migLow_genind_2 <- arp2genind(DNA_01pop_migLow.params, 2, marker = "DNA")
DNA_01pop_migLow_genind_3 <- arp2genind(DNA_01pop_migLow.params, 3, marker = "DNA")
DNA_01pop_migLow_genind_4 <- arp2genind(DNA_01pop_migLow.params, 4, marker = "DNA")
DNA_01pop_migLow_genind_5 <- arp2genind(DNA_01pop_migLow.params, 5, marker = "DNA")
DNA_01pop_migLow_genind <- list(DNA_01pop_migLow_genind_1, DNA_01pop_migLow_genind_2,
                                DNA_01pop_migLow_genind_3, DNA_01pop_migLow_genind_4,
                                DNA_01pop_migLow_genind_5)

DNA_01pop_migHigh_genind_1 <- arp2genind(DNA_01pop_migHigh.params, 1, marker = "DNA")
DNA_01pop_migHigh_genind_2 <- arp2genind(DNA_01pop_migHigh.params, 2, marker = "DNA")
DNA_01pop_migHigh_genind_3 <- arp2genind(DNA_01pop_migHigh.params, 3, marker = "DNA")
DNA_01pop_migHigh_genind_4 <- arp2genind(DNA_01pop_migHigh.params, 4, marker = "DNA")
DNA_01pop_migHigh_genind_5 <- arp2genind(DNA_01pop_migHigh.params, 5, marker = "DNA")
DNA_01pop_migHigh_genind <- list(DNA_01pop_migHigh_genind_1, DNA_01pop_migHigh_genind_2,
                                DNA_01pop_migHigh_genind_3, DNA_01pop_migHigh_genind_4,
                                DNA_01pop_migHigh_genind_5)

DNA_04pop_migLow_genind_1 <- arp2genind(DNA_04pop_migLow.params, 1, marker = "DNA")
DNA_04pop_migLow_genind_2 <- arp2genind(DNA_04pop_migLow.params, 2, marker = "DNA")
DNA_04pop_migLow_genind_3 <- arp2genind(DNA_04pop_migLow.params, 3, marker = "DNA")
DNA_04pop_migLow_genind_4 <- arp2genind(DNA_04pop_migLow.params, 4, marker = "DNA")
DNA_04pop_migLow_genind_5 <- arp2genind(DNA_04pop_migLow.params, 5, marker = "DNA")
DNA_04pop_migLow_genind <- list(DNA_04pop_migLow_genind_1, DNA_04pop_migLow_genind_2,
                                DNA_04pop_migLow_genind_3, DNA_04pop_migLow_genind_4,
                                DNA_04pop_migLow_genind_5)

DNA_04pop_migHigh_genind_1 <- arp2genind(DNA_04pop_migHigh.params, 1, marker = "DNA")
DNA_04pop_migHigh_genind_2 <- arp2genind(DNA_04pop_migHigh.params, 2, marker = "DNA")
DNA_04pop_migHigh_genind_3 <- arp2genind(DNA_04pop_migHigh.params, 3, marker = "DNA")
DNA_04pop_migHigh_genind_4 <- arp2genind(DNA_04pop_migHigh.params, 4, marker = "DNA")
DNA_04pop_migHigh_genind_5 <- arp2genind(DNA_04pop_migHigh.params, 5, marker = "DNA")
DNA_04pop_migHigh_genind <- list(DNA_04pop_migHigh_genind_1, DNA_04pop_migHigh_genind_2,
                                 DNA_04pop_migHigh_genind_3, DNA_04pop_migHigh_genind_4,
                                 DNA_04pop_migHigh_genind_5)

DNA_16pop_migLow_genind_1 <- arp2genind(DNA_16pop_migLow.params, 1, marker = "DNA")
DNA_16pop_migLow_genind_2 <- arp2genind(DNA_16pop_migLow.params, 2, marker = "DNA")
DNA_16pop_migLow_genind_3 <- arp2genind(DNA_16pop_migLow.params, 3, marker = "DNA")
DNA_16pop_migLow_genind_4 <- arp2genind(DNA_16pop_migLow.params, 4, marker = "DNA")
DNA_16pop_migLow_genind_5 <- arp2genind(DNA_16pop_migLow.params, 5, marker = "DNA")
DNA_16pop_migLow_genind <- list(DNA_16pop_migLow_genind_1, DNA_16pop_migLow_genind_2,
                                DNA_16pop_migLow_genind_3, DNA_16pop_migLow_genind_4,
                                DNA_16pop_migLow_genind_5)

DNA_16pop_migHigh_genind_1 <- arp2genind(DNA_16pop_migHigh.params, 1, marker = "DNA")
DNA_16pop_migHigh_genind_2 <- arp2genind(DNA_16pop_migHigh.params, 2, marker = "DNA")
DNA_16pop_migHigh_genind_3 <- arp2genind(DNA_16pop_migHigh.params, 3, marker = "DNA")
DNA_16pop_migHigh_genind_4 <- arp2genind(DNA_16pop_migHigh.params, 4, marker = "DNA")
DNA_16pop_migHigh_genind_5 <- arp2genind(DNA_16pop_migHigh.params, 5, marker = "DNA")
DNA_16pop_migHigh_genind <- list(DNA_16pop_migHigh_genind_1, DNA_16pop_migHigh_genind_2,
                                 DNA_16pop_migHigh_genind_3, DNA_16pop_migHigh_genind_4,
                                 DNA_16pop_migHigh_genind_5)

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

# 2. LOWER FST FOR SCENARIO WITH LOWER MIGRATION RATES ----
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
# QUESTION: when we calculate allele frequencies, do we divide my the number of individuals in the
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
