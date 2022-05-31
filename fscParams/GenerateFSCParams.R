# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FSC PARAMETERS USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(strataG)
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/")

# ----VARIABLES----
num_reps <- 5
# Demes
# 1 Population
demeA <- fscDeme(deme.size = 1200, sample.size = 1200)
# 4 Populations
demeB_1 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_2 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_3 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_4 <- fscDeme(deme.size = 300, sample.size = 300)
# 16 Populations
demeC_1 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_2 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_3 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_4 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_5 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_6 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_7 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_8 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_9 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_10 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_11 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_12 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_13 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_14 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_15 <- fscDeme(deme.size = 75, sample.size = 75)
demeC_16 <- fscDeme(deme.size = 75, sample.size = 75)
# Migration
low_mig <- 0.001
high_mig <- 0.01
mig.mat.4.Low <- matrix(low_mig, nrow=4, ncol = 4); diag(mig.mat.4.Low) <- 0
mig.mat.4.High <- matrix(high_mig, nrow=4, ncol = 4); diag(mig.mat.4.High) <- 0
mig.mat.4.Final <- matrix(0, nrow=4, ncol = 4)
mig.mat.16.Low <- matrix(low_mig, nrow=16, ncol = 16); diag(mig.mat.16.Low) <- 0
mig.mat.16.High <- matrix(high_mig, nrow=16, ncol = 16); diag(mig.mat.16.High) <- 0
mig.mat.16.Final <- matrix(0, nrow=16, ncol = 16)
# Historical events
hist.event0 <- fscEvent(event.time = 50000, source = 0, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event1 <- fscEvent(event.time = 50000, source = 1, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event2 <- fscEvent(event.time = 50000, source = 2, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event3 <- fscEvent(event.time = 50000, source = 3, sink = 0, prop.migrants = 1, migr.mat = 1)
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
# MSAT Genetic parameters
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)
# DNA Genetic parameters
dna <- fscBlock_dna(sequence.length = 25, mut.rate = 1e-5, chromosome = 4)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, dna, num.chrom = 5)

# ----MICROSATELLITE SIMULATIONS----
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/params_MSAT")
# 1 POPULATION----
demes <- fscSettingsDemes(demeA)
# Write parameter files. We make a mighHigh .par file as well, even though it's identical to migLow (with one population)
MSAT_01pop_migLow.params <- fscWrite(demes = demes, genetics = MSATgenetics, label = "MSAT_01pop_migLow", use.wd=TRUE)
MSAT_01pop_migHigh.params <- fscWrite(demes = demes, genetics = MSATgenetics, label = "MSAT_01pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps, exec = "fsc2702")
MSAT_01pop_migHigh.params <- fscRun(MSAT_01pop_migHigh.params, num.sims = num_reps, exec = "fsc2702")

# 4 POPULATIONS----
# Demes
demes <- fscSettingsDemes(demeB_1, demeB_2, demeB_3, demeB_4)
# Migration
mig4Low <- fscSettingsMigration(mig.mat.4.Low, mig.mat.4.Final)
mig4High <- fscSettingsMigration(mig.mat.4.High, mig.mat.4.Final)
# Historical events
histEvent <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
# Write parameter files
MSAT_04pop_migLow.params <- fscWrite(demes = demes, migration = mig4Low, events = histEvent, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
MSAT_04pop_migHigh.params <- fscWrite(demes = demes, migration = mig4High, events = histEvent, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_04pop_migLow.params <- fscRun(MSAT_04pop_migLow.params, num.sims = num_reps, exec = "fsc2702")
MSAT_04pop_migHigh.params <- fscRun(MSAT_04pop_migHigh.params, num.sims = num_reps, exec = "fsc2702")

# 16 POPULATIONS----
# Demes
demes <- fscSettingsDemes(demeC_1,demeC_2,demeC_3,demeC_4,demeC_5,demeC_6,demeC_7,demeC_8,demeC_9,demeC_10,demeC_11,demeC_12,
                          demeC_13,demeC_14,demeC_15,demeC_16)
# Migration
mig16Low <- fscSettingsMigration(mig.mat.16.Low, mig.mat.16.Final)
mig16High <- fscSettingsMigration(mig.mat.16.High, mig.mat.16.Final)
# Historical events
histEvent <- fscSettingsEvents(hist.event0,hist.event1,hist.event2,hist.event3,hist.event4,hist.event5,hist.event6,hist.event7,
                               hist.event8,hist.event9,hist.event10,hist.event11,hist.event12,hist.event13,hist.event14,hist.event15)
# Write parameter files
MSAT_16pop_migLow.params <- fscWrite(demes = demes, migration = mig16Low, events = histEvent, 
                                     genetics = MSATgenetics, label = "MSAT_16pop_migLow", use.wd=TRUE)
MSAT_16pop_migHigh.params <- fscWrite(demes = demes, migration = mig16High, events = histEvent, 
                                      genetics = MSATgenetics, label = "MSAT_16pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_16pop_migLow.params <- fscRun(MSAT_16pop_migLow.params, num.sims = num_reps, exec = "fsc2702")
MSAT_16pop_migHigh.params <- fscRun(MSAT_16pop_migHigh.params, num.sims = num_reps, exec = "fsc2702")

# ----DNA SIMULATIONS----
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/params_DNA")
# 1 POPULATION----
demes <- fscSettingsDemes(demeA)
# Write parameter files. We make a mighHigh .par file as well, even though it's identical to migLow (with one population)
DNA_01pop_migLow.params <- fscWrite(demes = demes, genetics = DNAgenetics, label = "DNA_01pop_migLow", use.wd=TRUE)
DNA_01pop_migHigh.params <- fscWrite(demes = demes, genetics = DNAgenetics, label = "DNA_01pop_migHigh", use.wd=TRUE)
# Run parameter files
DNA_01pop_migLow.params <- fscRun(DNA_01pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = "fsc2702")
DNA_01pop_migHigh.params <- fscRun(DNA_01pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = "fsc2702")

DNA_01pop_migLow <- fscReadArp(DNA_01pop_migLow.params)

DNA_01pop_migLow_genind <- fsc2gtypes(DNA_01pop_migLow.params,marker = "dna")
DNA_01pop_migLow_genind <- fsc2gtypes(DNA_01pop_migLow,marker = "dna")



test <- arlequinRead("DNA_01pop_migHigh/DNA_01pop_migHigh_1_3.arp")
DNA_01pop_migLow_genind <- arp2gtypes(test)
DNA_01pop_migLow_genind <- gtypes2genind(DNA_01pop_migLow_genind)
nLoc(DNA_01pop_migLow_genind)

# 4 POPULATIONS----
# Demes
demes <- fscSettingsDemes(demeB_1, demeB_2, demeB_3, demeB_4)
# Migration
mig4Low <- fscSettingsMigration(mig.mat.4.Low, mig.mat.4.Final)
mig4High <- fscSettingsMigration(mig.mat.4.High, mig.mat.4.Final)
# Historical events
histEvent <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
# Write parameter files
DNA_04pop_migLow.params <- fscWrite(demes = demes, migration = mig4Low, events = histEvent, 
                                     genetics = DNAgenetics, label = "DNA_04pop_migLow", use.wd=TRUE)
DNA_04pop_migHigh.params <- fscWrite(demes = demes, migration = mig4High, events = histEvent, 
                                      genetics = DNAgenetics, label = "DNA_04pop_migHigh", use.wd=TRUE)
# Run parameter files
DNA_04pop_migLow.params <- fscRun(DNA_04pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = "fsc2702")
DNA_04pop_migHigh.params <- fscRun(DNA_04pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = "fsc2702")

# 16 POPULATIONS----
# Demes
demes <- fscSettingsDemes(demeC_1,demeC_2,demeC_3,demeC_4,demeC_5,demeC_6,demeC_7,demeC_8,demeC_9,demeC_10,demeC_11,demeC_12,
                          demeC_13,demeC_14,demeC_15,demeC_16)
# Migration
mig16Low <- fscSettingsMigration(mig.mat.16.Low, mig.mat.16.Final)
mig16High <- fscSettingsMigration(mig.mat.16.High, mig.mat.16.Final)
# Historical events
histEvent <- fscSettingsEvents(hist.event0,hist.event1,hist.event2,hist.event3,hist.event4,hist.event5,hist.event6,hist.event7,
                               hist.event8,hist.event9,hist.event10,hist.event11,hist.event12,hist.event13,hist.event14,hist.event15)
# Write parameter files
DNA_16pop_migLow.params <- fscWrite(demes = demes, migration = mig16Low, events = histEvent, 
                                     genetics = DNAgenetics, label = "DNA_16pop_migLow", use.wd=TRUE)
DNA_16pop_migHigh.params <- fscWrite(demes = demes, migration = mig16High, events = histEvent, 
                                      genetics = DNAgenetics, label = "DNA_16pop_migHigh", use.wd=TRUE)
# Run parameter files
DNA_16pop_migLow.params <- fscRun(DNA_16pop_migLow.params, num.sims = num_reps, all.sites = TRUE, exec = "fsc2702")
DNA_16pop_migHigh.params <- fscRun(DNA_16pop_migHigh.params, num.sims = num_reps, all.sites = TRUE, exec = "fsc2702")
