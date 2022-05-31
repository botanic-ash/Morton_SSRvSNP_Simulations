# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FSC PARAMETERS USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

library(strataG)
setwd("~/Documents/SSRvSNP/Simulations/Code/strataG_Demo/")

# ----VARIABLES----
num_reps <- 5
low_mig <- 0.001
high_mig <- 0.01
# MSATS
# Genetic markers: same for all MSAT datasets
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
MSATgenetics <- fscSettingsGenetics(msats, num.chrom = 20)
# DNA

# ----EXAMPLE----
# Generate demes
deme0 <- fscDeme(deme.size = 1000, sample.size = 4)
demes <- fscSettingsDemes(deme0)
# Specify genetic markers
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 1e-3)
genetics <- fscSettingsGenetics(msats, num.chrom = 5)
# Write the parameters to a .par file (in the pwd)
ex1.params <- fscWrite(demes = demes, genetics = genetics, label = "ex1", use.wd=TRUE)
print(str(ex1.params))
# Run the parameters. Running creates an updated params object, hence attaching the output to the params variable
ex1.params <- fscRun(ex1.params, num.sims = 2)

# ----MICROSATELLITE SIMULATIONS----
# 1 POPULATION----
demeA <- fscDeme(deme.size = 1200, sample.size = 1200)
demes <- fscSettingsDemes(demeA)
# Write parameter files. We make a mighHigh .par file as well, even though it's identical to migLow (with one population)
MSAT_01pop_migLow.params <- fscWrite(demes = demes, genetics = MSATgenetics, label = "MSAT_01pop_migLow", use.wd=TRUE)
MSAT_01pop_migHigh.params <- fscWrite(demes = demes, genetics = MSATgenetics, label = "MSAT_01pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps)
MSAT_01pop_migHigh.params <- fscRun(MSAT_01pop_migHigh.params, num.sims = num_reps)

# 4 POPULATIONS----
# Demes
demeB_1 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_2 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_3 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_4 <- fscDeme(deme.size = 300, sample.size = 300)
demes <- fscSettingsDemes(demeB_1, demeB_2, demeB_3, demeB_4)
# Migration
mig.mat.Low <- matrix(low_mig, nrow=4, ncol = 4); diag(mig.mat.Low) <- 0
mig.mat.High <- matrix(high_mig, nrow=4, ncol = 4); diag(mig.mat.High) <- 0
mig.mat.1 <- matrix(0, nrow=4, ncol = 4)
migLow <- fscSettingsMigration(mig.mat.Low, mig.mat.1)
migHigh <- fscSettingsMigration(mig.mat.High, mig.mat.1)
# Historical events
hist.event0 <- fscEvent(event.time = 50000, source = 0, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event1 <- fscEvent(event.time = 50000, source = 1, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event2 <- fscEvent(event.time = 50000, source = 2, sink = 0, prop.migrants = 1, migr.mat = 1)
hist.event3 <- fscEvent(event.time = 50000, source = 3, sink = 0, prop.migrants = 1, migr.mat = 1)
histEvent <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
# Write parameter files
MSAT_04pop_migLow.params <- fscWrite(demes = demes, migration = migLow, events = histEvent, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
MSAT_04pop_migHigh.params <- fscWrite(demes = demes, migration = migHigh, events = histEvent, 
                                     genetics = MSATgenetics, label = "MSAT_04pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_04pop_migLow.params <- fscRun(MSAT_04pop_migLow.params, num.sims = num_reps)
MSAT_04pop_migHigh.params <- fscRun(MSAT_04pop_migHigh.params, num.sims = num_reps)

# 16 POPULATIONS----
# Demes
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
demes <- fscSettingsDemes(demeC_1,demeC_2,demeC_3,demeC_4,demeC_5,demeC_6,demeC_7,demeC_8,demeC_9,demeC_10,demeC_11,demeC_12,
                          demeC_13,demeC_14,demeC_15,demeC_16)
# Migration
mig.mat.Low <- matrix(low_mig, nrow=16, ncol = 16); diag(mig.mat.Low) <- 0
mig.mat.High <- matrix(high_mig, nrow=16, ncol = 16); diag(mig.mat.High) <- 0
mig.mat.1 <- matrix(0, nrow=16, ncol = 16)
migLow <- fscSettingsMigration(mig.mat.Low, mig.mat.1)
migHigh <- fscSettingsMigration(mig.mat.High, mig.mat.1)
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
histEvent <- fscSettingsEvents(hist.event0,hist.event1,hist.event2,hist.event3,hist.event4,hist.event5,hist.event6,hist.event7,
                               hist.event8,hist.event9,hist.event10,hist.event11,hist.event12,hist.event13,hist.event14,hist.event15)
# Write parameter files
MSAT_16pop_migLow.params <- fscWrite(demes = demes, migration = migLow, events = histEvent, 
                                     genetics = MSATgenetics, label = "MSAT_16pop_migLow", use.wd=TRUE)
MSAT_16pop_migHigh.params <- fscWrite(demes = demes, migration = migHigh, events = histEvent, 
                                      genetics = MSATgenetics, label = "MSAT_16pop_migHigh", use.wd=TRUE)
# Run parameter files
MSAT_16pop_migLow.params <- fscRun(MSAT_16pop_migLow.params, num.sims = num_reps)
MSAT_16pop_migHigh.params <- fscRun(MSAT_16pop_migHigh.params, num.sims = num_reps)
