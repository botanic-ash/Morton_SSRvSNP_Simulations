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
# 1 POPULATION
demeA <- fscDeme(deme.size = 1200, sample.size = 1200)
demes <- fscSettingsDemes(demeA)
# Specify genetic markers
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
genetics <- fscSettingsGenetics(msats, num.chrom = 20)
# Write the parameters to a .par file (in the pwd)
MSAT_01pop_migLow.params <- fscWrite(demes = demes, genetics = genetics, label = "MSAT_01pop_migLow", use.wd=TRUE)
print(str(MSAT_01pop_migLow.params))
# Run the parameters. Running creates an updated params object, hence attaching the output to the params variable
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps)

# 4 POPULATIONS
# Demes
demeB_1 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_2 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_3 <- fscDeme(deme.size = 300, sample.size = 300)
demeB_4 <- fscDeme(deme.size = 300, sample.size = 300)
demes <- fscSettingsDemes(demeB_1, demeB_2, demeB_3, demeB_4)
# Migration
mig.mat.0 <- matrix(low_mig, nrow=4, ncol = 4); diag(mig.mat) <- 0
mig.mat.1 <- matrix(0, nrow=4, ncol = 4)
mig <- fscSettingsMigration(mig.mat.0, mig.mat.1)
# Historical events
hist.event0 <- fscEvent(event.time = 50000, source = 0, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event1 <- fscEvent(event.time = 50000, source = 1, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event2 <- fscEvent(event.time = 50000, source = 2, sink = 0, prop.migrants = 0, migr.mat = 1)
hist.event3 <- fscEvent(event.time = 50000, source = 3, sink = 0, prop.migrants = 0, migr.mat = 1)
histEvent <- fscSettingsEvents(hist.event0, hist.event1, hist.event2, hist.event3)
# Genetics
msats <- fscBlock_microsat(num.loci = 1, mut.rate = 5e-4, range.constraint = 10)
genetics <- fscSettingsGenetics(msats, num.chrom = 20)
# Write parameters file
MSAT_04pop_migLow.params <- fscWrite(demes = demes, migration = mig, events = histEvent, 
                                     genetics = genetics, label = "MSAT_04pop_migLow", use.wd=TRUE)
# Run parameters file
MSAT_01pop_migLow.params <- fscRun(MSAT_01pop_migLow.params, num.sims = num_reps)
