# Parsing fastSimcoal2 outputs

library(adegenet)
library(hierfstat)
library(stringr)
setwd("/home/akoontz/Documents/SSRvSNP/Simulations/Code/")

# Specify /RAID1 directory, where outputs are being held
sim.directory <- "/RAID1/Simulations/draft_fscParams/"

# sim.directory <- "/RAID1/Simulations/draft_fscParams_lowMutation/"
# sim.directory <- "/RAID1/Simulations/draft_fscParams_highMutation/"
sim.directory <- "/RAID1/Simulations/draft_fscParams_DNA100/"

# Conversion using strataG package----
library(strataG)
# Wrapper function for converting Arlequin (.arp) to genepop (.gen), using strataG functions
strataG_arp2gen <- function(arp.path){
  arp.object <- arlequinRead(arp.path)
  gtype.object <- arp2gtypes(arp.object)
  genind.object <- gtypes2genind(gtype.object)
  return(genind.object)
}

# Function for converting all Arlequin files in a specified directory (using above strataG_arp2gen fxn)
convertAllArp <- function(arp.path){
  # Retrieve original working directory, to reset to after conversion
  original.WD <- getwd()
  # Navigate to the folder containing simulation outputs
  setwd(arp.path)
  # Create an empty list object to receive list of genind
  genind.list <- list(length=length(dir()[str_detect(dir(), pattern = ".arp")]))
  # Convert all Arlequin files to adegenet, creating a list of genind objects
  genind.list <- lapply(dir()[str_detect(dir(), pattern = ".arp")], strataG_arp2gen)
  # Reset to original working directory, and return a list of genind objects
  setwd(original.WD)
  return(genind.list)
}
# For all simulations, number of individuals = 1200
nInd <- 1200

# Microsatellite files
MSAT_1pop_migLow <- convertAllArp(paste0(sim.directory,"MSAT_01pops_migLow"))
MSAT_1pop_migHigh <- convertAllArp(paste0(sim.directory,"MSAT_01pops_migHigh"))

MSAT_4pop_migLow <- convertAllArp(paste0(sim.directory,"MSAT_04pops_migLow"))
MSAT_4pop_migHigh <- convertAllArp(paste0(sim.directory,"MSAT_04pops_migHigh"))

MSAT_16pop_migLow <- convertAllArp(paste0(sim.directory,"MSAT_16pops_migLow"))
MSAT_16pop_migHigh <- convertAllArp(paste0(sim.directory,"MSAT_16pops_migHigh"))

# DNA files
DNA_1pop_migLow <- convertAllArp(paste0(sim.directory,"DNA_01pops_migLow"))
DNA_1pop_migHigh <- convertAllArp(paste0(sim.directory,"DNA_01pops_migHigh"))

DNA_4pop_migLow <- convertAllArp(paste0(sim.directory,"DNA_04pops_migLow"))
DNA_4pop_migHigh <- convertAllArp(paste0(sim.directory,"DNA_04pops_migHigh"))

DNA_16pop_migLow <- convertAllArp(paste0(sim.directory,"DNA_16pops_migLow"))
DNA_16pop_migHigh <- convertAllArp(paste0(sim.directory,"DNA_16pops_migHigh"))

# Sense check----
# 1. More alleles in scenarios with a larger number of populations
# MSAT
mean(sapply(MSAT_1pop_migLow, function(x) ncol(x@tab)))
mean(sapply(MSAT_4pop_migLow, function(x) ncol(x@tab)))
mean(sapply(MSAT_16pop_migLow, function(x) ncol(x@tab)))

mean(sapply(MSAT_1pop_migHigh, function(x) ncol(x@tab)))
mean(sapply(MSAT_4pop_migHigh, function(x) ncol(x@tab)))
mean(sapply(MSAT_16pop_migHigh, function(x) ncol(x@tab)))

# DNA
mean(sapply(DNA_1pop_migLow, function(x) ncol(x@tab)))
mean(sapply(DNA_4pop_migLow, function(x) ncol(x@tab)))
mean(sapply(DNA_16pop_migLow, function(x) ncol(x@tab)))

mean(sapply(DNA_1pop_migHigh, function(x) ncol(x@tab)))
mean(sapply(DNA_4pop_migHigh, function(x) ncol(x@tab)))
mean(sapply(DNA_16pop_migHigh, function(x) ncol(x@tab)))

# 2. Lower Fst for scenarios with lower migration rates
# MSAT
sapply(MSAT_4pop_migLow, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_4pop_migHigh, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

sapply(MSAT_16pop_migLow, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(MSAT_16pop_migHigh, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# DNA
sapply(DNA_4pop_migLow, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))
sapply(DNA_4pop_migHigh, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

sapply(DNA_16pop_migLow, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE)) # Error
sapply(DNA_16pop_migHigh, function(x) mean(c(pairwise.neifst(genind2hierfstat(x))), na.rm=TRUE))

# 3. Allele frequency spectra
# MSAT
# DNA

# 1 population
unique(colSums(DNA_1pop_migLow[[1]]@tab))
nInd(DNA_1pop_migLow[[1]])

DNA_allFreq_1popLow_1 <- colSums(DNA_1pop_migLow[[1]]@tab)/(nInd*2)*100
DNA_allFreq_1popLow_2 <- colSums(DNA_1pop_migLow[[2]]@tab)/(nInd*2)*100
DNA_allFreq_1popLow_3 <- colSums(DNA_1pop_migLow[[3]]@tab)/(nInd*2)*100
DNA_allFreq_1popLow_4 <- colSums(DNA_1pop_migLow[[4]]@tab)/(nInd*2)*100
DNA_allFreq_1popLow_5 <- colSums(DNA_1pop_migLow[[5]]@tab)/(nInd*2)*100

unique(DNA_allFreq_1popLow_1)
unique(DNA_allFreq_1popLow_2)
unique(DNA_allFreq_1popLow_3)
unique(DNA_allFreq_1popLow_4)
unique(DNA_allFreq_1popLow_5)

hist(DNA_allFreq_1popLow_1)
hist(DNA_allFreq_1popLow_2)
hist(DNA_allFreq_1popLow_3)
hist(DNA_allFreq_1popLow_4)
hist(DNA_allFreq_1popLow_5)

# 16 populations
unique(colSums(DNA_16pop_migLow[[1]]@tab))
nInd(DNA_16pop_migLow[[1]])

DNA_allFreq_16popLow_1 <- colSums(DNA_16pop_migLow[[1]]@tab)/(nInd*2)*100
DNA_allFreq_16popLow_2 <- colSums(DNA_16pop_migLow[[2]]@tab)/(nInd*2)*100
DNA_allFreq_16popLow_3 <- colSums(DNA_16pop_migLow[[3]]@tab)/(nInd*2)*100
DNA_allFreq_16popLow_4 <- colSums(DNA_16pop_migLow[[4]]@tab)/(nInd*2)*100
DNA_allFreq_16popLow_5 <- colSums(DNA_16pop_migLow[[5]]@tab)/(nInd*2)*100

unique(DNA_allFreq_16popLow_1)
unique(DNA_allFreq_16popLow_2)
unique(DNA_allFreq_16popLow_3)
unique(DNA_allFreq_16popLow_4)
unique(DNA_allFreq_16popLow_5)

hist(DNA_allFreq_16popLow_1)
hist(DNA_allFreq_16popLow_2)
hist(DNA_allFreq_16popLow_3)
hist(DNA_allFreq_16popLow_4)
hist(DNA_allFreq_16popLow_5)

# # OUTDATED: Conversion using diveRsity package----
# library(diveRsity)
# 
# # DNA Markers
# arp2gen("DNA_01pops_migLow/DNA_01pops_migLow_1_1.arp")
# # Error: "Data are not in 'MICROSAT' format!"
# 
# # MSAT Markers
# # 1 population 
# arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_1.arp")
# # Error: "Error in substr(dat[sampSizeLine], start = sampNpos + 1, stop = nchar(dat[sampSizeLine])): object 'sampNpos' not found"
# # 4 populations
# arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_1.arp")
# MSAT_1pop_migHigh_diveRsity <- read.genepop("MSAT_04pops_migLow/MSAT_04pops_migLow_1_1.gen", ncode=3)
# nInd(MSAT_1pop_migHigh_diveRsity)