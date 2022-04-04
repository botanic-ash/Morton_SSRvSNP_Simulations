# Parsing fastSimcoal2 outputs

library(adegenet)
library(hierfstat)
library(stringr)
setwd("/home/akoontz/Documents/SSRvSNP/Simulations/Code/")

# Specify /RAID1 directory, where outputs are being held
sim.directory <- "/RAID1/Simulations/Even_20220329_FSCdraft/"
# Additional directory, for samples simulated using the -s flag (for DNA markers)
# sim.directory <- "/RAID1/Simulations/SNP_20220404_FSCdraft/"

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
# Exploring what DNA alleles matrices look like
str(DNA_1pop_migLow[[1]]@tab)
DNA_1pop_migLow[[1]]@tab[1:5,1:5]

str(DNA_4pop_migLow[[1]]@tab)
DNA_4pop_migLow[[1]]@tab[1:5,1:5]

unique(rowSums(DNA_1pop_migLow[[1]]@tab))
length(which(DNA_1pop_migLow[[1]]@tab[1,] == 0))
length(which(DNA_1pop_migLow[[1]]@tab[2,] == 0))

length(which(DNA_1pop_migLow[[1]]@tab[,1] == 0))
length(which(DNA_1pop_migLow[[1]]@tab[,2] == 0))

# Microsatellites
# Number of alleles for each simulation instance
lapply(MSAT_1pop_migLow, FUN = function(x) ncol(x@tab))
# Average
mean(as.numeric(lapply(MSAT_1pop_migLow, FUN = function(x) ncol(x@tab))))
# 4 populations
lapply(MSAT_4pop_migLow, FUN = function(x) ncol(x@tab))
mean(as.numeric(lapply(MSAT_4pop_migLow, FUN = function(x) ncol(x@tab))))
# 16 populations
lapply(MSAT_16pop_migLow, FUN = function(x) ncol(x@tab))
mean(as.numeric(lapply(MSAT_16pop_migLow, FUN = function(x) ncol(x@tab))))

mean(as.numeric(lapply(MSAT_1pop_migHigh, FUN = function(x) ncol(x@tab))))
mean(as.numeric(lapply(MSAT_4pop_migHigh, FUN = function(x) ncol(x@tab))))
mean(as.numeric(lapply(MSAT_16pop_migHigh, FUN = function(x) ncol(x@tab))))

# DNA
# Number of alleles for each simulation instance
lapply(DNA_1pop_migLow, FUN = function(x) ncol(x@tab))
lapply(DNA_4pop_migLow, FUN = function(x) ncol(x@tab))
lapply(DNA_16pop_migLow, FUN = function(x) ncol(x@tab))
lapply(DNA_16pop_migHigh, FUN = function(x) ncol(x@tab))

# 2. Lower Fst for scenarios with lower migration rates

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