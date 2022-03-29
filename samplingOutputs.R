# Parsing fastSimcoal2 Arlequin outputs

library(adegenet)
library(hierfstat)
# Specify /RAID1 directory, where outputs are being held
setwd("/RAID1/Simulations/Even_20220329_FSCdraft/")

# Conversion using strataG package----
library(strataG)
# Wrapper function for converting Arlequin (.arp) to genepop (.gen), using strataG functions
strataG_arp2gen <- function(arp.path){
  arp.object <- arlequinRead(arp.path)
  gtype.object <- arp2gtypes(arp.object)
  genind.object <- gtypes2genind(gtype.object)
  return(genind.object)
}

# Microsatellite files
MSAT_1pop_migHigh <- strataG_arp2gen("MSAT_01pops_migHigh/MSAT_01pops_migHigh_1_1.arp")

MSAT_4pops_migHigh <- strataG_arp2gen("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_1.arp")
nInd(MSAT_4pops_migHigh)

MSAT_16pops_migHigh <- strataG_arp2gen("MSAT_16pops_migHigh/MSAT_16pops_migHigh_1_1.arp")
nInd(MSAT_16pops_migHigh)

# DNA files
DNA_1pop_migHigh <- strataG_arp2gen("DNA_01pops_migHigh/DNA_01pops_migHigh_1_1.arp")
nInd(DNA_1pop_migHigh)

DNA_4pops_migHigh <- strataG_arp2gen("DNA_04pops_migHigh/DNA_04pops_migHigh_1_1.arp")
nInd(DNA_4pops_migHigh)

DNA_16pops_migHigh <- strataG_arp2gen("DNA_16pops_migHigh/DNA_16pops_migHigh_1_1.arp")
nInd(DNA_16pops_migHigh)

# OUTDATED: Conversion using diveRsity package----
library(diveRsity)

# DNA Markers
arp2gen("DNA_01pops_migLow/DNA_01pops_migLow_1_1.arp")
# Error: "Data are not in 'MICROSAT' format!"

# MSAT Markers
# 1 population 
arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_1.arp")
# Error: "Error in substr(dat[sampSizeLine], start = sampNpos + 1, stop = nchar(dat[sampSizeLine])): object 'sampNpos' not found"
# 4 populations
arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_1.arp")
MSAT_1pop_migHigh_diveRsity <- read.genepop("MSAT_04pops_migLow/MSAT_04pops_migLow_1_1.gen", ncode=3)
nInd(MSAT_1pop_migHigh_diveRsity)