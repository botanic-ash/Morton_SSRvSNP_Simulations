# Parsing fastSimcoal2 Arlequin outputs
library(adegenet)
library(hierfstat)
# Specify /RAID1 directory, where outputs are being held
setwd("/RAID1/Simulations/Original_fsc_draftParams/")

# Conversion using diveRsity package----
library(diveRsity)

# MSAT examples: conversions only need to be run once
# 1 population----
# Not working. Error: "Data are not in 'MICROSAT' format!". This is somehow caused by having 1 population?
# Low migration
arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_1.arp")
arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_2.arp")
arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_3.arp")
arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_4.arp")
arp2gen("MSAT_01pops_migLow/MSAT_01pops_migLow_1_5.arp")
# High migration
arp2gen("MSAT_01pops_migHigh/MSAT_01pops_migHigh_1_1.arp")
arp2gen("MSAT_01pops_migHigh/MSAT_01pops_migHigh_1_2.arp")
arp2gen("MSAT_01pops_migHigh/MSAT_01pops_migHigh_1_3.arp")
arp2gen("MSAT_01pops_migHigh/MSAT_01pops_migHigh_1_4.arp")
arp2gen("MSAT_01pops_migHigh/MSAT_01pops_migHigh_1_5.arp")

# 4 populations----
# Low migration
arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_1.arp")
arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_2.arp")
arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_3.arp")
arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_4.arp")
arp2gen("MSAT_04pops_migLow/MSAT_04pops_migLow_1_5.arp")
# High migration
arp2gen("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_1.arp")
arp2gen("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_2.arp")
arp2gen("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_3.arp")
arp2gen("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_4.arp")
arp2gen("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_5.arp")

# 04pops conversions generate the below warnings
# Warning messages:
# 1: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  :
#  data length [10593] is not a sub-multiple or multiple of the number of rows [625]
# 2: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] + 1)], split = "\\s+")),  :
#  data length [10578] is not a sub-multiple or multiple of the number of rows [625]
# 3: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  :
#  data length [10593] is not a sub-multiple or multiple of the number of rows [625]
# 4: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] + 1)], split = "\\s+")),  :
#  data length [10578] is not a sub-multiple or multiple of the number of rows [625]
# 5: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  :
#  data length [10593] is not a sub-multiple or multiple of the number of rows [625]
# 6: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] + 1)], split = "\\s+")),  :
#  data length [10578] is not a sub-multiple or multiple of the number of rows [625]
# 7: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  :
#  data length [5636] is not a sub-multiple or multiple of the number of rows [625]
# 8: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] + 1)], split = "\\s+")),  :
#  data length [5616] is not a sub-multiple or multiple of the number of rows [625]

# 16 populations----
# Low migration
arp2gen("MSAT_16pops_migLow/MSAT_16pops_migLow_1_1.arp")
arp2gen("MSAT_16pops_migLow/MSAT_16pops_migLow_1_2.arp")
arp2gen("MSAT_16pops_migLow/MSAT_16pops_migLow_1_3.arp")
arp2gen("MSAT_16pops_migLow/MSAT_16pops_migLow_1_4.arp")
arp2gen("MSAT_16pops_migLow/MSAT_16pops_migLow_1_5.arp")
# High migration
arp2gen("MSAT_16pops_migHigh/MSAT_16pops_migHigh_1_1.arp")
arp2gen("MSAT_16pops_migHigh/MSAT_16pops_migHigh_1_2.arp")
arp2gen("MSAT_16pops_migHigh/MSAT_16pops_migHigh_1_3.arp")
arp2gen("MSAT_16pops_migHigh/MSAT_16pops_migHigh_1_4.arp")
arp2gen("MSAT_16pops_migHigh/MSAT_16pops_migHigh_1_5.arp")

# 16pops conversions generate the below warnings
# Warning messages:
# 1: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 2: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 3: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 4: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 5: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 6: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 7: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 8: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 9: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 10: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 11: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 12: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 13: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 14: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 15: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 16: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 17: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 18: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 19: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 20: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 21: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 22: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 23: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 24: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 25: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 26: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 27: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 28: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 29: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [2637] is not a sub-multiple or multiple of the number of rows [157]
# 30: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [2622] is not a sub-multiple or multiple of the number of rows [157]
# 31: In matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]], split = "\\s+")),  ... :
# data length [1430] is not a sub-multiple or multiple of the number of rows [157]
# 32: In matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +  ... :
# data length [1410] is not a sub-multiple or multiple of the number of rows [157]

# Processing genind objects----
# Reading in the genind file generated by diveRsity
genind_MSAT04_migH_1_diveRsity <- read.genepop("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_1.gen")
# Error: "some alleles are not encoded with 2 characters. Check 'ncode' argument"

# Conversion using strataG package----
library(strataG)

arp_DNA04_migH_1 <- arlequinRead("DNA_04pops_migHigh/DNA_04pops_migHigh_1_1.arp")
gtype_DNA04_migH_1 <- arp2gtypes(arp_DNA04_migH_1)
genind_DNA04_migH_1 <- gtypes2genind(gtype_DNA04_migH_1)

str(genind_DNA04_migH_1)

# Reading in the genind file generated by strataG
arp_MSAT04_migH_1 <- arlequinRead("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_1.arp")
# This step errors--generates a null gtypes object
# The code for arp2gtypes doesn't seem to contain an option for the arp$profile.info$data.type=MSAT...
gtype_MSAT04_migH_1 <- arp2gtypes(arp_MSAT04_migH_1)
genind_MSAT04_migH_1 <- gtypes2genind(gtype_MSAT04_migH_1)
# Error: Error in (function (cond)  : 
# error in evaluating the argument 'x' in selecting a method for function 'as.data.frame': 
# Can't find column `id` in `.data`.

# Demo using scen1_migHigh_0 from Kaylee's GH
arp_MSAT05_migH_0 <- arlequinRead("/home/akoontz/Documents/SSRvSNP/Simulations/Demos/Kaylee_REU/scen1_highMig_0.arp")
gtype_MSAT05_migH_0 <- arp2gtypes(arp_MSAT05_migH_0) # Works!
genind_MSAT05_migH_0 <- gtypes2genind(gtype_MSAT05_migH_0)

# Demo using altered MSAT_04pops_mighHigh_1_1.arp, with two fields changed from 0 to 1:
# GenotypicData and GameticPhase
arp_MSAT04_migH_TEST <- arlequinRead("MSAT_04pops_migHigh/MSAT_04pops_migHigh_TEST.arp")
gtype_MSAT04_migH_TEST <- arp2gtypes(arp_MSAT04_migH_TEST) 
genind_MSAT04_migH_TEST <- gtypes2genind(gtype_MSAT04_migH_TEST)

# Demo using altered MSAT_04pops_mighHigh_1_1.arp, with even numbers of samples and loci:
setwd("/RAID1/Simulations/Even_fsc_draftParams/")
arp_MSAT04_migH_Even <- arlequinRead("MSAT_04pops_migHigh/MSAT_04pops_migHigh_1_1.arp")
gtype_MSAT04_migH_Even <- arp2gtypes(arp_MSAT04_migH_Even) 
genind_MSAT04_migH_Even <- gtypes2genind(gtype_MSAT04_migH_Even)
