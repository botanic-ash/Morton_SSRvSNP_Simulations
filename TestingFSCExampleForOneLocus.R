# Parsing fastSimcoal2 outputs

library(adegenet)
library(hierfstat)
library(stringr)
setwd("/home/akoontz/Software/Simulations/fsc27_linux64/exampleFiles/")

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
sim.directory <- "/home/akoontz/Software/Simulations/fsc27_linux64/exampleFiles/"

# 1PopDNA----
# Specify example directory to explore
Example_1PopDNA <- convertAllArp(paste0(sim.directory,"1PopDNA"))

# Number of loci
nLoc(Example_1PopDNA[[1]])
mean(sapply(Example_1PopDNA, function(x) nLoc(x)))
sd(sapply(Example_1PopDNA, function(x) nLoc(x)))

# 1PopDNArec----
# Specify example directory to explore
Example_1PopDNArec <- convertAllArp(paste0(sim.directory,"1PopDNArec"))

# Number of loci
nLoc(Example_1PopDNArec[[1]])

# 1PopDNArec10Mb----
# Specify example directory to explore
# Not run--.arp file is almost 1 GB in size, takes a long time to load
# Example_1PopDNArec10Mb <- convertAllArp(paste0(sim.directory,"1PopDNArec10Mb/"))

# Number of loci
# nLoc(Example_1PopDNArec10Mb[[1]])

# 1PopDNArec100Mb----
# Specify example directory to explore
# Error: Error in (start + 1):(end - 1) : result would be too long a vector
# Example_1PopDNArec100Mb <- convertAllArp(paste0(sim.directory,"1PopDNArec100Mb/"))

# Number of loci
# nLoc(Example_1PopDNArec10Mb[[1]])

# 1PopDNAnoRec10Mb----
# Specify example directory to explore
# Not run--.arp file is almost 1 GB in size, takes a long time to load
# Example_1PopDNAnoRec10Mb <- convertAllArp(paste0(sim.directory,"1PopDNAnoRec10Mb"))

# Number of loci
# nLoc(Example_1PopDNAnoRec10Mb[[1]])

# 1PopDNAnoRec100Mb----
# Specify example directory to explore
# Not run--.arp file is 9.4 GB in size, takes a long time to load
# Example_1PopDNAnoRec100Mb <- convertAllArp(paste0(sim.directory,"1PopDNAnoRec100Mb"))

# Number of loci
# nLoc(Example_1PopDNAnoRec100Mb[[1]])

# 1PopMultiLocus----
# Specify example directory to explore
# Error in arp2gen (strataG): 'SampleData' must have an even number of rows if 'GenotypicData=1'
# Example_1PopMultiLocus <- convertAllArp(paste0(sim.directory,"1PopMultiLocus"))

# Number of loci
# nLoc(Example_1PopMultiLocus[[1]])

# 2PopDNArec----
# Specify example directory to explore
# Error in arp2gen (strataG): 'SampleData' must have an even number of rows if 'GenotypicData=1'
# Example_2PopDNArec <- convertAllArp(paste0(sim.directory,"2PopDNArec"))

# Number of loci
# nLoc(Example_2PopDNArec[[1]])

# 3PopDNAInbreeding----
# Specify example directory to explore
# Error in arp2gen (strataG): 'SampleData' must have an even number of rows if 'GenotypicData=1'
Example_3PopDNAInbreeding <- convertAllArp(paste0(sim.directory,"3PopDNAInbreeding"))

# Number of loci
nLoc(Example_3PopDNAInbreeding[[1]])

# 3PopDNASerial----
# Specify example directory to explore
# Error in arp2gen (strataG): 'SampleData' must have an even number of rows if 'GenotypicData=1'
Example_3PopDNASerial <- convertAllArp(paste0(sim.directory,"3PopDNASerial"))

# Number of loci
nLoc(Example_3PopDNASerial[[1]])

# 3PopDNASFS----
# Specify example directory to explore
# Error in arp2gen (strataG): 'SampleData' must have an even number of rows if 'GenotypicData=1'
Example_3PopDNASFS <- convertAllArp(paste0(sim.directory,"3PopDNASFS"))

# Number of loci
nLoc(Example_3PopDNASFS[[1]])
