# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR SSRvSNP SIMULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used in the Simulations component of the SSRvSNP comparison project.
# Primarily, these functions are used to convert the Arlequin outputs from fastSimcoal2 into
# genind objects that can be processed by the R package adegenet

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

# ---- FUNCTIONS ----
# Function converting Arlequin output to a single genind object (through gtypes format)
strataG_arp2gen <- function(params, repNumber){
  # Extract marker type from params argument
  marker <- params$settings$genetics$fsc.type
  # Read in the Arlequin file, convert it to a gtype object, then to a genind object
  arp <- fscReadArp(params, sim=c(1,repNumber), marker = marker)
  gtype <- df2gtypes(arp, ploidy = 2)
  genind <- gtypes2genind(gtype)
  return(genind)
}

# Function for converting all of the Arlequin files in a directory to genind, generating a list of genind objects
convertAllArp <- function(arp.path, params){
  # Retrieve original working directory, to reset to after conversion
  original.WD <- getwd()
  # Navigate to the folder containing simulation outputs
  setwd(arp.path)
  # Create an empty list object to receive list of genind.
  # The length of this list is the number of replicates, which is specified as a numeric vector
  genind.list <- vector("list",length=length(dir()[str_detect(dir(), pattern = ".arp")]))
  fscReps <- seq(1, length(genind.list))
  # Move up one directory, in order for the fscReadArp command (within strataG_arp2gen) to work
  setwd("..")
  # Convert all Arlequin files to a list of genind objects
  for(i in 1:length(genind.list)){
    genind.obj <- strataG_arp2gen(params, rep=i)
    genind.list[[i]] <- genind.obj
  }
  # Reset to original working directory, and return a list of genind objects
  setwd(original.WD)
  return(genind.list)
}
