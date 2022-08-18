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

# Function for randomly assigning a proportion of a genind matrix to a "garden" population (the rest get "wild")
assignGardenSamples <- function(genind.obj, proportion=0.2){
  # Create a vector to store population names (start with all samples being "wild")
  popIDs <- rep("wild",length=(nInd(genind.test)))
  # Get the names of randomly sampled rows of the genind matrix, based on the proportion argument
  gardenSamples <- rownames(genind.obj@tab[sample(nrow(genind.obj@tab), 
                                                  size=nInd(genind.obj)*proportion, replace = FALSE),])
  # Assign randomly sampled rows as "garden"
  popIDs[which(rownames(genind.test@tab) %in% gardenSamples)] <- "garden"
  # Assign pop values and return genind object
  pop(genind.obj) <- popIDs
  return(genind.obj)
}

# Function for reporting ex situ representation rates, using a single genind object
exSituRepresentation <- function(gen.obj){
  # Generate numerical vectors corresponding to garden and wild rows
  garden.Rows <- which(pop(gen.obj)=="garden")
  wild.Rows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector
  wildFreqs <- colSums(gen.obj@tab[wild.Rows,], na.rm = TRUE)/(length(wild.Rows)*2)*100
  # Calculate representation rates
  # Total
  total <- length(which(names(which(wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 0))*100
  # Very common
  veryCommon <- length(which(names(which(wildFreqs > 10)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 10))*100
  # Common
  common <- length(which(names(which(wildFreqs > 5)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 5))*100
  # Low frequency
  lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
  # Rare
  rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[garden.Rows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
  # Build list of rates
  repRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(repRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print representation rates and return
  return(repRates)
}
