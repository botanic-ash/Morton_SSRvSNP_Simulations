# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR SSRvSNP SIMULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used in the Simulations component of the SSRvSNP comparison project.
# Functions are split into two sections. The first consists of commands used to process fastSimcoal2
# and strataG output files. The second consists of commands for measuring ex situ representation from
# the simulated genind files.

library(strataG)
library(adegenet)
library(stringr)
library(hierfstat)

# ---- FUNCTIONS ----
# PROCESSING ARLEQUIN/STRATAG FILES ----
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
  original.wd <- getwd()
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
  setwd(original.wd)
  return(genind.list)
}

# Hard-coded function for reading in strataG MSAT params files, in the specified directory
readParams_MSAT <- function(params.wd){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing strataG params objects
  setwd(params.wd)
  # The ^ character allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple params objects are present in the directory, read the most recent one
  # <<- means pass the variable to the global environment (otherwise, variable will be kept local, and unseen)
  MSAT_01pop_migLow.params <<- readRDS(
    dir(pattern = "^params.MSAT_01pop_migLow")[length(dir(pattern = "^params.MSAT_01pop_migLow"))])
  MSAT_01pop_migHigh.params <<- readRDS(
    dir(pattern = "^params.MSAT_01pop_migHigh")[length(dir(pattern = "^params.MSAT_01pop_migHigh"))])
  MSAT_04pop_migLow.params <<- readRDS(
    dir(pattern = "^params.MSAT_04pop_migLow")[length(dir(pattern = "^^params.MSAT_04pop_migLow"))])
  MSAT_04pop_migHigh.params <<- readRDS(
    dir(pattern = "^params.MSAT_04pop_migHigh")[length(dir(pattern = "^params.MSAT_04pop_migHigh"))])
  MSAT_16pop_migLow.params <<- readRDS(
    dir(pattern = "^params.MSAT_16pop_migLow")[length(dir(pattern = "^params.MSAT_16pop_migLow"))])
  MSAT_16pop_migHigh.params <<- readRDS(
    dir(pattern = "^params.MSAT_16pop_migHigh")[length(dir(pattern = "^params.MSAT_16pop_migHigh"))])
  # Reset to original working directory
  setwd(original.wd)
}

# Hard-coded function for reading in strataG DNA params files, in the specified directory
readParams_DNA <- function(params.wd){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing strataG params objects
  setwd(params.wd)
  # The ^ character allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple params objects are present in the directory, read the most recent one
  # <<- means pass the variable to the global environment (otherwise, variable will be kept local, and unseen)
  DNA_01pop_migLow.params <<- readRDS(
    dir(pattern = "^params.DNA_01pop_migLow")[length(dir(pattern = "^params.DNA_01pop_migLow"))])
  DNA_01pop_migHigh.params <<- readRDS(
    dir(pattern = "^params.DNA_01pop_migHigh")[length(dir(pattern = "^params.DNA_01pop_migHigh"))])
  DNA_04pop_migLow.params <<- readRDS(
    dir(pattern = "^params.DNA_04pop_migLow")[length(dir(pattern = "^^params.DNA_04pop_migLow"))])
  DNA_04pop_migHigh.params <<- readRDS(
    dir(pattern = "^params.DNA_04pop_migHigh")[length(dir(pattern = "^params.DNA_04pop_migHigh"))])
  DNA_16pop_migLow.params <<- readRDS(
    dir(pattern = "^params.DNA_16pop_migLow")[length(dir(pattern = "^params.DNA_16pop_migLow"))])
  DNA_16pop_migHigh.params <<- readRDS(
    dir(pattern = "^params.DNA_16pop_migHigh")[length(dir(pattern = "^params.DNA_16pop_migHigh"))])
  # Reset to original working directory
  setwd(original.wd)
}

# EX SITU REPRESENTATION/RESAMPLING ----
# Function for randomly assigning a proportion of a genind matrix to a "garden" population (the rest get "wild")
assignGardenSamples <- function(genind.obj, proportion=0.2){
  # Create a vector to store population names (start with all samples being "wild")
  popIDs <- rep("wild",length=(nInd(genind.obj)))
  # Get the names of randomly sampled rows of the genind matrix, based on the proportion argument
  gardenSamples <- rownames(genind.obj@tab[sample(nrow(genind.obj@tab), 
                                                  size=nInd(genind.obj)*proportion, replace = FALSE),])
  # Assign randomly sampled rows as "garden"
  popIDs[which(rownames(genind.obj@tab) %in% gardenSamples)] <- "garden"
  # Assign pop values and return genind object
  pop(genind.obj) <- popIDs
  return(genind.obj)
}

# Function for generating a vector of allele frequencies from a genind object
getWildFreqs <- function(gen.obj){
  # Build a vector of rows corresponding to wild individuals (those that do not have a population of "garden")
  wildRows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)*100
  return(wildFreqs)
}

# Exploratory function for reporting the proprtion of alleles of each category, from a frequency vector
getAlleleFreqProportions <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  wildFreqs <- getWildFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- wildFreqs[which(wildFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(wildFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- wildFreqs[which(wildFreqs < 10 & wildFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(wildFreqs))*100
  # Rare
  rareAlleles <- wildFreqs[which(wildFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(wildFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# Function for reporting ex situ representation rates, using a single genind object
exSituRepresentation <- function(gen.obj){
  # Generate numerical vectors corresponding to garden and wild rows
  gardenRows <- which(pop(gen.obj)=="garden")
  wildRows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector, using the getWildFreqs function
  wildFreqs <- getWildFreqs(gen.obj)
  # Calculate representation rates
  # Total
  total <- length(which(names(which(wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 0))*100
  # Very common
  veryCommon <- length(which(names(which(wildFreqs > 10)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 10))*100
  # Common
  common <- length(which(names(which(wildFreqs > 5)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 5))*100
  # Low frequency
  lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
  # Rare
  rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
  # Build list of rates
  repRates <- c(total,veryCommon,common,lowFrequency,rare)
  names(repRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Print representation rates and return
  return(repRates)
}
