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

# TO DO: INCORPORATE PREFIX ARGUMENTS TO ALLOW FOR DIFFERENT PARAMS/GENIND OBJECTS TO BE READ IN

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

# Function for reading in MSAT strataG params files, in specified directory. Prefix specifies how to name the params variables
readParams_MSAT <- function(params.wd, prefix="MSAT"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing strataG params objects
  setwd(params.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple genind objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.MSAT_01pop_migLow")[length(dir(pattern = "^genind.MSAT_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.MSAT_01pop_migHigh")[length(dir(pattern = "^genind.MSAT_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.MSAT_04pop_migLow")[length(dir(pattern = "^genind.MSAT_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.MSAT_04pop_migHigh")[length(dir(pattern = "^genind.MSAT_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.MSAT_16pop_migLow")[length(dir(pattern = "^genind.MSAT_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.MSAT_16pop_migHigh")[length(dir(pattern = "^genind.MSAT_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# Function for reading in MSAT genind files, in specified directory. Prefix specifies how to name the genind variables
readGeninds_MSAT <- function(geninds.wd, prefix="MSAT"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing genind objects
  setwd(geninds.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple genind objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.MSAT_01pop_migLow")[length(dir(pattern = "^genind.MSAT_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.MSAT_01pop_migHigh")[length(dir(pattern = "^genind.MSAT_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.MSAT_04pop_migLow")[length(dir(pattern = "^genind.MSAT_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.MSAT_04pop_migHigh")[length(dir(pattern = "^genind.MSAT_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.MSAT_16pop_migLow")[length(dir(pattern = "^genind.MSAT_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.MSAT_16pop_migHigh")[length(dir(pattern = "^genind.MSAT_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# Function for reading in DNA strataG params files, in specified directory. Prefix specifies how to name the params variables
readParams_DNA <- function(params.wd, prefix="DNA"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing strataG params objects
  setwd(params.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple genind objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.DNA_01pop_migLow")[length(dir(pattern = "^genind.DNA_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.DNA_01pop_migHigh")[length(dir(pattern = "^genind.DNA_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.DNA_04pop_migLow")[length(dir(pattern = "^genind.DNA_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.DNA_04pop_migHigh")[length(dir(pattern = "^genind.DNA_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.DNA_16pop_migLow")[length(dir(pattern = "^genind.DNA_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.DNA_16pop_migHigh")[length(dir(pattern = "^genind.DNA_16pop_migHigh"))]), pos = 1)
  # Reset to original working directory
  setwd(original.wd)
}

# Function for reading in DNA genind files, in specified directory. Prefix specifies how to name the genind variables
readGeninds_DNA <- function(geninds.wd, prefix="DNA"){
  # Retrieve original working directory, to reset to after conversion
  original.wd <- getwd()
  # Navigate to the folder containing genind objects
  setwd(geninds.wd)
  # Use assign to assign the value (2nd argument) to the string (1st argument) specified. Prefix argument
  # allows variable names to be specifiable, but all variables will end with the same string
  # The pos=1 argument allows the variable to be passed to the global environment (rather than being kept locally)
  # The ^ character in dir allows any file with the given suffix to be read in, allowing for file name flexibility
  # [length(dir(pattern))] means if multiple genind objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.DNA_01pop_migLow")[length(dir(pattern = "^genind.DNA_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.DNA_01pop_migHigh")[length(dir(pattern = "^genind.DNA_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.DNA_04pop_migLow")[length(dir(pattern = "^genind.DNA_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.DNA_04pop_migHigh")[length(dir(pattern = "^genind.DNA_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.genind"), readRDS(
    dir(pattern = "^genind.DNA_16pop_migLow")[length(dir(pattern = "^genind.DNA_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.genind"), readRDS(
    dir(pattern = "^genind.DNA_16pop_migHigh")[length(dir(pattern = "^genind.DNA_16pop_migHigh"))]), pos = 1)
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

# Function for summarizing allele frequency proportions across replicates
summarize_alleleFreqProportions <- function(freqProportions){
  # Calculate the mean allele frequency proportions across replicates using apply
  # (Rows are allele frequency categories, columns are replicates. So margin value is 1)
  means <- apply(freqProportions, 1, mean)
  # Calculate standard deviations
  stdevs <- apply(freqProportions, 1, sd)
  # Combine statistics into a matrix, and return
  freqPropStats <- cbind(means, stdevs)
  return(freqPropStats)
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
  # Return vector of rates
  return(repRates)
}

# Function for summarizing ex situ representation rates across replicates
summarize_exSituRepresentation <- function(repRates){
  # Calculate the mean ex situ representation rate across replicates using apply
  # (Rows are rate categories, columns are replicates. So margin value is 1)
  means <- apply(repRates, 1, mean)
  # Calculate standard deviations
  stdevs <- apply(repRates, 1, sd)
  # Combine statistics into a matrix, and return
  repStats <- cbind(means, stdevs)
  return(repStats)
}
