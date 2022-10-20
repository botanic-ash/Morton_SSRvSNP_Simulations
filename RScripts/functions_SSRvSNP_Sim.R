# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% FUNCTIONS FOR SSRvSNP SIMULATIONS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script declares the functions used in the Simulations component of the SSRvSNP comparison project.
# Functions are split into three sections:
#   1. Functions used for procesing Arlequin/strataG files. These include functions
#      for converting Arlequin outputs to genind, and for reading in existing strataG
#      params objects and genind objects (so simulations don't need to be run multiple times)
#   2. Commands for measuring ex situ respresentation in the simulated genind files
#   3. Commands for running the resampling analyses, in the simulated genind files


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
  # [length(dir(pattern))] means if multiple params objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.params"), readRDS(
    dir(pattern = "^params.MSAT_01pop_migLow")[length(dir(pattern = "^params.MSAT_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.params"), readRDS(
    dir(pattern = "^params.MSAT_01pop_migHigh")[length(dir(pattern = "^params.MSAT_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.params"), readRDS(
    dir(pattern = "^params.MSAT_04pop_migLow")[length(dir(pattern = "^params.MSAT_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.params"), readRDS(
    dir(pattern = "^params.MSAT_04pop_migHigh")[length(dir(pattern = "^params.MSAT_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.params"), readRDS(
    dir(pattern = "^params.MSAT_16pop_migLow")[length(dir(pattern = "^params.MSAT_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.params"), readRDS(
    dir(pattern = "^params.MSAT_16pop_migHigh")[length(dir(pattern = "^params.MSAT_16pop_migHigh"))]), pos = 1)
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
  # [length(dir(pattern))] means if multiple params objects are present in the directory, read the most recent one
  assign(paste0(prefix, "_01pop_migLow.params"), readRDS(
    dir(pattern = "^params.DNA_01pop_migLow")[length(dir(pattern = "^params.DNA_01pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_01pop_migHigh.params"), readRDS(
    dir(pattern = "^params.DNA_01pop_migHigh")[length(dir(pattern = "^params.DNA_01pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migLow.params"), readRDS(
    dir(pattern = "^params.DNA_04pop_migLow")[length(dir(pattern = "^params.DNA_04pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_04pop_migHigh.params"), readRDS(
    dir(pattern = "^params.DNA_04pop_migHigh")[length(dir(pattern = "^params.DNA_04pop_migHigh"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migLow.params"), readRDS(
    dir(pattern = "^params.DNA_16pop_migLow")[length(dir(pattern = "^params.DNA_16pop_migLow"))]), pos = 1)
  assign(paste0(prefix, "_16pop_migHigh.params"), readRDS(
    dir(pattern = "^params.DNA_16pop_migHigh")[length(dir(pattern = "^params.DNA_16pop_migHigh"))]), pos = 1)
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

# EX SITU REPRESENTATION ----
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

# Function for generating a vector of wild allele frequencies from a genind object
getWildFreqs <- function(gen.obj, wholeValues=TRUE){
  # Build a vector of rows corresponding to wild individuals (those that do not have a population of "garden")
  wildRows <- which(pop(gen.obj)!="garden")
  # Build the wild allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  # wholeValues argument determines whether to return whole percentages or fractions
  if(wholeValues==TRUE){
    wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)*100
  } else{
    # (Fractions are more useful when trying to generate histograms)
    wildFreqs <- colSums(gen.obj@tab[wildRows,], na.rm = TRUE)/(length(wildRows)*2)
  }
  return(wildFreqs)
}

# Function for generating a vector of total allele frequencies from a genind object
getTotalFreqs <- function(gen.obj){
  # Build allele frequency vector: colSums of alleles (removing NAs), divided by number of haplotypes (Ne*2)
  totalFreqs <- colSums(gen.obj@tab, na.rm = TRUE)/(nInd(gen.obj)*2)*100
  return(totalFreqs)
}

# Exploratory function for reporting the proprtion of alleles of each category, from a (wild) frequency vector
getWildAlleleFreqProportions <- function(gen.obj){
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

# Exploratory function for reporting the proprtion of alleles of each category, from a frequency vector (of ALL alleles--garden AND wild)
getTotalAlleleFreqProportions <- function(gen.obj){
  # Build the total allele frequency vector, using the getTotalFreqs function
  totalFreqs <- getTotalFreqs(gen.obj)
  # Very common
  veryCommonAlleles <- totalFreqs[which(totalFreqs > 10)]
  veryCommon_prop <- (length(veryCommonAlleles)/length(totalFreqs))*100
  # Low frequency
  lowFrequencyAlleles <- totalFreqs[which(totalFreqs < 10 & totalFreqs > 1)]
  lowFrequency_prop <- (length(lowFrequencyAlleles)/length(totalFreqs))*100
  # Rare
  rareAlleles <- totalFreqs[which(totalFreqs < 1)]
  rare_prop <- (length(rareAlleles)/length(totalFreqs))*100
  # Build list of proportions, and return
  freqProportions <- c(veryCommon_prop, lowFrequency_prop, rare_prop)
  names(freqProportions) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  return(freqProportions)
}

# Function for summarizing allele frequency proportions across replicates
summarize_alleleFreqProportions <- function(freqProportions){
  # Calculate the mean allele frequency proportions across replicates using apply
  # (Rows are allele frequency categories, columns are replicates. So margin value is 1)
  means <- apply(freqProportions, 1, mean, na.rm=TRUE)
  # Calculate standard deviations
  stdevs <- apply(freqProportions, 1, sd, na.rm=TRUE)
  # Combine statistics into a matrix, and return
  freqPropStats <- cbind(means, stdevs)
  return(freqPropStats)
}

# Function for reporting representation rates, using a sample matrix and a vector of allele frequencies
# This function assumes that the freqVector represents the absolute allele frequencies
# for the entire population of interest. Allele names between the frequency vector and the sample matrix
# must correspond in order for values to be comparable.
getAlleleCategories <- function(freqVector, sampleMat){
  # Determine how many total alleles in the sample (i.e. greater than 0) are found in the frequency vector 
  total <- length(which(names(which(freqVector > 0)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector > 0))*100
  # Very common alleles (greater than 10%)
  v_com <- length(which(names(which(freqVector > 10)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector > 10))*100
  # Common alleles (greater than 5%)
  com <- length(which(names(which(freqVector > 5)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector > 5))*100
  # Low frequency alleles (between 1% and 10%)
  low_freq <- length(which(names(which(freqVector < 10 & freqVector > 1)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector < 10 & freqVector > 1))*100
  # Rare alleles (less than 1%)
  rare <- length(which(names(which(freqVector < 1 & freqVector > 0)) %in% names(which(colSums(sampleMat, na.rm = TRUE) > 0))))/length(which(freqVector < 1 & freqVector > 0))*100
  # Concatenate values to a vector, name that vector, and return
  exSituRepRates <- c(total, v_com, com, low_freq, rare) 
  names(exSituRepRates) <- c("Total","Very common (>10%)","Common (>5%)",
                             "Low frequency (1% -- 10%","Rare (<1%)")
  return(exSituRepRates)
}

# # Function for reporting ex situ representation rates, using a single genind object
# exSituRepresentation_OLD <- function(gen.obj){
#   # Generate numerical vectors corresponding to garden and wild rows
#   gardenRows <- which(pop(gen.obj)=="garden")
#   wildRows <- which(pop(gen.obj)!="garden")
#   # Build the wild allele frequency vector, using the getWildFreqs function
#   wildFreqs <- getWildFreqs(gen.obj)
#   # Calculate representation rates
#   # Total
#   total <- length(which(names(which(wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 0))*100
#   # Very common
#   veryCommon <- length(which(names(which(wildFreqs > 10)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 10))*100
#   # Common
#   common <- length(which(names(which(wildFreqs > 5)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs > 5))*100
#   # Low frequency
#   lowFrequency <- length(which(names(which(wildFreqs < 10 & wildFreqs > 1)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 10 & wildFreqs > 1))*100
#   # Rare
#   rare <- length(which(names(which(wildFreqs < 1 & wildFreqs > 0)) %in% names(which(colSums(gen.obj@tab[gardenRows,], na.rm = TRUE) > 0))))/length(which(wildFreqs < 1 & wildFreqs > 0))*100
#   # Build list of rates
#   repRates <- c(total,veryCommon,common,lowFrequency,rare)
#   names(repRates) <- c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
#   # Return vector of rates
#   return(repRates)
# }

# Wrapper function for reporting ex situ representation rates from genind object
exSituRepresentation <- function(gen.obj){
  # Build the wild allele frequency vector, using the getWildFreqs function
  wildFreqs <- getWildFreqs(gen.obj)
  # Build the matrix of garden samples
  gardenMat <- gen.obj@tab[which(pop(gen.obj)=="garden"),]
  # Calculate representation rates, using getAlleleCategories function
  repRates <- getAlleleCategories(freqVector = wildFreqs, sampleMat = gardenMat)
  # Return vector of rates
  return(repRates)
}

# Function for summarizing ex situ representation rates across replicates
summarize_exSituRepresentation <- function(repRates){
  # Calculate the mean ex situ representation rate across replicates using apply
  # (Rows are rate categories, columns are replicates. So margin value is 1)
  means <- apply(repRates, 1, mean, na.rm=TRUE)
  # Calculate standard deviations
  stdevs <- apply(repRates, 1, sd, na.rm=TRUE)
  # Combine statistics into a matrix, and return
  repStats <- cbind(means, stdevs)
  return(repStats)
}

# Wrapper function, which generates both allele frequency proportions and ex situ representation rates,
# for a list of genind objects
summarize_simulations <- function(genind.list, gardenRate=0.05){
  # Build array to capture allele frequency proportions
  alleleFreqSummaries <- array(dim = c(3, 2, length(genind.list)))
  rownames(alleleFreqSummaries) <- c("Very common (>10%)","Low frequency (1% -- 10%)","Rare (<1%)")
  # Build array to capture ex situ representation rates
  repRateSummaries <- array(dim = c(5, 2, length(genind.list)))
  rownames(repRateSummaries) <- 
    c("Total","Very common (>10%)","Common (>5%)","Low frequency (1% -- 10%)","Rare (<1%)")
  colnames(alleleFreqSummaries) <- colnames(repRateSummaries) <-c("mean", "sd")
  
  # Loop through list of genind objects, calculating metrics for each item
  for (i in 1:length(genind.list)){
    # Calculate and summarize allele frequency scenarios. Each array slot is a different scenario
    alleleFrequencies <- sapply(genind.list[[i]], getWildAlleleFreqProportions)
    alleleFreqSummaries[,,i] <- summarize_alleleFreqProportions(alleleFrequencies)
    # Assign individuals to garden population
    genind.list[[i]] <- lapply(genind.list[[i]], assignGardenSamples, proportion=gardenRate)
    # Calculate and summarize ex situ representation rates. Each array slot is a different scenario
    representationRates <- sapply(genind.list[[i]], exSituRepresentation)
    repRateSummaries[,,i] <- summarize_exSituRepresentation(representationRates)
  }
  # Round results to 2 digits
  alleleFreqSummaries <- round(alleleFreqSummaries, 2)  
  repRateSummaries <- round(repRateSummaries, 2)
  # Generate a list of the two arrays, and return
  return(list("alleleFrequencyProportions"=alleleFreqSummaries, "representationRates"=repRateSummaries))
}

# Exploratory function, which creates a histogram of allele frequencies from a genind object
makeAlleleFreqHist <- function(gen.obj, title="Allele frequency histogram"){
  # Make a vector of allele frequency values, from the genind object
  wildAlleleFreqs <- getWildFreqs(gen.obj, wholeValues = FALSE)
  # Set the break values according to the maximum and minimum allele
  max.freq <- max
  # Specify the break values to use for the histogram conditionally, 
  # based on min/max frequencies
  if(max(wildAlleleFreqs == 1)){
    if(min(wildAlleleFreqs == 0)){
      b.vector <- seq(0, 1.0, 0.01)
    } else{
      b.vector <- seq(0.01, 1.0, 0.01)
    }
  } else{
    if(min(wildAlleleFreqs == 0)){
      b.vector <- seq(0, 1.0, 0.01)
    } else{
      b.vector <- seq(0.01, 1.0, 0.01)
    }
  }
  # Generate histogram
  hist(wildAlleleFreqs, breaks=b.vector, main=title, freq=FALSE)
}

# RESAMPLING FUNCTIONS ----
# Ex situ sample function, which finds the level of ex situ representation of a sample of individuals
# (using the getAlleleCategories function above)
exSitu_Sample <- function(wildMat, numSamples){
  # Build the wild allele frequency vector, using the getWildFreqs function
  freqVector <- getWildFreqs(gen.obj)
  # From a matrix of individuals, select a set of random individuals (rows)
  samp <- wildMat[sample(nrow(wildMat), size=numSamples, replace = FALSE),]
  # Calculate how many alleles (of each category) that sample captures, and return
  repRates <- getAlleleCategories(freqVector, samp)
  return(repRates)
}

# Wrapper for the exSitu_Sample function, iterating that function over the entire sample matrix
exSitu_Resample <- function(gen.obj){
  # Create a matrix of wild individuals (those with population "wild") from genind object
  wildMat <- gen.obj@tab[which(pop(gen.obj) == "wild"),]
  # Apply the exSituSample function to all rows of the sample matrix
  # (except row 1, because we need at least 2 individuals to sample)
  # Resulting matrix needs to be transposed in order to keep columns as different allele categories
  representationMatrix <- t(sapply(2:nrow(wildMat), function(x) exSitu_Sample(wildMat, x)))
  # Name columns according to categories of allelic representation, and return matrix
  colnames(representationMatrix) <- c("Total","Very common","Common","Low frequency","Rare")
  return(representationMatrix)
}

# !!! WORKING !!!
# Wrapper for exSitu_Resample, which iterates the function for a list of genind objects
Resample_Test <- function(gen.list){
  sapply( )
  lapply(gen.list, exSitu_Resample)
}
