# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DEMONSTRATION OF SINGLE LOCUS ISSUE %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script is meant to demonstrate the single locus issue seen in fastSimcoal2
# Arlequin output files using the DNA marker type.
# The issue is first illustrated by reading in simulation outputs NOT processed
# using strataG. Then, the same outputs, processed with strataG, are analyzed/printed

library(adegenet)
library(hierfstat)
library(stringr)
library(strataG)

# %%% FUNCTIONS %%% ----
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

# ANALYZE FSC OUTPUTS USING ARLEQUINREAD (SINGLE DNA LOCUS ISSUE) ----
setwd("/home/akoontz/Software/Simulations/fsc27_linux64/exampleFiles/")
sim.directory <- getwd()
# Specify 1PopDNA example directory to pull Arlequin outputs from
Example_1PopDNA_3 <- arlequinRead(paste0(sim.directory,"/DemoBroken_1PopDNA/DemoBroken_1PopDNA_1_3.arp"))
Example_1PopDNA_3 <- arp2gtypes(Example_1PopDNA_3)
Example_1PopDNA_3 <- gtypes2genind(Example_1PopDNA_3)

# Number of loci is one
nLoc(Example_1PopDNA_3)
# Locus is a long string of nucleotides (should be only 10 bp, but is 30 bp)
colnames(Example_1PopDNA_3@tab)[[1]]

# ANALYZE FSC OUTPUTS USING FSCREADARP (ISSUE SOLVED) ----
# Write deme parameters
Example_1Pop_Deme <- fscSettingsDemes(fscDeme(10000,5))
# Write genetics parameters
DNA_block <- fscBlock_dna(sequence.length = 10,mut.rate = 2e-4, transition.rate = 0.33)
Example_1Pop_DNA <- fscSettingsGenetics(DNA_block, DNA_block, DNA_block, num.chrom = 1)
# Create a fsc parameters file
Example_1PopDNA.params <- fscWrite(demes = Example_1Pop_Deme, 
                                    genetics = Example_1Pop_DNA, 
                                   label = "DemoFix_1PopDNA", use.wd=TRUE)
# Run Example_1PopDNA using strataG, and capture the parameters object output
Example_1PopDNA.params <- fscRun(Example_1PopDNA.params, num.sims = 5, 
                                 all.sites = TRUE, exec = "fsc2709")
# Read in Arlequin outputs and convert to genind
setwd("/home/akoontz/Software/Simulations/fsc27_linux64/exampleFiles/")

Example_1PopDNA_arpPath <- "/home/akoontz/Software/Simulations/fsc27_linux64/exampleFiles/DemoFix_1PopDNA/"

Example_1PopDNA_genind <- convertAllArp(arp.path = Example_1PopDNA_arpPath, params=Example_1PopDNA.params)

# Number of loci is three
nLoc(Example_1PopDNA_genind[[1]])
# Locus is 10 nucleotides long, as specified
colnames(Example_1PopDNA_genind[[1]]@tab)[1:3]
