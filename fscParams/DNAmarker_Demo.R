# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DNA MARKER DEMONSTRATION USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses strataG to create the fastsimcoal2 (fsc) parameter files with DNA markers
# It then demonstrates the issue of single loci, using these DNA markers
# Multiple loci appear in the DNAmarker_Demo object, but only 1 in the DNAmarker_Demo_genind object

library(strataG)
library(adegenet)
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/")

# ----VARIABLES----
num_reps <- 5
fscVersion <- "fsc2702"
# Single population
demeA <- fscDeme(deme.size = 60, sample.size = 60)
# DNA Genetic parameters
# This should generate 8 blocks (4 loci, for a diploid individual) of DNA sequences of length 15 bp
dna <- fscBlock_dna(sequence.length = 15, mut.rate = 1e-4)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, num.chrom = 1)

# Generate a parameter file
DNAmarker_Demo.params <- fscWrite(demes = fscSettingsDemes(demeA), genetics = DNAgenetics, 
                                  label = "DNAmarker_Demo", use.wd=TRUE)
# Run the parameters file
DNAmarker_Demo.params <- fscRun(DNAmarker_Demo.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)

# Read the Arlequin outputs using strataG
DNAmarker_Demo <- fscReadArp(DNAmarker_Demo.params, sim = c(1,1), marker = "dna")
dim(DNAmarker_Demo)

# Can't convert fsc params to gtype: "Error: the number of genes in 'sequences' is not equal to the number of loci"
DNAmarker_Demo_gtype <- fsc2gtypes(DNAmarker_Demo.params,marker = "dna")

# When converting Arlequin file to genind, same issue: single locus
DNAmarker_Demo_Arlequin <- arlequinRead("DNAmarker_Demo/DNAmarker_Demo_1_3.arp")
DNAmarker_Demo_gtype <- arp2gtypes(DNAmarker_Demo_Arlequin)
DNAmarker_Demo_genind <- gtypes2genind(DNAmarker_Demo_gtype)
nLoc(DNAmarker_Demo_genind)
