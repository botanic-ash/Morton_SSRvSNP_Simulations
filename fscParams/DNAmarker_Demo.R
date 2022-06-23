# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DNA MARKER DEMONSTRATION USING STRATAG %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script uses strataG to create the fastsimcoal2 (fsc) parameter files with DNA markers
# It then demonstrates the issue of a single locus, using these DNA markers (in ANALYZING OUTPUTS):
# Multiple loci appear in the DNAmarker_Demo object (correctly), but trying to convert this to a genind object errors
# When Arlequin outputs are converted to genind objects, only 1 locus is seen in the DNAmarker_Demo_genind object

library(strataG)
library(adegenet)
# Set below line to whatever working directory is necessary
setwd("~/Documents/SSRvSNP/Simulations/Code/fscParams/")

# ----VARIABLES----
num_reps <- 5
fscVersion <- "fsc2702"
# Single population
demeA <- fscDeme(deme.size = 10, sample.size = 10)
# DNA Genetic parameters
# This should generate 6 blocks (3 loci, for a diploid individual) of DNA sequences of length 15 bp
dna <- fscBlock_dna(sequence.length = 15, mut.rate = 1e-3)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, num.chrom = 1)

# ----RUNNING FASTSIMCOAL2----
# Generate a parameter file
DNAmarker_Demo.params <- fscWrite(demes = fscSettingsDemes(demeA, ploidy = 2), genetics = DNAgenetics, 
                                  label = "DNAmarker_Demo_Diploid", use.wd=TRUE)
# Run the parameters file
DNAmarker_Demo.params <- fscRun(DNAmarker_Demo.params, num.sims = num_reps, all.sites = TRUE, exec = fscVersion)

# ----ANALYZING OUTPUTS----
# To convert to gtype, then genind
DNAmarker_Demo_gtype <- fsc2gtypes(DNAmarker_Demo.params,marker = "dna")
DNAmarker_Demo_genind <- gtypes2genind(DNAmarker_Demo_gtype)
DNAmarker_Demo_genind@tab
