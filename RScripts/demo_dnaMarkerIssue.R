# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% DNA MARKER ISSUE DEMONSTRATION %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script demonstrates a bug when simulating DNA marker types using 
# fastSimcoal2 (fsc), and how that bug is solved by running fsc through
# the R package strataG (rather than running fsc standalone). 

# It does this by first generating a single parameters file (via strataG), 
# used in both approaches. This file is reference in the bash script called
# below, and the Arlequin outputs are read into R. The result are loci that
# are concatenated together (30 bp long, rather than 10 bp).

# The same parameters file is then used to run simulations using strataG. The
# results params object is converted to a genind object, and the resulting loci
# are demonstrated (10 bp long as specified).

# In order to test future version of fsc, the code is copied into multiple chunks,
# with each chunk using a different fscVersion variable

library(strataG)
library(adegenet)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs;
# this is a file located on the RAID1 drive, for space reasons, and linked in the home directory)
# The demo_dnaMarkerIssue_fsc.sh BASH script must be present within this folder
sim.wd <- "~/Shared/SSRvSNP_Sim/Code/SimulationOutputs/demo_dnaMarkerIssue/"
setwd(sim.wd)

# ----VARIABLES----
num_reps <- 5
# Single population
demeA <- fscDeme(deme.size = 10, sample.size = 5)
# DNA Genetic parameters
# This should generate 3 blocks (3 loci, for a diploid individual) 
# of DNA sequences of length 10 bp
dna <- fscBlock_dna(sequence.length = 10, mut.rate = 0.0015)
DNAgenetics <- fscSettingsGenetics(dna, dna, dna, num.chrom = 1)

# %%% FASTSIMCOAL VERSION 2709 %%% ----
fscVersion <- "fsc2709"

# %%% GENERATE A PARAMETERS FILE, USED BY BOTH APPROACHES ----
demo_dnaMarkerIssue.params <- fscWrite(demes = fscSettingsDemes(demeA, ploidy = 2), 
                                       genetics = DNAgenetics, 
                                       label = "demo_dnaMarkerIssue", use.wd=TRUE)

# %%% PROBLEM: RUNNING FASTSIMCOAL2 STANDALONE (OUTSIDE STRATAG) ----
# First, run demo_dnaMarkerIssue_fsc.sh script (i.e. bash demo_dnaMarkerIssue_fsc.sh)
bashCall <- "bash demo_dnaMarkerIssue_fsc.sh"
system(command = bashCall)
# Convert the fsc Arlequin output to gtypes object, then genind
fsc_Output <- arlequinRead(paste0(
  sim.wd,"demo_dnaMarkerIssue/demo_dnaMarkerIssue_1_5.arp"))
fsc_Output <- arp2gtypes(fsc_Output)
fsc_Output <- gtypes2genind(fsc_Output)
# fsc output: concatenates 10 bp segments into 30 bp segments
colnames(fsc_Output@tab)

# %%% SOLUTION: RUNNING FASTSIMCOAL2 USING STRATAG ----
# Run fastSimcoal2 from stataG, capturing the output in the params object
demo_dnaMarkerIssue.params <- fscRun(demo_dnaMarkerIssue.params, num.sims = num_reps, 
                                     all.sites = TRUE, exec = fscVersion)

# Convert the strataG params output to gtypes object, then genind
strataG_Output <- fsc2gtypes(demo_dnaMarkerIssue.params, marker = "dna")
strataG_Output <- gtypes2genind(strataG_Output)
# strataG output: properly separates loci into 10 bp segments
colnames(strataG_Output@tab)

# %%% FASTSIMCOAL VERSION 27093 %%% ----
fscVersion <- "fsc27093"

# %%% GENERATE A PARAMETERS FILE, USED BY BOTH APPROACHES ----
demo_dnaMarkerIssue.params <- fscWrite(demes = fscSettingsDemes(demeA, ploidy = 2), 
                                       genetics = DNAgenetics, 
                                       label = "demo_dnaMarkerIssue", use.wd=TRUE)

# %%% PROBLEM: RUNNING FASTSIMCOAL2 STANDALONE (OUTSIDE STRATAG) ----
# First, run demo_dnaMarkerIssue_fsc.sh script (i.e. bash demo_dnaMarkerIssue_fsc.sh)
bashCall <- "bash demo_dnaMarkerIssue_fsc.sh"
system(command = bashCall)
# Convert the fsc Arlequin output to gtypes object, then genind
fsc_Output <- arlequinRead(paste0(
  sim.wd,"demo_dnaMarkerIssue/demo_dnaMarkerIssue_1_5.arp"))
fsc_Output <- arp2gtypes(fsc_Output)
fsc_Output <- gtypes2genind(fsc_Output)
# fsc output: concatenates 10 bp segments into 30 bp segments
colnames(fsc_Output@tab)

# %%% SOLUTION: RUNNING FASTSIMCOAL2 USING STRATAG ----
# Run fastSimcoal2 from stataG, capturing the output in the params object
demo_dnaMarkerIssue.params <- fscRun(demo_dnaMarkerIssue.params, num.sims = num_reps, 
                                     all.sites = TRUE, exec = fscVersion)

# Convert the strataG params output to gtypes object, then genind
strataG_Output <- fsc2gtypes(demo_dnaMarkerIssue.params, marker = "dna")
strataG_Output <- gtypes2genind(strataG_Output)
# strataG output: properly separates loci into 10 bp segments
colnames(strataG_Output@tab)
