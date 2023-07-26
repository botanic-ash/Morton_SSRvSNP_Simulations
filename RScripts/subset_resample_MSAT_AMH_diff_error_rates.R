# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE MSAT SIMULATED DATASETS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Some major edits had to be made to this script because strataG would not download correctly onto Ash's computer. Aditionally, to make the data more human readable (tables/dfs instead of genenid files) and easier to manipulate, some code was copied and edited from Kaylee Rosenberger's Pollen_dispersal_sims github repository.

# This script reads in genind files generated from previously run fastSimcoal simulations, then subsets each file to specify a group of "ex situ" (garden) individuals. The ex situ representation (i.e. how well do these garden individuals represent total allelic diversity) is calculated and visualized. Then, the remaining individuals ("wild") are resampled iteratively, and the allelic diversity of sample subsets (in comparison to the whole of wild allelic diversity) is calculated, then plotted

# In order to function iteratively over large objects (i.e. lists of genind objects), the steps in this script use many apply family functions

#library(strataG) #for SNPS, doesn't work on Ash's computer
# Austin's packages
library(adegenet)
library(stringr)
library(parallel)
library(RColorBrewer)
library(scales)

# Kaylee's packages
library(poppr)
library(tidyr)
library(dplyr)
library(data.table)

# Ash's packages
library(ggplot2)
library(forcats)
library(ggnewscale)
library(MuMIn)
library(betareg)
library(parallelly)
library(future.apply)

plan(multisession, workers = detectCores() - 1) # A part of the futures package which will enable my apply functions to be paralleleized regardless of the OS I am using 

source("/Users/Ashley/Desktop/Fall\ 2019/stats/diagPlots.R") # This is from my stats prof at TAMUCC, it provides the SW p value to tell if things are significantly non-normal

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
sim.wd <- "/Users/Ashley/Desktop/Hoban/Hoban_rotation_2023/Morton_SSRvSNP_Simulations/"
arq.wd <- "/Users/Ashley/Desktop/Hoban/Hoban_rotation_2023/Morton_SSRvSNP_Simulations/SimulationOutputs/MSAT_marker/MSAT_ArlequinFiles"
setwd(sim.wd)

# Load custom functions
source(paste0(sim.wd, "RScripts/necessary_functions.R"))


####Loading the already simulated population data####
num_loci = 20 #number of loci simulated, needed to make a data frame to save the data, starts with 20 --> gets cut down to 10 later

loci <- seq(1, num_loci)

# Make a vector of all of the loci names
loci_names = c()
for(i in 1:num_loci){
  loci_names = c(loci_names, paste("locus", i, "a", sep=""))
  loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}

# Importing and converting arlequin files to genepop files
import_arp2gen_files(arq.wd,".arp$")

# Importing and converting genepop files to genalex
import_gen2genalex_files(arq.wd, ".gen$")

#list of genalex files for all simulation replicates--genalex files end in .csv
genalex_list = list.files(arq.wd, ".csv$")

#read in and then cut off first 2 rows in each data frame -- this is the population data, which is not required for our purposes
i = genalex_list[[1]]

# NEED TO START AGAIN FROM UP HERE AND MAKE THE BELOW CODE COMPATIBLE TO RUN WITH BOTH DATASET TYPES --> WILL PROB MAKE ALL FIGURES 2X AND THEN ADD SCENARIO TYPE TO MY ANOVAS
final_output <- lapply(genalex_list, function(i){
  genetic_data <- read.csv(paste(arq.wd, "/", i, sep=""), header=FALSE)[-c(1,2),]
  #giving the data frame columns new names
  names(genetic_data) = c("Ind", "Pop", loci_names)
  
#this genetic_data object now contains all of the data of 20 MSAT loci from each of individuals (in a single simulation) as simulated from SIM2COAL
  
  genetic_data = as.data.frame(genetic_data[-1,]) #removing the first row and converting to a dataframe
  }
)


# Make current genetic data from a single genind file
genetic_data <- read.csv(paste(arq.wd, "/", genalex_list[1], sep=""), header=FALSE)[-c(1,2),]

names(genetic_data) = c("Ind", "Pop", loci_names) #Give the data frame columns new name
genetic_data = as.data.frame(genetic_data[-1,]) # Remove the first row and converting to a dataframe


#### All things below are done to a single entry in the genetic_data list ####

#plan is to make the simulating error code output a list of dataframes so that everything else can be wrapped into an lapply call and then at the end I can rbind my data together, will just need to make sure I have a column which specifies which replicate it is and can then just group_by replicate and get medians   

### Adding error to the data

# Set error rates, NEED TO GIVE DIFF NAMES
syst_err_rate = c(0, .01, .02, .04, .08, .16) # Based off of rates of errored loci found in Austin's real SNP datasets where source of error is unknown
stoch_err_rate = c(0, .01, .02, .04, .08, .16)

# Make a list of all unique error rate combinations
rate_list <- expand.grid(syst_err_rate, stoch_err_rate) %>%
  rename(syst_err_rate = Var1, 
         stoch_err_rate = Var2) %>%
  filter(rowSums(across(everything())) > 0)

rate_combo_count <- seq(1, nrow(rate_list), 1) # Make a list corresponding to the number of unique error rate combinations

# Set number of reps per unique error rate combination
rep_num = 5 
reps = as.list(seq(1, rep_num, by = 1)) # Make a list corresponding to the number of reps per error rate combination

# Apply the simulating_error function to the input genetic data across the list of unique replicates and then across the list of unique rate value combinations
# Output is a list of lists of lists (3 nested lists):
#   1st list contains all of the below at each unique error rate combo
#   2nd list contains all of the below at each unique error replicate (with the given error rate value specified by above list)
#   3rd list contains 2 data frames: the dataset of every individuals' alleles and a dataset that indicates where errors have resulted in changed alleles for every individual
genetic_data_error <- future_lapply(rate_combo_count, future.seed=TRUE, function(x)
  lapply(reps, function(i)
       simulating_error(data_for_editing = genetic_data, syst_err_rate = rate_list$syst_err_rate[x], stoch_err_rate = rate_list$stoch_err_rate[x], msat_length = 1, replicate_number = i)))


#### Summarizing error added to the data ####
# Make a data frame from the tracked error dataset from each rep from each unique rate value combo and clean the dataframe  
genetic_data_track_error <- as.data.frame(do.call(rbind, lapply(genetic_data_error, function(x)
  do.call(rbind, lapply(x, function(i) 
  colSums(i[[2]])/(nrow(i[[2]])*2)))))) %>% # Calc total error rate at each locus
  # Need to multiply the last few columns by 2 bc divided by 2 x the number of rows to calculate the true loci error
  mutate(rep = as.factor(rep *2), 
         syst_err_rate = as.factor(syst_err_rate *2), 
         stoch_err_rate = as.factor(stoch_err_rate *2)) %>%
  pivot_longer(cols = starts_with("locus"), names_to = "locus", values_to = "error_rate") # Make it so each locus is a unique row instead of a unique column

# Summarize the error data at each unique error rate combo
# pull the first ind for each rep when 95% diversity is captured
obs_error_rate <- genetic_data_track_error %>%
  group_by(locus, syst_err_rate, stoch_err_rate) %>%
  group_by(syst_err_rate, stoch_err_rate) %>%
  summarize(medians = median(error_rate))

#Plot the summarized error data at each locus: faceted by unique error rate combo, replicates condensed into boxplots, all 20 loci on the x axis
genetic_data_track_error %>%
  ggplot() +
  geom_boxplot(aes(x = locus, y = error_rate), outlier.shape = NA) + 
  geom_hline(data = obs_error_rate, aes(yintercept = medians, group = interaction(syst_err_rate, stoch_err_rate)), alpha = .3, color = "red") +
  geom_text(data = obs_error_rate, aes(y = medians + .05, label = round(medians, 2), x = 10, angle = 0, vjust = -0.2, group = interaction(syst_err_rate, stoch_err_rate)), size = 2.5, color = "red") +
  facet_wrap(~ syst_err_rate * stoch_err_rate) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Plot the summarized error data at each locus: one boxplot represents all replicates at all loci at a single unique error rate combo
genetic_data_track_error %>%
  ggplot() +
  geom_boxplot(aes(y = error_rate, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  ylab("observed error rate") +
  xlab("stochastic error rate") +
  theme_classic()


#### Comparing whole error dataset to non-errored dataset ####

# Get the frequencies of each allele from each locus from the "real" dataset
allele_freqs_real <- get_alleles_and_freqs(input_data = genetic_data)

# Get the number of alleles for each locus into a data frame
allele_counts_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_real, function(x) lengths(x)[1])))

# Make the allele_counts_df into a df that can be compared to each rep of each unique rate combo (aka the same format as the data in allele_counts_error_df)
comparison_allele_count <- do.call(rbind, replicate(n= nrow(rate_list)*length(reps), expr = allele_counts_df, simplify = FALSE))

# Pull out just the errored reads data frames but keep the nested list structure the same
genetic_data_sim_error <- lapply(rate_combo_count, function(x)
  lapply(reps, function(i) genetic_data_error[[x]][[i]][[1]]))

# Get the frequency of each allele from each locus from the errored datasets
allele_freqs_error <- lapply(genetic_data_sim_error, function(x)
  lapply(x, get_alleles_and_freqs))

# Get the number of alleles for each locus for each rep for each error rate combo and make into a data frame similar to allele_counts_df
allele_counts_error_df <- as.data.frame(do.call(rbind, lapply(rate_combo_count, function(x)
  as.data.frame(do.call(rbind, lapply(reps, function(i)
    as.data.frame(do.call(rbind, lapply(loci, function(z)
    return(lengths(allele_freqs_error[[x]][[i]][[z]])[1])))) %>%
      mutate(rep = i, 
             syst_err_rate = rate_list$syst_err_rate[x], 
             stoch_err_rate = rate_list$stoch_err_rate[x]))))))) %>%
  rename(allele_count = V1) %>%
  mutate(locus = rep(seq(1, num_loci, by = 1), nrow(rate_list)*length(reps))) 


# Calculate the proportion of loci for each unique replicate and error rate combo that have a different number of alleles in the error data compared to the original data
prop_loci_diff <- allele_counts_error_df %>%
  mutate(equal_to_orig = allele_count == comparison_allele_count) %>%
  group_by(across(all_of(c("rep", "syst_err_rate", "stoch_err_rate")))) %>%
  summarise(prop_diff = 1 - (sum(equal_to_orig)/num_loci))


# Make a boxplot of from all replicates within each unique error rate combo that shows the proportion of loci where the errored dataset had a different number of alleles than the non-errored dataset
prop_loci_diff %>%
  ggplot(aes(group = interaction(syst_err_rate, stoch_err_rate))) +
  geom_boxplot(aes(y = prop_diff, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  xlab("Stochastic error rate") +
  ylab("Proportion of loci where error caused change in allele #") +
  theme_classic()


# Get the proportion of the 3 least common alleles and identity of those alleles for all loci from the non-errored dataset (using custom get_df_of_sm_alleles function)
rarest_allele_prop_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_real, get_df_of_sm_alleles))) %>%
    mutate(locus = loci, 
           rep = "no error", 
           syst_err_rate = 0, 
           stoch_err_rate = 0)


# Get the proportion of the 3 least common alleles and identity of those alleles for all loci from all of the errored datasets 
rarest_allele_prop_df_error <- do.call(rbind, lapply(rate_combo_count, function(i) {
  single_rate_rarest_alleles <- do.call(rbind, lapply(reps, function(x) {
  # Make a data frame for each rep that consists of the info for the rarest alleles at each locus and notes the rep and locus
    single_rep_rarest_alleles <- as.data.frame(do.call(rbind, lapply(allele_freqs_error[[i]][[x]], get_df_of_sm_alleles))) %>%
      mutate(locus = loci, 
             rep = x, 
             syst_err_rate = rate_list$syst_err_rate[i], 
             stoch_err_rate = rate_list$stoch_err_rate[i])
    return(single_rep_rarest_alleles)
  }))
return(single_rate_rarest_alleles)
}))


# Merge the real and error rare allele dfs for visualization purposes
rarest_allele_prop_df %>% 
  rbind(rarest_allele_prop_df_error) %>%
  mutate(rep = as.factor(rep)) %>%
  # Plot the freq of the 3 rarest alleles w/ and without error
  ggplot(aes(x = locus)) +
  geom_point(aes(y = smallest_prop, color = rep), alpha = .4) + 
  geom_point(aes(y = sec_smallest_prop, color = rep), alpha = .4) + 
  geom_point(aes(y = third_smallest_prop, color = rep), alpha = .4) + 
  #scale_color_manual(values = c(rep("red", length(reps)), "black")) +
  facet_wrap(~ syst_err_rate * stoch_err_rate) +
  ylab("Frequency of alleles") +
  labs(title = "Frequency of 3 rarest alleles in the population with and without error") +
  theme_classic()


# Merge the real and error rare allele dfs for visualization purposes
rarest_allele_prop_df_error %>% 
  group_by(locus, syst_err_rate, stoch_err_rate) %>%
  summarise_at(vars("smallest_prop", "sec_smallest_prop", "third_smallest_prop"), median) %>%
  rbind(subset(rarest_allele_prop_df, select = -c(rep, smallest_allele_num, sec_smallest_allele_num, third_smallest_allele_num))) %>%
  # Plot the freq of the 3 rarest alleles w/ and without error
  ggplot(aes(x =  interaction(syst_err_rate, stoch_err_rate))) +
  geom_point(aes(y = smallest_prop, color =  interaction(syst_err_rate, stoch_err_rate))) + 
  geom_point(aes(y = sec_smallest_prop, color =  interaction(syst_err_rate, stoch_err_rate))) + 
  geom_point(aes(y = third_smallest_prop, color =  interaction(syst_err_rate, stoch_err_rate))) + 
  facet_wrap(~locus) +
  theme_classic() +
  scale_color_manual(values = c("red", rep("black", length(rate_combo_count)))) +
  ylab("Frequency of alleles") +
  labs(title = "Frequency of 3 rarest alleles in the population with and without error by locus") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())




#### Subsetting the initial dataset to get wild + garden dfs ####

# Percent of the wild population to pull for the garden df
perc_sampled <- .2

pops <- c(unique(genetic_data$Pop)) # Get the number of populations in the data_frame
genetic_data$location <- "wild" # Add a location column to the genetic data df

# Get a vector containing the IDs of the individuals which will be labelled as garden individuals
garden_inds <- unlist(lapply(pops, function(x) {
  inds <- genetic_data$Ind[genetic_data$Pop == x]
  garden <- sample(inds, size = perc_sampled * length(inds), replace = F) # Samples the specified percent of each pop
  return(garden)
}))

# Use the garden_inds vector to overwrite the location of the inds that have been pulled to be in the garden
genetic_data[genetic_data$Ind %in% garden_inds, c("location")] <- "garden"

save(genetic_data,file=paste0(sim.wd, "RObjects/genetic_data.RData"))

# Making the all of error datasets have all the info present in the real dataset but keep the data in the nested lists format
genetic_data_error <- lapply(genetic_data_sim_error, function(i) {
  all_genetic_data_reps <- lapply(i, function(x) {
  single_genetic_data_error_rep <- x %>%
  mutate(location = genetic_data$location, 
         Ind = genetic_data$Ind, 
         Pop = genetic_data$Pop)
  return(single_genetic_data_error_rep)})
  return(all_genetic_data_reps)})

save(genetic_data_error,file=paste0(sim.wd, "RObjects/genetic_data_error.RData"))

#### Comparing the wild data and the garden data (allele count, allele frequency in wild, heterozygosity at each allele in population) across all datasets) ####


# Make lists of all alleles and frequencies from the error and non-errored data 
w_g_real_alleles_and_freqs <- get_alleles_and_freqs_w_g(genetic_data)
w_g_error_alleles_and_freqs <- lapply(genetic_data_error, function(x)
  lapply(x, get_alleles_and_freqs_w_g))


# Make the full data frames that hold all summary info across all alleles/loci
full_dfs_in_lists <- future_lapply(w_g_error_alleles_and_freqs, future.seed=TRUE, function(x) {
  lapply(x, function(i) {
  lapply(loci, compare_w_g, w_g_real_alleles_and_freqs, i)})})


# Make data frame with all of the info across all rate combos and replicates at the allele level
allele_info_df <- do.call(rbind, lapply(rate_combo_count, function(z){
  all_allele_info_per_rate_combo <- do.call(rbind, lapply(reps, function(i){
    all_allele_info_per_rep <- do.call(rbind, lapply(loci, function(x) {
    allele_info_per_rep <- full_dfs_in_lists[[z]][[i]][[x]][[1]] %>%
      mutate(rep = i, 
             syst_err_rate = rate_list$syst_err_rate[z], 
             stoch_err_rate = rate_list$stoch_err_rate[z]) 
    return(allele_info_per_rep)
    }))
  return(all_allele_info_per_rep)
  }))
return(all_allele_info_per_rate_combo)
}))


# Make the allele dataset nicer for graphing
allele_info_df_for_graphing <- allele_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error"), 
         syst_err_rate = ifelse(dataset %in% c("wild_real", "garden_real"), 0, syst_err_rate), 
         stoch_err_rate = ifelse(dataset %in% c("wild_real", "garden_real"), 0, stoch_err_rate),
         dataset_rates = ifelse(dataset %notin% c("wild_error" , "garden_error"), 
                                paste0(dataset),
                                paste0(dataset, "_", syst_err_rate, "_", stoch_err_rate)), 
         w_g = ifelse(dataset %in% c("garden_error", "garden_real"), "g", "w"), 
         rep= as.factor(rep)) %>%
  mutate(dataset_rates = fct_relevel(dataset_rates, "garden_real", "wild_real", after = 0))

save(allele_info_df_for_graphing,file=paste0(sim.wd, "RObjects/allele_info_df_for_graphing.RData"))


# Make data frame with all of the info across all rate combos and replicates at the locus level
loci_info_df <- do.call(rbind, lapply(rate_combo_count, function(z){
  all_loci_info_per_rate_combo <- do.call(rbind, lapply(reps, function(i){
    all_loci_info_per_rep <- do.call(rbind, lapply(loci, function(x) {
      loci_info_per_rep <- full_dfs_in_lists[[z]][[i]][[x]][[2]] %>%
        mutate(rep = i, 
               syst_err_rate = rate_list$syst_err_rate[z], 
               stoch_err_rate = rate_list$stoch_err_rate[z])
      return(loci_info_per_rep)
    }))
    return(all_loci_info_per_rep)
  }))
  return(all_loci_info_per_rate_combo)
}))

# Make the loci dataset nicer for graphing
loci_info_df_for_graphing <- loci_info_df%>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error"), 
         syst_err_rate = ifelse(dataset %in% c("wild_real", "garden_real"), 0, syst_err_rate), 
         stoch_err_rate = ifelse(dataset %in% c("wild_real", "garden_real"), 0, stoch_err_rate),
         dataset_rates = ifelse(dataset %notin% c("wild_error" , "garden_error"), 
                                paste0(dataset),
                                paste0(dataset, "_", syst_err_rate, "_", stoch_err_rate)), 
         w_g = ifelse(dataset %in% c("garden_error", "garden_real"), "g", "w"), 
         rep= as.factor(rep)) %>%
  mutate(dataset_rates = fct_relevel(dataset_rates, "garden_real", "wild_real", after = 0)) %>%
  distinct(.keep_all = TRUE)

save(loci_info_df_for_graphing,file=paste0(sim.wd, "RObjects/loci_info_df_for_graphing.RData"))




#### Resampling all wild datasets and comparing error to non-error data

# Resample the non-errored data 100 times
resampling_reps = 100

resample_df <- resampling_diversity(ind_data = genetic_data, reps = resampling_reps)  %>%
 mutate(error_rep = NA, 
        syst_err_rate = as.factor(0), 
        stoch_err_rate = as.factor(0)) %>%
  filter(locus == "all") # Line to hopefully reduce the amount of time/effort and aborted sessions since Sean says he doesn't feel like we need to look at all loci


# Resample the errored data sets 20 times (since 5 reps of added error this = 100 resampling reps per error rate combo)
resampling_reps = 20


resample_df_error_list <- future_lapply(genetic_data_error, future.seed=TRUE, function(i){
  lapply(i, resampling_diversity, reps = resampling_reps)
})
save(resample_df_error_list,file=paste0(sim.wd, "RObjects/20_resamples_df_error.RData"))


# Bind the errored data together to make one huge dataframe
resample_df_error <- do.call(rbind, future_lapply(rate_combo_count, future.seed=TRUE, function(z){
  resample_df_info_per_rate_combo <- do.call(rbind, lapply(reps, function(i){
      resample_df_info_per_rep <- resample_df_error_list[[z]][[i]] %>%
        mutate(error_rep = i, 
               syst_err_rate = as.factor(rate_list$syst_err_rate[z]), 
               stoch_err_rate = as.factor(rate_list$stoch_err_rate[z])) %>%
        filter(locus == "all") # Line to hopefully reduce the amount of time/effort and aborted sessions since Sean says he doesn't feel like we need to look at all loci
      return(resample_df_info_per_rep)
  }))
  return(resample_df_info_per_rate_combo)
}))

resample_df_full <- rbind(resample_df, resample_df_error)
save(resample_df_full,file=paste0(sim.wd, "RObjects/100_resamples_all_datasets_df.RData"))
