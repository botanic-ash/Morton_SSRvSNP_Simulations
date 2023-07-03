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


source("/Users/Ashley/Desktop/Fall\ 2019/stats/diagPlots.R") # This is from my stats prof at TAMUCC, it provides the SW p value to tell if things are significantly non-normal

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
#sim.wd <- "/home/akoontz/Shared/SSRvSNP_Sim/Code/"
sim.wd <- "/Users/Ashley/Desktop/Hoban/Hoban_rotation_2023/Morton_SSRvSNP_Simulations/"
arq.wd <- "/Users/Ashley/Desktop/Hoban/Hoban_rotation_2023/Morton_SSRvSNP_Simulations/SimulationOutputs/MSAT_marker/MSAT_ArlequinFiles"
setwd(sim.wd)

# Parallelism: specify number of cores to use; NOT USED RIGHT NOW
# num_cores <- detectCores() - 8
# num_cores <- 4


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
syst_err_rate = c(.01, .05, .1, .2) # Based off of rates of errored loci found in Austin's real SNP datasets where source of error is unknown
stoch_err_rate = c(.01, .05, .1, .2)

# Make a list of all unique error rate combinations
rate_list <- expand.grid(syst_err_rate, stoch_err_rate) %>%
  rename(syst_err_rate = Var1, 
         stoch_err_rate = Var2)

rate_combo_count <- seq(1, nrow(rate_list), 1) # Make a list corresponding to the number of unique error rate combinations

# Set number of reps per unique error rate combination
rep_num = 5 
reps = as.list(seq(1, rep_num, by = 1)) # Make a list corresponding to the number of reps per error rate combination


# Below code is pretty clunky and slow but it get's the job done, could be parallelized fairly easily though I think (with mclapply)

# Apply the simulating_error function to the input genetic data across the list of unique replicates and then across the list of unique rate value combinations
# Output is a list of lists of lists (3 nested lists):
#   1st list contains all unique error rate combos
#   2nd list contains each replicate error dataset with the given error rate values
#   3rd list contains 2 dataframes: the error dataset and the tracked error dataset
genetic_data_error <- lapply(rate_combo_count, function(x)
  lapply(reps, function(i)
       simulating_error(data_for_editing = genetic_data, syst_err_rate = rate_list$syst_err_rate[x], stoch_err_rate = rate_list$stoch_err_rate[x], msat_length = 1, replicate_number = i)))


#### Summarizing error added to the data ####
# Make a dataframe from the tracked error dataset from each rep from each unique rate value combo and clean the dataframe  
genetic_data_track_error <- as.data.frame(do.call(rbind, lapply(genetic_data_error, function(x)
  do.call(rbind, lapply(x, function(i) 
  colSums(i[[2]])/(nrow(i[[2]])*2)))))) %>% # Calc total error rate at each locus
  # Need to multiply the last few columns by 2 bc divided by 2 x the number of rows to calculate the true loci error
  mutate(rep = as.factor(rep *2), 
         syst_err_rate = as.factor(syst_err_rate *2), 
         stoch_err_rate = as.factor(stoch_err_rate *2)) %>%
  pivot_longer(cols = starts_with("locus"), names_to = "locus", values_to = "error_rate") # Make it so each locus is a unique row instead of a unique column

#Plot the summarized error data at each locus: faceted by unique error rate combo, replicates condensed into boplots, all 20 loci on the x axis
genetic_data_track_error %>%
  ggplot() +
  geom_hline(aes(yintercept = median(error_rate)), alpha = .3) +
  geom_boxplot(aes(x = locus, y = error_rate), outlier.shape = NA) + 
  #geom_point(aes(x = locus, y = error_rate, color = rep)) +
  facet_wrap(~ syst_err_rate * stoch_err_rate) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Plot the summarized error data at each locus: one boxplot represents all replicates at all loci at a single unique error rate combo
genetic_data_track_error %>%
  ggplot() +
  geom_hline(aes(yintercept = median(error_rate)), alpha = .3) +
  geom_boxplot(aes(y = error_rate, group = interaction(syst_err_rate, stoch_err_rate), color = interaction(syst_err_rate, stoch_err_rate))) +
  ylab("real error rate") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())


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
  geom_boxplot(aes(y = prop_diff, x = interaction(syst_err_rate, stoch_err_rate), color = interaction(syst_err_rate, stoch_err_rate))) +
  #scale_color_viridis(discrete = TRUE) +
  xlab("type 1 rate x type 2 rate") +
  ylab("Proportion of loci where error caused change in allele #") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #facet_wrap(~ syst_err_rate * stoch_err_rate)


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


#### Want to make a graph that shows the number of reps that have changed the identity of the 1st 2nd and 3rd least frequent allele !####



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

# Making the all of error datasets have all the info present in the real dataset but keep the data in the nested lists format
genetic_data_error <- lapply(genetic_data_sim_error, function(i) {
  all_genetic_data_reps <- lapply(i, function(x) {
  single_genetic_data_error_rep <- x %>%
  mutate(location = genetic_data$location, 
         Ind = genetic_data$Ind, 
         Pop = genetic_data$Pop)
  return(single_genetic_data_error_rep)})
  return(all_genetic_data_reps)})


#### Comparing the wild data and the garden data (allele count, allele frequency in wild, heterozygosity at each allele in population) across all datasets) ####


# Make lists of all alleles and frequencies from the error and non-errored data 
w_g_real_alleles_and_freqs <- get_alleles_and_freqs_w_g(genetic_data)
w_g_error_alleles_and_freqs <- lapply(genetic_data_error, function(x)
  lapply(x, get_alleles_and_freqs_w_g))


# Make the full data frames that hold all summary info across all alleles/loci
full_dfs_in_lists <- lapply(w_g_error_alleles_and_freqs, function(x) {
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




#### Making some graphs####

# Plot of a single locus
plot_locus = 1

# Plot of all alleles of a single locus at each error rate combo
allele_info_df_for_graphing  %>%
  group_by(all_alleles, locus, syst_err_rate, stoch_err_rate, dataset) %>%
  summarise_at(vars("allele_freq"), median) %>% # Summarize at median of replicates at same error rates
  filter(locus == plot_locus) %>%
  ggplot(alpha = .7) +
  geom_point(aes(x = all_alleles, y = allele_freq, color = dataset), shape = 1) +
  scale_color_manual(values = c( "black", "gray",  "red", "indianred")) +
  facet_wrap(~ syst_err_rate * stoch_err_rate) +
  ylab("Allele frequency in the given dataset") +
  xlab("Allele identity") +
  theme_classic()


# Boxplot of frequencies of all alleles across all loci at each error rate combo
allele_info_df_for_graphing  %>%
  mutate(w_g = fct_relevel(w_g, "w", "g")) %>%
  ggplot() +
  #geom_point(aes(y = allele_freq, x = w_g, color = dataset), alpha = .1) +
  geom_boxplot(aes(y = allele_freq, x = w_g, color = dataset),outlier.shape = NA) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_grid(~ syst_err_rate * stoch_err_rate) +
  theme_classic()



# Needed for graph below to ensure colors follow the correct order since colors coming from 2 different columns I can't relevel to make them correct
my_colors <- c("garden_real" = "gray", 
               "wild_real" = "black",
               "garden_error" = "indianred",
               "wild_error" = "red", 
               "TRUE" = "gold", 
               "FALSE" = "brown")

# Plot differences between expected and observed het for each dataset split by each locus (all rate combos on x axis), color of connecting line shows whether obs het (brown) or exp het (gold) is higher.. not overly important
loci_info_df_for_graphing  %>%
  group_by(locus, syst_err_rate, stoch_err_rate, dataset, dataset_rates) %>%
  summarise_at(vars("exp_het", "obs_het"), median) %>%
  mutate(line_color = exp_het >= obs_het) %>%
  ggplot() +
  geom_point(aes(y = obs_het, x = dataset_rates, color = dataset), shape = 19) + # Obs het = closed circle
  geom_point(aes(y = exp_het, x = dataset_rates, color = dataset), shape = 1) + # Exp het = closed circle
  geom_segment(aes(y = obs_het, yend = exp_het, x = dataset_rates, xend= dataset_rates, color = line_color)) +
  scale_color_manual(values = c(my_colors)) +
  facet_wrap(~ locus) +
  theme_classic()+
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  ylab("Heterozygosity estimates") +
  xlab("Systematic error rate * stochastic error rate")
#Interpretation: at lower systematic error rates, expected het is often higher than obs het but never by very much and as the systematic rate increases the difference between the expect and observed het gets larger


# Plot of the comparison of observed heterozygosity of all loci in each dataset 
loci_info_df_for_graphing  %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = as.factor(syst_err_rate), color = as.factor(stoch_err_rate))) +
  facet_wrap(~ w_g) +
  theme_classic()
#Interpretation: at higher stochastic error rates, obs het increases.... makes sense bc literally making more heterozygotes, at higher systematic error rates, obs het decreases makes sense bc removing hets


loci_info_df  %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Observed heterozygosity") +
  xlab("Stochastic error rate") +
  theme_classic()
#Interpretation: much more variance in obs het than in the expected het estimates, at higher systematic error rates a bit lower obs het, at higher stochastic error erates higher obs het


# Plot of the comparison of expected heterozygosity of all loci in each dataset
loci_info_df_for_graphing  %>%
  mutate(w_g = fct_relevel(w_g, "w", "g")) %>%
  ggplot() +
  geom_boxplot(aes(y = exp_het, x = w_g, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_grid(~ syst_err_rate * stoch_err_rate) +
  theme_classic()
#Interpretation: at higher stochastic error rates, exp het increases, there doesn't seem to be much effect of systematic error on expected heterozygosity


# Plot expected het from the wild datasets of all combos of error rates against next to the single real
loci_info_df  %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = exp_het, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Expected heterozygosity") +
  xlab("Stochastic error rate") +
  theme_classic()
#Interpretation: at higher stochastic error rates, exp het increases, there doesn't seem to be much effect of systematic error on expected heterozygosity


# Plot of the comparison of inbreeding coefficient of all loci in each dataset 
loci_info_df_for_graphing  %>%
  mutate(w_g = fct_relevel(w_g, "w", "g")) %>%
  ggplot() +
  geom_boxplot(aes(y = inbreeding_coeff, x = w_g, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_grid(~ syst_err_rate * stoch_err_rate) +
  theme_classic()


# Plot inbreeding coefficient from the wild datasets of all combos of error rates against next to the single real
loci_info_df  %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = inbreeding_coeff, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Inbreeding coefficient (F)") +
  xlab("Stochastic error rate") +
  theme_classic()
#Interpretation: at higher stochastic error rates, inbreeding coeff increases, but at higher systematc error rates inbreeding coeff decreases








# Comparison of number of across all loci in each dataset 
loci_info_df_for_graphing  %>%
  ggplot() +
  geom_boxplot(aes(y = num_alleles, x = interaction(syst_err_rate, stoch_err_rate), color = as.factor(stoch_err_rate)), outlier.shape = NA) +
  geom_jitter(aes(y = num_alleles, x = interaction(syst_err_rate, stoch_err_rate), color = as.factor(stoch_err_rate)), alpha = .3) +
  facet_wrap(~ w_g) +
  ylab("Total number of alleles across all loci") +
  xlab("Systematic error rate * stochastic error rate") +
  labs(colour = "Stochastic error rate") +
  theme_classic()
# Interpretation: wild dataset has consistently more alleles than garden


# Plot of the proportion of wild rare alleles captured in garden samples in error vs real dataset across all loci,
loci_info_df_for_graphing %>%
  mutate(stoch_err_rate = as.factor(stoch_err_rate)) %>%
  ggplot() +
  geom_boxplot(aes(y = prop_rare_wild_alleles_captured, x = stoch_err_rate, color = interaction(syst_err_rate, stoch_err_rate)), outlier.shape = NA) +
  # geom_jitter(aes(y = prop_rare_wild_alleles_captured, x = stoch_err_rate, color = interaction(syst_err_rate, stoch_err_rate)), alpha = .3) + # add points onto the boxplots
  ylab("Proportion of \"rare\" wild alleles in the garden sample") +
  xlab("stochastic error rate") +
  labs(colour = "Systematic error rate * stochastic error rate") +
  theme_classic()
# Interpretation: rare loci are less likely to be represented in the garden subsample in the error dataset


# Plot comparing the number of rare alleles across all loci in the wild dataset without error and the wild datasets with increasing amounts of error
loci_info_df_for_graphing %>%
  filter(dataset %in% c("wild_real", "wild_error")) %>%
  group_by(locus, syst_err_rate, stoch_err_rate, dataset, dataset_rates) %>%
  summarise_at(vars("rare_allele_count"), median) %>%
  ggplot(aes(x = rare_allele_count)) +
  geom_bar(aes(fill = dataset), position=position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("black","red")) +
  ylab("Number of loci with given amount of \"rare\" alleles") +
  xlab("Number of \"rare\" alleles") +
  theme_classic() +
  facet_wrap(~ syst_err_rate * stoch_err_rate, scales = 'free_x') +
  scale_x_continuous(limits = c(min(loci_info_df_for_graphing$rare_allele_count - 1), max(loci_info_df_for_graphing$rare_allele_count + 1)))
# Interpretation: as you increase error, there are more loci with rare alleles


# Plot comparing the number of rare alleles across all loci in the garden dataset without error and the garden datasets with increasing amounts of error
loci_info_df_for_graphing %>%
  filter(dataset %in% c("garden_real", "garden_error")) %>%
  group_by(locus, syst_err_rate, stoch_err_rate, dataset, dataset_rates) %>%
  summarise_at(vars("rare_allele_count"), median) %>%
  ggplot(aes(x = rare_allele_count)) +
  geom_bar(aes(fill = dataset), position=position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("gray","indianred")) +
  ylab("Number of loci with given amount of \"rare\" alleles") +
  xlab("Number of \"rare\" alleles") +
  theme_classic() +
  facet_wrap(~ syst_err_rate * stoch_err_rate, scales = 'free_x') +
  scale_x_continuous(limits = c(min(loci_info_df_for_graphing$rare_allele_count - 1), max(loci_info_df_for_graphing$rare_allele_count + 1)))
# Interpretation: as you increase error, there are more rare loci


# Plot comparing the number of rare alleles per locus all loci in the garden dataset without error and the garden datasets with increasing amounts of error
loci_info_df_for_graphing %>%
  ggplot(aes(y = rare_allele_count)) +
  geom_boxplot(aes(x = as.factor(stoch_err_rate), color = interaction(syst_err_rate, stoch_err_rate)), outlier.shape = NA) +
  #geom_jitter(aes(x = as.factor(stoch_err_rate), color = interaction(syst_err_rate, stoch_err_rate)), alpha = .3) +
  ylab("Number of \"rare\" alleles per locus") +
  xlab("Stochastic error rate") +
  labs(colour = "Systematic error rate * stochastic error rate") +
  theme_classic() +
  facet_wrap(~ w_g)
# Interpretation: there are less rare alleles in the garden subset than the wild dataset (as expected)



#### Statistical tests####

### Analyzing differences in allele number and rarity ###

## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the number of alleles per locus 
allele_number_aov_full <- aov(num_alleles ~ dataset + locus + syst_err_rate + stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(allele_number_aov_full), main = "Residuals") 
# Interpretation:  a little widdly but doens't seem like too big of a deal
dredge(allele_number_aov_full, rank = "AICc") 
# Interpretation: 2 best models include all but one of the variables (systematic error rate) or all of the variables

allele_number_top_model <- summary(lm(data = loci_info_df_for_graphing, num_alleles ~ dataset + locus + syst_err_rate + stoch_err_rate))
allele_number_top_model
# Interpretation: significant effect of all dataset types except garden_error, presumably because most of the differences are also explain by garden_real?, significant effect of locus and stochastic error rate

allele_number_second_model <- summary(lm(data = loci_info_df_for_graphing, num_alleles ~ dataset + locus + stoch_err_rate))
allele_number_second_model
# Interpretation: almost the exact same results as the above model except garden_error has a very slght negative estimate instead of a very slight positive estimate

#Summary: the garden_real is significantly correlated to less alleles per locus, the wild_real and the wild_real datasets are significantly correlated to more alleles per locus, the garden_error dataset is not a significant explanatory variable and is not consistently correlated to more or less alleles per locus, the stochastic error rate is significantly correlated to an increase in the number of alleles after accounting for variability due to dataset type, the systematic error rate is not a significant explanatory variable but is correlated to less rare alleles per locus


## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the number of rare alleles per locus 
rare_allele_number_aov_full <- aov(rare_allele_count ~ dataset + locus + syst_err_rate + stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(rare_allele_number_aov_full), main = "Residuals") 
#Interpretation:  pretty non normal I think bc it's not really a continuous indep variable?
dredge(rare_allele_number_aov_full, rank = "AICc") 
#Interpretation: first 4 models are all very similar, all include the dataset and stochastic error rate (2nd best includes only those 2 variables), will look at the best 2 and then decide if I want to see the other 2 best ones

rare_allele_number_top_model <- summary(lm(data = loci_info_df_for_graphing, rare_allele_count ~ dataset + syst_err_rate + stoch_err_rate))
rare_allele_number_top_model
# Interpretation: significant effect of all dataset types, significant effect of stochastic error rate

rare_allele_number_second_model <- summary(lm(data = loci_info_df_for_graphing, rare_allele_count ~ dataset + stoch_err_rate))
rare_allele_number_second_model
# Interpretation: all significant variables as the top model and nearly identical estimates

# Summary: the garden_real and garden_error datasets are significantly correlated to less rare alleles per locus, the wild_real and the wild_real datasets are significantly correlated to more rare alleles per locus, the stochastic error rate is significantly correlated to an increase in the number of rare alleles after accounting for variability due to dataset type, the systematic error rate is not a significant explanatory variable but is correlated to less rare alleles per locus, R2 around .20


## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the proportion of rare wild alleles captured per locus
prop_rare_alleles_captured_aov_full <- aov(prop_rare_wild_alleles_captured ~ dataset + locus + syst_err_rate + stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(prop_rare_alleles_captured_aov_full), main = "Residuals") 
#Interpretation: wildly non-normal presumably because all values are between 0 and 1 and most seem to be 0 or 1 --> probably need to do a beta regression?
dredge(prop_rare_alleles_captured_aov_full, rank = "AICc") 
#Interpretation: top 2 models are fairly comparable and include either all factors or all but systematic error rate 

prop_rare_alleles_captured_top_model <- summary(lm(data = loci_info_df_for_graphing, prop_rare_wild_alleles_captured ~ dataset + locus + stoch_err_rate))
prop_rare_alleles_captured_top_model
# Interpretation: significant effect of all included explanatory variables, note: since this metric is a ratio of garden to wild datasets, this is looking at the difference between real and error datasets

prop_rare_alleles_captured_second_model <- summary(lm(data = loci_info_df_for_graphing, prop_rare_wild_alleles_captured ~ dataset + locus + syst_err_rate + stoch_err_rate))
prop_rare_alleles_captured_second_model
# Interpretation: same significant variables as the top model (aka systematic error is not significant) and nearly identical estimates

# Summary: the real dataset is significantly correlated to a higher proportion of rare wild alleles being captured across all loci, the error dataset is significantly correlated to more a lower proportion of rare wild alleles being captured across all loci, there is a significant effect of locus on the number of rare wild alleles being captured, the stochastic error rate is significantly correlated to higher proportion of rare wild alleles being captured across all loci, the systematic error rate is not a significant explanatory variable but is correlated to lower proportion of rare wild alleles being captured across all loci


# Since the anova residuals were so non-normal, trying a beta regression on a z transformed version of prop rare wild alleles

# Make the transformed column
trans_loci_info_df_for_graphing <- loci_info_df_for_graphing %>%
  filter(dataset %notin% c("garden_real", "garden_error"))  %>% # Remove columns with na so that the betregression won't fail
  filter(!is.nan(prop_rare_wild_alleles_captured))  %>%
  mutate(trans_prop_rare_wild_alleles_captured = (prop_rare_wild_alleles_captured*(nrow(.) - 1) + .5)/nrow(.))


betareg_prop_rare_wild_alleles_captured_full <- betareg(trans_prop_rare_wild_alleles_captured ~ dataset + locus + syst_err_rate + stoch_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
qqnorm.with.sim.bounds(resid(betareg_prop_rare_wild_alleles_captured_full), main = "Residuals")
#Interpretation: Still not normal at all...
dredge(betareg_prop_rare_wild_alleles_captured_full, rank = "AICc") 
#Interpretation: top 4 models are best and all within 2 AICc points, fairly comparable and include either all factors or all but systematic error rate, best 2 don't include systematic error rate, 2nd and 4th don't include locus

betareg_prop_rare_wild_alleles_captured_top_model <- betareg(trans_prop_rare_wild_alleles_captured ~ dataset + locus + stoch_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
summary(betareg_prop_rare_wild_alleles_captured_top_model)
# Interpretation: all included explanatory variables except locus are significant, pretty low R2

betareg_prop_rare_wild_alleles_captured_second <- betareg(trans_prop_rare_wild_alleles_captured ~ dataset + stoch_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
summary(betareg_prop_rare_wild_alleles_captured_second)
# Interpretation: all included explanatory variables are significant, pretty low R2


# Summary: similar to the non beta-regression --> the real dataset is significantly correlated to a higher proportion of rare wild alleles being captured across all loci, the error dataset is significantly correlated to a lower proportion of rare wild alleles being captured across all loci, and  the stochastic error rate is significantly correlated to higher proportion of rare wild alleles being captured across all loci while the systematic error rate is not a significant explanatory variable. However unlike the non beta-regression, there is not a significant effect of locus on the number of rare wild alleles being captured



### Analyzing differences in heterozygosity ###

## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the expected heterozygosity
expected_het_aov_full <- aov(exp_het ~ dataset + locus + syst_err_rate + stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(expected_het_aov_full), main = "Residuals") 
#Interpretation: verry close to normal!
dredge(expected_het_aov_full, rank = "AICc") 
#Interpretation: only one good model, has all variables

expected_het_top_model <- summary(lm(data = loci_info_df_for_graphing, exp_het ~ dataset + locus + syst_err_rate + stoch_err_rate))
expected_het_top_model
# Interpretation: significant effect of both the wild_real and wild_error datasets where wild_error dataset correlates with a significantly lower expected heterozygosity, both error types have a significant effect on expected heterozygosity

# Summary: the wild datasets are both significantly correlated eith the expected heterozygosity, but the wild_error dataset is significantly correlated to a lower expected heterozygosity, while the real dataset is significantly correlated to a higher expected heterozygosity, the stochastic error rate is significantly correlated to a higher expected heterozygosity while the systematic error rate is significantly correlated to a lower expected heterozygosity


## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the observed heterozygosity
observed_het_aov_full <- aov(obs_het ~ dataset + locus + syst_err_rate + stoch_err_rate, data = loci_info_df_for_graphing) 
qqnorm.with.sim.bounds(resid(observed_het_aov_full), main = "Residuals") 
#Interpretation: vnot super normal...
dredge(observed_het_aov_full, rank = "AICc") 
#Interpretation: 2 top models and neither include dataset type, the top model also doens't include the systematic error rate


observed_het_top_model <- summary(lm(data = loci_info_df_for_graphing, obs_het ~ locus + stoch_err_rate))
observed_het_top_model
# Interpretation: significant effect of all included explanatory variables

observed_het_second_model <- summary(lm(data = loci_info_df_for_graphing, obs_het ~ locus + syst_err_rate + stoch_err_rate))
observed_het_second_model
# Interpretation: same significant variables as the top model (aka systematic error is not significant) and nearly identical estimates


# Summary: datasets type is not important in the explanation of differences in observed heterozygosity, the stochastic error rate is significantly correlated to a higher observed heterozygosity while the systematic error rate is not a significant explanatory variable but does correlate to a slightly higher observed heterozygosity


## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the inbreeding coefficient (F)
inbreeding_coeff_aov_full <- aov(inbreeding_coeff ~ dataset + locus + syst_err_rate + stoch_err_rate, data = loci_info_df_for_graphing) 
qqnorm.with.sim.bounds(resid(inbreeding_coeff_aov_full), main = "Residuals") 
#Interpretation: looks pretty good!!
dredge(inbreeding_coeff_aov_full, rank = "AICc") 
#Interpretation: 2 top models, one includes all explanatory variables and the other includes all but locus


inbreeding_coeff_top_model <- summary(lm(data = loci_info_df_for_graphing, inbreeding_coeff ~ dataset + syst_err_rate+ stoch_err_rate))
inbreeding_coeff_top_model
# Interpretation: significant effect of all included explanatory variables except the wild_real dataset

inbreeding_coeff_second_model <- summary(lm(data = loci_info_df_for_graphing, inbreeding_coeff ~ dataset + locus + syst_err_rate + stoch_err_rate))
inbreeding_coeff_second_model
# Interpretation: same significant variables as the top model (aka locus is not significant) and nearly identical estimates


# Summary: datasets type is important in the explanation of differences in the inbreeding coefficient with the real datasets correlating to lower F values (significant only for the garden_real dataset while both the wild and garden error datasets significantly correlate to slightly higher F values, the stochastic error rate is significantly correlated to a lower F value while the systematic error rate is significantly correlated to a higher F value


# Questions for Sean:
# Do I want to add any other interaction terms? 
# Is it not okay to be treating my error values as continous when they are discrete (atm)?
  # Is it worth making more error values so these variables are closer to continuous?
