#### Packages and objects to load ####
library(ggplot2)
library(forcats)
library(ggnewscale)
library(MuMIn)
library(betareg)
library(parallelly)
library(future.apply)
library(tidyverse)


loci_info_df_for_graphing <- load(file=paste0(sim.wd, "RObjects/loci_info_df_for_graphing.RData"))
allele_info_df_for_graphing <- load(file=paste0(sim.wd, "RObjects/allele_info_df_for_graphing.RData"))
resample_df_full <- load(file=paste0(sim.wd, "RObjects/100_resamples_all_datasets_df.RData"))

#### Observed heterozygosity ####

# Plot of the comparison of observed heterozygosity of all loci in each dataset 
loci_info_df_for_graphing  %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = as.factor(syst_err_rate), color = as.factor(stoch_err_rate))) +
  facet_wrap(~ w_g) +
  theme_classic()
#Interpretation: at higher stochastic error rates, obs het increases.... makes sense bc literally making more heterozygotes, at higher systematic error rates, obs het decreases makes sense bc removing hets


loci_info_df_for_graphing  %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Observed heterozygosity") +
  xlab("Stochastic error rate") +
  theme_classic()
#Interpretation: much more variance in obs het than in the expected het estimates, lower obs het at higher systematic error rates, higher obs het at higher stochastic error erates--> both error types interact to effect obs het
#ggsave(paste0(sim.wd, "/Figs/obs_het_boxplots.png"), width = 7, height = 3)


### ANOVA(s) ###
## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the observed heterozygosity
observed_het_aov_full <- aov(obs_het ~ w_g + syst_err_rate + stoch_err_rate + syst_err_rate:stoch_err_rate, data = loci_info_df_for_graphing) 
qqnorm.with.sim.bounds(resid(observed_het_aov_full), main = "Residuals") 
#Interpretation: not terrible
dredge(observed_het_aov_full, rank = "AICc") 
#Interpretation: only 1 top model with all explanatory variables


observed_het_top_model <- summary(lm(data = loci_info_df_for_graphing, obs_het ~ w_g + syst_err_rate * stoch_err_rate ))
observed_het_top_model
# Interpretation: significant effect of all included explanatory variables


# Summary: wild datasets have slightly higher obs heterozygosity than garden datasets (?) differences in observed heterozygosity, the stochastic error rate is significantly correlated to a higher observed heterozygosity while the systematic error rate is significantly  correlated to a lower observed heterozygosity, and the interaction of systematic and stochastic error rates is significantly correlated to a higher obs het, R2 of .38



#### Expected heterozygosity ####
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
loci_info_df_for_graphing  %>%
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


### ANOVA(s) ###

## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the expected heterozygosity
expected_het_aov_full <- aov(exp_het ~ w_g + syst_err_rate + stoch_err_rate + syst_err_rate:stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(expected_het_aov_full), main = "Residuals") 
#Interpretation: decently normal
dredge(expected_het_aov_full, rank = "AICc") 
#Interpretation: 4 top models, best 2 include just stochastic error rate and stochastic error rate and w_g

expected_het_top_model <- summary(lm(data = loci_info_df_for_graphing, exp_het ~ stoch_err_rate))
expected_het_top_model
# Interpretation: significant positive effect of stochastic error

# NEED TO RUN THE OTHER MODELS

# Summary: the wild datasets are both significantly correlated with the expected heterozygosity, but the wild_error dataset is significantly correlated to a lower expected heterozygosity, while the real dataset is significantly correlated to a higher expected heterozygosity, the stochastic error rate is significantly correlated to a higher expected heterozygosity while the systematic error rate is significantly correlated to a lower expected heterozygosity


#### Inbreeding coefficient ####

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


### ANOVA(s) ###

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



#### Number of alleles ####

# Comparison of number of alleles across all loci in each dataset 
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



# Compare number of alleles across all loci only in wild
loci_info_df_for_graphing  %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = num_alleles, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Total number of alleles per locus") +
  xlab("Stochastic error rate") +
  theme_classic()
ggsave(paste0(sim.wd, "/Figs/num_alleles_boxplots.png"), width = 7, height = 3)



### ANOVA(s) ###

# Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the number of alleles per locus 

allele_number_aov_full <- aov(num_alleles ~ w_g + syst_err_rate + stoch_err_rate + syst_err_rate:stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(allele_number_aov_full), main = "Residuals") # Unclear why this is failing... probably bc this is a 2 way ANOVA?
# Interpretation:  a little wiggly but doesn't seem like too big of a deal
dredge(allele_number_aov_full, rank = "AICc") 
# Interpretation: 3 best models, top 2 don't include the interaction term, top doesn't include systematic error rate

allele_number_top_model <- summary(lm(data = loci_info_df_for_graphing, num_alleles ~ w_g + stoch_err_rate))
allele_number_top_model
# Interpretation: significant effect of all included variables, number of alleles per locus in wild data is slightly higher than in the garden data, stochastic error rate is significantly correlated to a larger number of alleles per locus

allele_number_second_model <- summary(lm(data = loci_info_df_for_graphing, num_alleles ~ w_g + stoch_err_rate + syst_err_rate))
allele_number_second_model
# Interpretation: almost the exact same results as the above (systematic error rate does not have a significant effect on number of alleles per locus)

allele_number_third_model <- summary(lm(data = loci_info_df_for_graphing, num_alleles ~ w_g + stoch_err_rate * syst_err_rate))
allele_number_third_model
# Interpretation: almost the exact same results as the above (interaction of error types does not have a significant effect on number of alleles per locus)

#Summary: dataset type (wild or garden) is significantly correlated to the number of alleles per locus with the wild datasets having slightly more alleles per locus than the garden datasets, the stochastic error rate is significantly correlated to an increase in the number of alleles after accounting for variability due to dataset (wild or garden), neither the systematic error rate or the interaction of the two error types are significant explanatory variables, R2 around .11


#### Number of rare alleles ####

# Plot comparing the number of rare alleles per locus all loci in the garden dataset without error and the garden datasets with increasing amounts of error
loci_info_df_for_graphing %>%
  ggplot(aes(y = rare_allele_count)) +
  geom_boxplot(aes(x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate)), outlier.shape = NA) +
  ylab("Number of \"rare\" alleles per locus") +
  xlab("Stochastic error rate") +
  labs(colour = "Systematic error rate") +
  theme_classic() +
  facet_wrap(~ w_g)
# Interpretation: there are less rare alleles in the garden subset than the wild dataset (as expected)


### ANOVA(s) ###

## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the number of rare alleles per locus 
rare_allele_number_aov_full <- aov(rare_allele_count ~ w_g + syst_err_rate + stoch_err_rate + syst_err_rate:stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(rare_allele_number_aov_full), main = "Residuals") 
#Interpretation:  pretty non normal I think bc it's not really a continuous indep variable?
dredge(rare_allele_number_aov_full, rank = "AICc") 
#Interpretation: first 4 models are all very similar, all include the dataset and stochastic error rate (2nd best includes only those 2 variables), will look at the best 2 and then decide if I want to see the other 2 best ones

rare_allele_number_top_model <- summary(lm(data = loci_info_df_for_graphing, rare_allele_count ~ w_g + syst_err_rate + stoch_err_rate))
rare_allele_number_top_model
# Interpretation: significant effect of wild and garden dataset types, significant effect of stochastic error rate and not systematic error rate

rare_allele_number_second_model <- summary(lm(data = loci_info_df_for_graphing, rare_allele_count ~ w_g + stoch_err_rate))
rare_allele_number_second_model
# Interpretation: all significant variables as the top model and nearly identical estimates

rare_allele_number_third_model <- summary(lm(data = loci_info_df_for_graphing, rare_allele_count ~ w_g + syst_err_rate * stoch_err_rate))
rare_allele_number_third_model
# Interpretation: all significant variables as the top model and nearly identical estimates, additionally stochastic error rate is nearly significant (p value of .0559)

# Summary: the garden datasets are significantly correlated to more (?) rare alleles per locus than the wild datasets (both are significantly correlated with higher rare alleles per locus but the wild has a smaller estimate), the stochastic error rate is significantly correlated to an increase in the number of rare alleles after accounting for variability due to dataset type, the systematic error rate is generally not a significant explanatory variable but is correlated to less rare alleles per locus, R2 around .14

#### Proportion of rare wild alleles in garden ####
# Plot of the proportion of wild rare alleles captured in garden samples in error vs real dataset across all loci,
loci_info_df_for_graphing %>%
  mutate(stoch_err_rate = as.factor(stoch_err_rate)) %>%
  ggplot() +
  geom_boxplot(aes(y = prop_rare_wild_alleles_captured, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate)), outlier.shape = NA) +
  # geom_jitter(aes(y = prop_rare_wild_alleles_captured, x = stoch_err_rate, color = interaction(syst_err_rate, stoch_err_rate)), alpha = .3) + # add points onto the boxplots
  ylab("Proportion of \"rare\" wild alleles in the garden sample") +
  xlab("stochastic error rate") +
  labs(colour = "Systematic error rate * stochastic error rate") +
  theme_classic()
# Interpretation: rare loci are less likely to be represented in the garden subsample in the error dataset


### ANOVA(s) ###

## Looking for the potential effect of locus, dataset, systematic error rate, and/or stochastic error rate on the proportion of rare wild alleles captured per locus
prop_rare_alleles_captured_aov_full <- aov(prop_rare_wild_alleles_captured ~ syst_err_rate + stoch_err_rate + syst_err_rate:stoch_err_rate, data = loci_info_df_for_graphing)
qqnorm.with.sim.bounds(resid(prop_rare_alleles_captured_aov_full), main = "Residuals") 
#Interpretation: wildly non-normal presumably because all values are between 0 and 1 and most seem to be 0 or 1 --> probably need to do a beta regression?
dredge(prop_rare_alleles_captured_aov_full, rank = "AICc") 
#Interpretation: top 4 models are all within 3 AICc points, note that best model is the null model and 2nd best model contains only systematic error rate

# Unsure how to code null model
prop_rare_alleles_captured_top_model <- summary(lm(data = loci_info_df_for_graphing, prop_rare_wild_alleles_captured))
prop_rare_alleles_captured_top_model
# Interpretation: significant effect of all included explanatory variables, note: since this metric is a ratio of garden to wild datasets, this is looking at the difference between real and error datasets

prop_rare_alleles_captured_second_model <- summary(lm(data = loci_info_df_for_graphing, prop_rare_wild_alleles_captured ~ syst_err_rate))
prop_rare_alleles_captured_second_model
# Interpretation: systematic error rate is not significant

prop_rare_alleles_captured_third_model <- summary(lm(data = loci_info_df_for_graphing, prop_rare_wild_alleles_captured ~ stoch_err_rate))
prop_rare_alleles_captured_third_model
# Interpretation: stochastic error rate is not significant

prop_rare_alleles_captured_fourth_model <- summary(lm(data = loci_info_df_for_graphing, prop_rare_wild_alleles_captured ~ stoch_err_rate + syst_err_rate))
prop_rare_alleles_captured_fourth_model
# Interpretation: neither stochastic error rate or systematic error rate is not significant

# Summary: the null model is clearly the best but will double check results with beta regression


# Since the anova residuals were so non-normal, trying a beta regression on a z transformed version of prop rare wild alleles

# Make the transformed column
trans_loci_info_df_for_graphing <- loci_info_df_for_graphing %>%
  filter(dataset %notin% c("garden_real", "garden_error"))  %>% # Remove columns with na so that the betregression won't fail
  filter(!is.nan(prop_rare_wild_alleles_captured))  %>%
  mutate(trans_prop_rare_wild_alleles_captured = (prop_rare_wild_alleles_captured*(nrow(.) - 1) + .5)/nrow(.))


betareg_prop_rare_wild_alleles_captured_full <- betareg(trans_prop_rare_wild_alleles_captured ~ syst_err_rate + stoch_err_rate + syst_err_rate:stoch_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
qqnorm.with.sim.bounds(resid(betareg_prop_rare_wild_alleles_captured_full), main = "Residuals")
#Interpretation: Still not normal at all...
dredge(betareg_prop_rare_wild_alleles_captured_full, rank = "AICc") 
#Interpretation: 4 top models within 3 AICc points, only model that is not includes all explanatory variables

# still not sure how to code null model but it is again the best model
betareg_prop_rare_wild_alleles_captured_top_model <- betareg(trans_prop_rare_wild_alleles_captured ~ , data = trans_loci_info_df_for_graphing, na.action = "na.fail")
summary(betareg_prop_rare_wild_alleles_captured_top_model)
# Interpretation: 

betareg_prop_rare_wild_alleles_captured_second <- betareg(trans_prop_rare_wild_alleles_captured ~ stoch_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
summary(betareg_prop_rare_wild_alleles_captured_second)
# Interpretation: stochastic error rate is not significant 

betareg_prop_rare_wild_alleles_captured_third <- betareg(trans_prop_rare_wild_alleles_captured ~ syst_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
summary(betareg_prop_rare_wild_alleles_captured_third)
# Interpretation: systematic error rate is not significant 

betareg_prop_rare_wild_alleles_captured_fourth <- betareg(trans_prop_rare_wild_alleles_captured ~ syst_err_rate + stoch_err_rate, data = trans_loci_info_df_for_graphing, na.action = "na.fail")
summary(betareg_prop_rare_wild_alleles_captured_fourth)
# Interpretation: neither systematic error rate nor stochastic error rate is significant 


# Summary: similar to the non beta-regression --> null model is clearly the best model



#### Allele frequencies ####

# Compare allele frequencies across all loci only in wild
allele_info_df_for_graphing %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = allele_freq, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Allele frequency") +
  xlab("Stochastic error rate") +
  theme_classic()



#### Resampling figures ####


# pull the first ind for each rep when 95% diversity is captured
ind_w_95_perc_diversity_captured <- resample_df_full %>%
  filter(prop_alleles_captured >= .95) %>%
  group_by(rep, locus, error_rep, syst_err_rate, stoch_err_rate) %>%
  slice(which.min(prop_alleles_captured)) %>%
  filter(locus == "all") %>%
  group_by(syst_err_rate, stoch_err_rate) %>%
  summarize(means = ceiling(mean(inds_sampled)), 
            sds = sd(inds_sampled)) %>%
  mutate(high_95 = ceiling(means + 2*sds), 
         low_95 = ceiling(means - 2*sds))



# make it so I can plot lines instead of so many damn points
lines_for_plot<- resample_df_full %>%
  group_by(inds_sampled, syst_err_rate, stoch_err_rate) %>%
  summarize(mean = mean(prop_alleles_captured), 
            sd = sd(prop_alleles_captured)) %>%
  mutate(high_95 = mean + 2*sd, 
         low_95 = mean - 2*sd) # These look completely wrong for reasons that are unclear to me and I feel I've wasted enough time on this for now that I'm going to move on and just not plot these


# Plot all of the proportion of alleles captured by error rates
resample_df_full %>%
  filter(locus == "all") %>%
  ggplot(group = interaction(syst_err_rate, stoch_err_rate)) +
  facet_wrap(~syst_err_rate * stoch_err_rate) +
  #geom_point(aes(x = inds_sampled, y = prop_alleles_captured), alpha = .025, color = "light gray") +
  geom_line(data = lines_for_plot, aes(x = inds_sampled, y = mean)) +
  geom_vline(data = ind_w_95_perc_diversity_captured, aes(xintercept = low_95, group = interaction(syst_err_rate, stoch_err_rate)), linetype = 2, color = "gray",) +
  geom_vline(data = ind_w_95_perc_diversity_captured, aes(xintercept = high_95, group = interaction(syst_err_rate, stoch_err_rate)), linetype = 2, color = "gray", ) +
  geom_vline(data = ind_w_95_perc_diversity_captured, aes(xintercept = means, group = interaction(syst_err_rate, stoch_err_rate)), color = "gray", ) +
  geom_text(data = ind_w_95_perc_diversity_captured, aes(x = low_95, label = low_95, y = 0.5, angle = 90, vjust = -0.2, group = interaction(syst_err_rate, stoch_err_rate)), size = 2.5) +
  geom_text(data = ind_w_95_perc_diversity_captured, aes(x = high_95, label = high_95, y = 0.5, angle = 90, vjust = -0.2, group = interaction(syst_err_rate, stoch_err_rate)), size = 2.5) +
  geom_text(data = ind_w_95_perc_diversity_captured, aes(x = means, label = means, y = 0.5, angle = 90, vjust = -0.2, group = interaction(syst_err_rate, stoch_err_rate)), size = 2.5) +
  geom_hline(yintercept = .95, alpha = .5, color = "red") +
  ylab("Prop of alleles captured") +
  xlab("Number of individuals sampled") +
  theme_classic() 


resample_df_full %>%
  filter(prop_alleles_captured >= .95) %>%
  group_by(rep, locus, error_rep, syst_err_rate, stoch_err_rate) %>%
  slice(which.min(prop_alleles_captured)) %>%
  filter(locus == "all") %>%
  #filter(stoch_err_rate %in% c(0, .01, .1)) %>%
  #filter(syst_err_rate %in% c(0, .01, .1)) %>%
  ggplot() +
  geom_boxplot(aes(y = inds_sampled, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  #geom_point(aes(y = inds_sampled, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate), shape = as.factor(error_rep)), alpha = .3) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("# of individuals needed to sample 95% allelic diversity") +
  xlab("Stochastic error rate") +
  theme_classic()
#ggsave(paste0(sim.wd, "/Figs/resample_95_perc_boxplots_smaller.png"), width = 9, height = 5)


#### Figs for Botany presentation #### 

loci_info_df_for_graphing  %>%
  filter(stoch_err_rate %in% c(0, .01, .1)) %>%
  filter(syst_err_rate %in% c(0, .01, .1)) %>%
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
ggsave(paste0(sim.wd, "/Figs/obs_het_boxplots_smaller.png"), width = 7, height = 3)



loci_info_df_for_graphing  %>%
  filter(stoch_err_rate %in% c(0, .01, .1)) %>%
  filter(syst_err_rate %in% c(0, .01, .1)) %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = num_alleles, x = as.factor(stoch_err_rate), color = as.factor(syst_err_rate))) +
  facet_grid(~ syst_err_rate) +
  labs(colour = "Systematic error rate") +
  ylab("Total number of alleles per locus") +
  xlab("Stochastic error rate") +
  theme_classic()
ggsave(paste0(sim.wd, "/Figs/num_alleles_boxplots_smaller.png"), width = 7, height = 3)


#### Less important graphs####

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




#### Want to make a graph that shows the number of reps that have changed the identity of the 1st 2nd and 3rd least frequent allele!####

