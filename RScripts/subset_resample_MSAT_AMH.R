# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE MSAT AND SNP SIMULATED DATASETS %%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# This script reads in genind files generated from previously run fastSimcoal simulations,
# then subsets each file to specify a group of "ex situ" (garden) individuals. The ex situ
# representation (i.e. how well do these garden individuals represent total allelic diversity) is calculated.

# Then, the remaining individuals ("wild") are resampled iteratively, and the allelic diversity
# of sample subsets (in comparison to the whole of wild allelic diversity) is calculated, then plotted

# In order to function iteratively over large objects (i.e. lists of genind objects), the steps
# in this script use many apply family functions

#library(strataG) #for SNPS
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

# My packages
library(ggplot2)
library(forcats)

# Set working directory to the GitHub repo (containing scripts and fastSimcoal outputs)
#sim.wd <- "/home/akoontz/Shared/SSRvSNP_Sim/Code/"
sim.wd <- "/Users/Ashley/Desktop/Hoban/Hoban_rotation_2023/Morton_SSRvSNP_Simulations/"
arq.wd <- "/Users/Ashley/Desktop/Hoban/Hoban_rotation_2023/Morton_SSRvSNP_Simulations/SimulationOutputs/MSAT_marker/MSAT_ArlequinFiles"
setwd(sim.wd)

# Parallelism: specify number of cores to use
num_cores <- detectCores() - 8
num_cores <- 4

####Creating functions for importing initial MSAT data####

# Create a function that is the opposite of %in%
`%notin%` <- Negate(`%in%`)

#A function that converts an arlequin simulation file to a genepop file
arp2gen<- function (infile){
  #browser()
  flForm <- strsplit(infile, split = "\\.")[[1]]
  if (substr(infile, 1, 2) == "./") {
    flForm <- flForm[-1]
  }
  else if (substr(infile, 1, 3) == "../") {
    flForm <- flForm[-(1:2)]
  }
  if (length(flForm) > 3) {
    stop("There were multiple '.' characters in your file name!")
  }
  tstfile <- paste(flForm[1], ".gen", sep = "")
  if (!file.exists(tstfile)) {
    fastScan <- function(fname) {
      s <- file.info(fname)$size
      buf <- readChar(fname, s, useBytes = TRUE)
      if (length(grep("\r", buf)) != 0L) {
        buf <- gsub("\r", "\n", buf)
        buf <- gsub("\n\n", "\n", buf)
      }
      return(strsplit(buf, "\n", fixed = TRUE, useBytes = TRUE)[[1]])
    }
    dat <- fastScan(infile)
    dat <- gsub("^\\s+|\\s+$", "", dat)
    dataType <- grep("*datatype=*", tolower(dat))
    if (strsplit(dat[dataType], "=")[[1]][2] != "MICROSAT") {
      stop("Data are not in 'MICROSAT' format!")
    }
    missDataLine <- grep("*missingdata=*", tolower(dat))
    missData <- noquote(substr(dat[missDataLine], nchar(dat[missDataLine]) -
                                 1, nchar(dat[missDataLine]) - 1))
    sampSizeLine <- grep("*samplesize=*", tolower(dat))
    # if (length(sampSizeLine) > 1) {
    sampNpos <- sapply(sampSizeLine, function(i) {
      return(regexpr("=", dat[i])[1])
    })
    # }
    popSizes <- as.numeric(substr(dat[sampSizeLine], start = sampNpos +
                                    1, stop = nchar(dat[sampSizeLine])))/2 #EDITED: dividing by 2 bc Austins popSize is half of the number of lines in the file?
    npops <- length(popSizes)
    sampStrt <- grep("*sampledata=*", tolower(dat))
    strts <- sapply(sampStrt, function(x) {
      if (dat[(x + 1)] == "") {
        return(x + 2)
      }
      else {
        return(x + 1)
      }
    })
    ends <- strts + ((popSizes * 2) - 1)
    nloci <- length(strsplit(dat[strts[1]], split = "\\s+")[[1]]) - 2
    popGeno <- lapply(seq_along(strts), function(i) {
      return(dat[strts[i]:ends[i]])
    })
    popSzcheck <- sapply(popGeno, function(x) length(x))/2 #EDITED: divide by 2 here bc popgeno needs to be double the size of the pops bc these are diploid 
    if (!all(identical(popSzcheck, popSizes))) {
      stop("Failed! Please make sure that your file is formatted correctly.")
    }
    popIdx <- lapply(popGeno, function(x) {
      return(seq(1, length(x), 2)) #EDITED: turned this seq to add by 2 across the length of the genotypes in the pop, makes this count only hold every other line, enables the concatenation of every 2 lines to make a single ind (as they are formatted by Austin's output)
    })
    #somewhere in here an error is getting thrown
    popList <- lapply(seq_along(popGeno), function(i) {
      al1 <- matrix(unlist(strsplit(popGeno[[i]][popIdx[[i]]],
                                    split = "\\s+")), nrow = popSizes[i], ncol = nloci + 2, byrow = TRUE)[, -(1:2)] #EDITED: had to edit here to tell it how many columns this matrix should be... not sure why it didn't want to infer that for Austin's files but it didn't
      al2 <- matrix(unlist(strsplit(popGeno[[i]][(popIdx[[i]] +
                                                    1)], split = "\\s+")), nrow = popSizes[i], ncol = nloci + 2, byrow = TRUE)[, -(1:2)] #EDITED: had to edit here to tell it how many columns this matrix should be and delete the first two columns bc each line has that extra individual info for Austin's files
      tst <- matrix(paste(al1, al2, sep = ""), nrow = popSizes[i])
      tst <- cbind(paste(rep("pop", nrow(tst)), i, " ,",
                         sep = ""), tst)
      rm(al1, al2)
      z <- gc()
      rm(z)
      if (nchar(tst[1, 2]) == 4) {
        tst[tst == paste(missData, missData, sep = "")] <- "0000"
      }
      else {
        tst[tst == paste(missData, missData, sep = "")] <- "000000"
      }
      out <- apply(tst, 1, function(x) {
        return(paste(x, collapse = "\t"))
      })
      out <- c("POP", out)
      rm(tst)
      z <- gc()
      rm(z)
      return(out)
    })
    outfile <- strsplit(infile, "\\.")[[1]]
    if (length(outfile) >= 2) {
      outfile <- paste(outfile[-length(outfile)], collapse = ".")
    }
    else {
      outfile <- outfile[1]
    }
    loci <- paste("locus", 1:nloci, sep = "")
    loci <- c(paste(outfile, "_gen_converted", sep = ""),
              loci)
    of <- c(loci, unlist(popList))
    out <- file(paste(outfile, ".gen", sep = ""), "w")
    for (i in 1:length(of)) {
      cat(of[i], "\n", file = out, sep = "")
    }
    close(out)
    return(TRUE)
  }
  else {
    return(NULL)
  }
}

#A function that converts all arlequin simulation files in a folder to genepop files
import_arp2gen_files = function(mypath, mypattern) {
  setwd(mypath) #directory containing files
  temp_list_1 = list.files(mypath, mypattern) #list of files with .arp extension
  temp_list_2 = list(length = length(temp_list_1)) #empty list with the same length as the first list
  for(z in 1:length(temp_list_1)){temp_list_2[[z]]=arp2gen(temp_list_1[z])} # for every item in list 1, convert it to .gen file and save to list 2
  temp_list_2
} 

#A function that converts genepop files to genalex files in 2 main steps:
#   first you need to convert genepop files to genind objects
#   then you convert genind objects to genalex files
import_gen2genalex_files = function(mypath, mypattern) {
  setwd(mypath) #directory containing files
  temp_list_1 = list.files(mypath, mypattern) #list of all files with .gen extension
  temp_list_2 = list(length = length(temp_list_1)) #empty list with same length as the first
  temp_list_3 = list(length = length(temp_list_1))#another empty list with same length as the first
  for(z in 1:length(temp_list_1)){
    temp_list_2[[z]]=read.genepop(temp_list_1[z], ncode = 3)#converting each genepop file to a genind object
    #getting the name and filepath of each file so that they are preserved when converted to csv
    flForm <- strsplit(temp_list_1[z], split = "\\.")[[1]]
    if (substr(temp_list_1[z], 1, 2) == "./") {
      flForm <- flForm[-1]
    }
    else if (substr(temp_list_1[z], 1, 3) == "../") {
      flForm <- flForm[-(1:2)]
    }
    if (length(flForm) > 3) {
      stop("There were multiple '.' characters in your file name!")
    }
    temp_list_3[[z]]=genind2genalex(temp_list_2[[z]], filename = (paste(flForm[[1]],".csv", sep="")), overwrite = TRUE) #converting each genind object to a genalex file
    #will overwrite this so that Austin's pop names are preserved
  }
  temp_list_3 #retun
}

###functions for adding error and obtained allele frequencies

# Function for obtaining allele frequencies at each locus
# Input: Desired dataframe for which you will obtain allele frequency data 
# Output: A list of lists, each entry in list = info for that locus number
#                   Element 1) vector of allele frequencies at that locus
#                   Element 2) matrix with 2 rows, row 1= allele, row 2 = frequency of alleles
get_alleles_and_freqs <- function(input_data){
  
  # Pulling just the loci from the df
  locus_data <- input_data[,grepl( "locus" , colnames(input_data))]
  loci <- seq(1, num_loci)
  alleles_and_freqs_list <- lapply(loci, function(locus){
    # Pulling all alleles at a single locus from all parents in the population
    all_alleles <- sort(unique(c(locus_data[, c(2* locus -1) ], locus_data[, (2* locus)])))
    
    # Getting the frequency of each allele by counting the number of occurrences of the allele at each of the 2 loci and dividing by the number of individuals in the input dataframe
    allele_freqs <- sapply(all_alleles, simplify = TRUE, USE.NAMES = TRUE, function(x){
      return((sum(locus_data[,c(2*locus - 1)] %in% x) + sum(locus_data[,c(2*locus)] %in% x)) / (nrow(locus_data)*2))})
  
    return(list(all_alleles,allele_freqs))})
  
  return(alleles_and_freqs_list)
}

# Fixed so more realistic type 2 errors = slippage or mis type of single msat_length
simulating_error <- function(data_for_editing, type1_rate, type2_rate, msat_length){
  #browser()
  
  # Pulling just the loci from the df
  locus_data <- data_for_editing[,grepl( "locus" , colnames(data_for_editing))]
  
  # Getting loci number from data 
  num_loci <- ncol(locus_data)/2
  
  # Making the new matrix with the same shape as the old matrix
  locus_data_new <- data.frame(matrix(ncol = ncol(locus_data), nrow = nrow(locus_data)))
  diffs_in_genos <- data.frame(matrix(ncol = ncol(locus_data)/2, nrow = nrow(locus_data)))
  
  for(locus in 1:num_loci){
    alleles <- get_alleles_and_freqs(input_data = data_for_editing)[[locus]][[1]] # Want to get all possible alleles from the parental data
    k = length(alleles) # k = number of alleles at designated locus
    e_1 = type1_rate / (1 + type1_rate)
    #e_2 = type2_rate / (k - 1) # Old error rate, when type 2 error = change to any other already known geno
    e_2 = type2_rate / 2 # Dividing by 2 bc in the new version of type 2 error the real geno can become either of 2 genos
    genos <- locus_data[, c(2*locus - 1, 2*locus)] # Keep only the genotypes at the designated locus
    
    for(ind in 1:nrow(genos)){
      
      old_geno <- paste0(genos[ind,1], "-", genos[ind,2])
      current_geno <- old_geno
      #Ask if genotype is het, if yes, then add chance of allelic dropout error
      if(genos[ind,1] != genos[ind,2]){
        het <- paste0(genos[ind,1],"-", genos[ind,2])
        homo1 <- paste0(genos[ind,1],"-", genos[ind,1])
        homo2 <- paste0(genos[ind,2],"-", genos[ind,2])
        possible_genos =  c(het, homo1, homo2)
        new_geno <- sample(possible_genos, size = 1, prob = c(1 - (2*e_1), e_1, e_1))
        current_geno <- new_geno
      }
      # Getting the unique alleles into separate vectors again
      allele1 <- as.numeric(unlist(strsplit(current_geno, "-"))[1])
      allele2 <- as.numeric(unlist(strsplit(current_geno, "-"))[2])
      
      # Adding the error that isn't due to alleleic dropout, this class of error typically occurs at lower rate (effects less inds) but will turn the allele to any other allele at an equal rate --> assume these happen after class 1 errors
      
      # Old version of adding type 2 error
      #new_allele1 <- sample(c(allele1, alleles[alleles %notin% allele1]), size = 1, prob = c(1-e_2, rep(e_2/(k-1), k - 1)))
      #new_allele2 <- sample(c(allele2, alleles[alleles %notin% allele2]), size = 1, prob = c(1-e_2, rep(e_2/(k-1), k - 1)))
      
      # New version of adding type 2 error
      new_allele1 <- sample(c(allele1, allele1 + msat_length, allele1 - msat_length),  size = 1, prob = c(1- (2 * e_2), rep(e_2, 2)))
      new_allele2 <- sample(c(allele2, allele2 + msat_length, allele2 - msat_length),  size = 1, prob = c(1- (2 * e_2), rep(e_2, 2)))
      
      # Putting the new genos into a data set
      locus_data_new[ind, 2*locus - 1] <- as.character(new_allele1)
      locus_data_new[ind, 2*locus] <- as.character(new_allele2)
      
      # Checking if the new genos are *functionally* diff from the old genos
      new_geno_v1 <- paste0(new_allele1, "-", new_allele2)
      new_geno_v2 <- paste0(new_allele1, "-", new_allele2)
      
      # Noting whether or not the genotype has changed and putting answer in the diffs_in_genos matrix
      # Want to change this so it notes
      geno_changed <- "F"
      if (old_geno %notin% c(new_geno_v1, new_geno_v2)) {
        geno_changed <- "T"
      }
      diffs_in_genos[ind, locus] <- geno_changed
      
    }
  }
  diffs_in_genos<- as.data.frame(lapply(diffs_in_genos, as.logical)) # Making this output as actual TRUE/FALSE values that R recognizes so I can easily count the number of occurrences
  
  # Fixing column names of both output dfs so they are compatible with my current naming conventions 
  colnames(locus_data_new) <- colnames(locus_data)
  colnames(diffs_in_genos) <- unique(str_sub(colnames(locus_data_new), start = 1, end = -2))
  
  return(list(locus_data_new, diffs_in_genos))     
}


# Function to get the smallest alleles and their proportions from the current nested list format
get_df_of_sm_alleles <- function(x) {   
  prop_df <- x[[2]]
  # Get smallest allele info
  smallest_prop <- min(prop_df)
  smallest_allele <- which.min(prop_df)
  smallest_allele_num <- x[[1]][smallest_allele]
  
  # Get 2nd smallest allele info
  sec_smallest_prop <- min(prop_df[prop_df!= smallest_prop])
  sec_smallest_allele <- which.min(prop_df[prop_df!= smallest_prop])
  sec_smallest_allele_num <- x[[1]][sec_smallest_allele]

  # Get 3rd smallest allele info
  third_smallest_prop <- min(prop_df[prop_df!= smallest_prop & prop_df!= sec_smallest_prop])
  third_smallest_allele <- which.min(prop_df[prop_df!= smallest_prop & prop_df!= sec_smallest_prop])
  third_smallest_allele_num <- x[[1]][third_smallest_allele]
  
  # Put all in one df
  df <- data.frame(smallest_allele_num, smallest_prop, sec_smallest_allele_num, sec_smallest_prop, third_smallest_allele_num, third_smallest_prop)
  return(df)
}

# Very similar to get_alleles_and_freqs function except that the input df should have a location column specifying whether the individual is in the garden or the wild and then outputs a list of lists containing relevant summary information at each locus
get_alleles_and_freqs_w_g <- function(input_data){
  
  # Pulling just the loci and location data from the df
  locus_data <- input_data[,grepl( "locus" , colnames(input_data)) | grepl( "location" , colnames(input_data))]
  loci <- seq(1, num_loci)
  alleles_and_freqs_list <- lapply(loci, function(locus){
  # Separate out wild and garden inds
   wild_locus_data <- subset(locus_data, locus_data$location == "wild")
   garden_locus_data <- subset(locus_data, locus_data$location == "garden")
   
   # Get alleles and freqs from the wild inds
    # Pulling all alleles at a single locus from all parents in the population
    all_wild_alleles <- sort(unique(c(wild_locus_data[, c(2* locus -1) ], wild_locus_data[, (2* locus)])))
    
    # Getting the frequency of each allele by counting the number of occurrences of the allele at each of the 2 loci and dividing by the number of individuals in the input dataframe
    allele_wild_freqs <- sapply(all_wild_alleles, simplify = TRUE, USE.NAMES = TRUE, function(x){
      return((sum(wild_locus_data[,c(2*locus - 1)] %in% x) + sum(wild_locus_data[,c(2*locus)] %in% x)) / (nrow(wild_locus_data)*2))})
    
    # Getting obs # of hets and divide by total # of inds to get obs het 
    obs_wild_het <- (nrow(wild_locus_data) - sum(wild_locus_data[,c(2*locus - 1)] == wild_locus_data[,c(2*locus)]))/nrow(wild_locus_data)
    
    # Get alleles and freqs from the garden inds
    # Pulling all alleles at a single locus from all parents in the population
    all_garden_alleles <- sort(unique(c(garden_locus_data[, c(2* locus -1) ], garden_locus_data[, (2* locus)])))
    
    # Getting the frequency of each allele by counting the number of occurrences of the allele at each of the 2 loci and dividing by the number of individuals in the input dataframe
    allele_garden_freqs <- sapply(all_garden_alleles, simplify = TRUE, USE.NAMES = TRUE, function(x){
      return((sum(garden_locus_data[,c(2*locus - 1)] %in% x) + sum(garden_locus_data[,c(2*locus)] %in% x)) / (nrow(garden_locus_data)*2))})
    
    # Getting obs # of hets and divide by total # of inds to get obs het 
    obs_garden_het <- (nrow(garden_locus_data) - sum(garden_locus_data[,c(2*locus - 1)] == garden_locus_data[,c(2*locus)]))/nrow(garden_locus_data) 
    
    
    return(list(all_wild_alleles,allele_wild_freqs, obs_wild_het, all_garden_alleles, allele_garden_freqs, obs_garden_het))})
  
  return(alleles_and_freqs_list)
}


# A function that takes in a locus number and a dataframe that contains garden and wild alleles (output of get_alleles_and_freqs_w_g)
compare_w_g<- function(locus, input_real_data, input_error_data){
  
  # Pull out the list (which has 6 subset lists) that corresponds to the input locus
  locus_real_data <- input_real_data[[locus]]
  wild_real_alleles <- locus_real_data[[1]] # First item in locus_data list is a list of unique alleles in the wild 
  garden_real_alleles <-locus_real_data[[4]] # Fourth item in locus_data list is a list of unique alleles in the garden 
  
  locus_error_data <- input_error_data[[locus]]
  wild_error_alleles <- locus_error_data[[1]] # First item in locus_data list is a list of unique alleles in the wild 
  garden_error_alleles <-locus_error_data[[4]] # Fourth item in locus_data list is a list of unique alleles in the garden 
  
  # this will need to include all alleles across real and fake runs
  all_alleles <- unique(c(wild_real_alleles, garden_real_alleles, wild_error_alleles, garden_error_alleles)) # Make a vector of all alleles
  
  # Df with allele specific info
  alleles_df <- data.frame(locus = rep(locus, length(all_alleles))) %>%
    mutate(all_alleles = all_alleles,
           in_wild_real = all_alleles %in% wild_real_alleles, # If allele is in wild in real dataset, returns TRUE
           in_garden_real = all_alleles %in% garden_real_alleles, # If allele is in garden in real dataset, returns TRUE
           in_wild_error = all_alleles %in% wild_error_alleles, # If allele is in wild in real dataset, returns TRUE
           in_garden_error = all_alleles %in% garden_error_alleles, # If allele is in garden in real dataset, returns TRUE
           wild_real = ifelse(all_alleles %in% wild_real_alleles, locus_real_data[[2]][all_alleles], 0), # Second item is a list of unique alleles in the wild 
           garden_real = ifelse(all_alleles %in% garden_real_alleles, locus_real_data[[5]][all_alleles], 0), # Returns freq of allele in garden, if not present in garden then 0
           wild_error = ifelse(all_alleles %in% wild_error_alleles, locus_error_data[[2]][all_alleles], 0), # Second item is a list of unique alleles in the wild 
           garden_error = ifelse(all_alleles %in% garden_error_alleles, locus_error_data[[5]][all_alleles], 0)) %>% # Returns freq of allele in garden, if not present in garden then 0
    pivot_longer(cols = c(wild_real, garden_real, wild_error, garden_error), names_to = "dataset", values_to = "allele_freq") %>%
    # These bins apparently come from Hoban et al 2020
    mutate(allele_freq_cat = ifelse(allele_freq == 0, "Not present", 
                                    ifelse(allele_freq <= .01, "Rare", 
                                           ifelse(allele_freq > .01 & allele_freq <= .1, "Low Frequency", "Very Common")))) %>%
    mutate(common_allele = ifelse(allele_freq > .05, T, F))
    
  
  # Df with info across allleles summarized
  loci_df <- data.frame(locus) %>%
    # Makec columns with observed heterozygosity
    mutate(wild_real = locus_real_data[[3]],  # Third item in locus_data list is a the number observed heterozygotes in the wild
           garden_real = locus_real_data[[6]], # Sixth item in locus_data list is a the number observed heterozygotes in the garden
           wild_error = locus_error_data[[3]],  # Third item in locus_data list is a the number observed heterozygotes in the wild
           garden_error = locus_error_data[[6]]) %>% # Sixth item in locus_data list is a the number observed heterozygotes in the garden 
    # Make above tidy
    pivot_longer(cols = c(wild_real, garden_real, wild_error, garden_error), names_to = "dataset", values_to = "exp_het") %>%
    
    # Calc expected heterozygosity with 1 - (exp frequencies of all homozygotes based on frequencies of each allele)
    mutate(obs_het = ifelse(dataset == "wild_real", 1-sum(na.omit(subset(alleles_df, dataset == "wild_real", select = c(allele_freq)))^2),
                            ifelse(dataset == "garden_real", 1-sum(na.omit(subset(alleles_df, dataset == "garden_real", select = c(allele_freq)))^2),
                                   ifelse(dataset == "wild_error", 1-sum(na.omit(subset(alleles_df, dataset == "wild_error", select = c(allele_freq)))^2), 1-sum(na.omit(subset(alleles_df, dataset == "garden_error", select = c(allele_freq)))^2))))) %>%
    
    # Add column with number of alleles per dataset
    mutate(num_alleles = ifelse(dataset == "wild_real", length(wild_real_alleles), 
                         ifelse(dataset == "garden_real", length(garden_real_alleles),  
                         ifelse(dataset == "wild_error", length(wild_error_alleles), length(garden_error_alleles)))),
           prop_rare_wild_alleles_captured = ifelse(dataset == "wild_real", sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Rare" & alleles_df$in_garden_real == T)/sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Rare"),
                                             ifelse(dataset == "wild_error", sum(grepl("wild_error", alleles_df$dataset) & alleles_df$allele_freq_cat == "Rare" & alleles_df$in_garden_error == T)/sum(grepl("wild_error", alleles_df$dataset) & alleles_df$allele_freq_cat == "Rare"), NA)), 
           prop_lowfreq_wild_alleles_captured = ifelse(dataset == "wild_real", sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Low Frequency" & alleles_df$in_garden_real == T)/sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Low Frequency"),
                                             ifelse(dataset == "wild_error", sum(grepl("wild_error", alleles_df$dataset) & alleles_df$allele_freq_cat == "Low Frequency" & alleles_df$in_garden_error == T)/sum(grepl("wild_error", alleles_df$dataset) & alleles_df$allele_freq_cat == "Low Frequency"), NA)), 
           prop_vcommon_wild_alleles_captured = ifelse(dataset == "wild_real", sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Very Common" & alleles_df$in_garden_real == T)/sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Very Common"),
                                              ifelse(dataset == "wild_error", sum(grepl("wild_error", alleles_df$dataset) & alleles_df$allele_freq_cat == "Very Common" & alleles_df$in_garden_error == T)/sum(grepl("wild_error", alleles_df$dataset) & alleles_df$allele_freq_cat == "Very Common"), NA)), 
           rare_allele_count = ifelse(dataset == "wild_real", sum(alleles_df$allele_freq_cat == "Rare" & alleles_df$dataset == "wild_real"),
                               ifelse(dataset == "garden_real", sum(alleles_df$allele_freq_cat == "Rare" & alleles_df$dataset == "garden_real"), 
                               ifelse(dataset == "wild_error", sum(alleles_df$allele_freq_cat == "Rare" & alleles_df$dataset == "wild_error"), 
                                      sum(alleles_df$allele_freq_cat == "Rare" & alleles_df$dataset == "garden_error")))))
  
  return(list(alleles_df, loci_df))
  
}


####Loading the already simulated population data####
num_loci = 20 #number of loci simulated, needed to make a data frame to save the data, starts with 20 --> gets cut down to 10 later

# Make a vector of all of the loci names
loci_names = c()
for(i in 1:num_loci){
  loci_names = c(loci_names, paste("locus", i, "a", sep=""))
  loci_names = c(loci_names, paste("locus", i, "b", sep=""))
}

#importing and converting arlequin files to genepop files
import_arp2gen_files(arq.wd,".arp$")


# this command throws any error claiming "ncode = 3 isn't true (some alleles don't have 3 characters)" --> works for the first 4 files but then doesn't??
#importing and converting genepop files to genalex
import_gen2genalex_files(arq.wd, ".gen$")

#list of genalex files for all simulation replicates--genalex files end in .csv
genalex_list = list.files(arq.wd, ".csv$")

#read in and then cut off first 2 rows in each data frame -- this is the population data, which is not required for our purposes
i = genalex_list[[1]]
final_output <- lapply(genalex_list, function(i){
  genetic_data <- read.csv(paste(arq.wd, "/", i, sep=""), header=FALSE)[-c(1,2),]
  #giving the data frame columns new names
  names(genetic_data) = c("Ind", "Pop", loci_names)
  
#this genetic_data object now contains all of the data of 20 MSAT loci from each of individuals (in a single simulation) as simulated from SIM2COAL
  
  genetic_data = as.data.frame(genetic_data[-1,]) #removing the first row and converting to a dataframe
  }
)


#### All things below are done to a single entry in the genetic_data list####


# Adding error to the data
genetic_data_error <- simulating_error(data_for_editing = genetic_data, type1_rate = .1, type2_rate = .05, msat_length = 1)

genetic_data_sim_error <- genetic_data_error[[1]] # Error dataset to use

genetic_data_track_error <- genetic_data_error[[2]] # For error summarization purposes
colSums(genetic_data_track_error)/(nrow(genetic_data_sim_error)*2) # Calc total error rate at each locus



# Getting the frequencies of each allele from each locus from the "real" dataset
allele_freqs_real <- get_alleles_and_freqs(input_data = genetic_data)

# Get the number of alleles for each locus into a data frame
allele_counts_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_real, function(x) lengths(x)[1])))



# Getting the frequency of each allele from each locus from the dataset with error added to it
allele_freqs_error <- get_alleles_and_freqs(input_data = genetic_data_sim_error)

# Get the number of alleles for each locus into a data frame
allele_counts_error_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_error, function(x) lengths(x)[1])))

sum(allele_counts_error_df == allele_counts_df)/num_loci #Gets proportion of loci where alleles were lost/gained by addition of errors


# Get the proportion of the 3 least common alleles and identity of those alleles for all loci using get_df_of_sm_alleles
rarest_allele_prop_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_real, get_df_of_sm_alleles)))

rarest_allele_prop_df_error <- as.data.frame(do.call(rbind, lapply(allele_freqs_error, get_df_of_sm_alleles)))

# Merge the real and error dfs for visualization purposes
comparsion_df <- rarest_allele_prop_df %>% 
  rename(first_allele_real= smallest_allele_num, 
         first_prop_real = smallest_prop,
         sec_allele_real= sec_smallest_allele_num, 
         sec_prop_real = sec_smallest_prop,
         third_allele_real= third_smallest_allele_num, 
         third_prop_real = third_smallest_prop) %>%
  cbind(rarest_allele_prop_df_error) %>%
  rename(first_allele_error= smallest_allele_num, 
         first_prop_error = smallest_prop,
         sec_allele_error= sec_smallest_allele_num, 
         sec_prop_error = sec_smallest_prop,
         third_allele_error= third_smallest_allele_num, 
         third_prop_error = third_smallest_prop) %>%
  mutate(locus = row_number(),
         first_higher_prop = ifelse(0 <= first_prop_error - first_prop_real, "error", "real"), 
         sec_higher_prop = ifelse(0 <= sec_prop_error - sec_prop_real, "error", "real"), 
         third_higher_prop = ifelse(0 <= third_prop_error - third_prop_real, "error", "real")) 


# Make the graph of the freq of the 3 rarest alleles w/ and without error
comparsion_df %>%
  ggplot(aes(x = locus, y = first_prop_real)) +
  geom_point(color = "black") +
  geom_point(aes(y = first_prop_error), color = "red") + 
  geom_point(aes(y = sec_prop_real), color = "black") +
  geom_point(aes(y = sec_prop_error), color = "red") + 
  geom_point(aes(y = third_prop_real), color = "black") +
  geom_point(aes(y = third_prop_error), color = "red") + 
  geom_segment(aes(x = locus, xend = locus, y= first_prop_real, yend = first_prop_error, color = first_higher_prop)) +
  geom_segment(aes(x = locus, xend = locus, y= sec_prop_real, yend = sec_prop_error, color = sec_higher_prop)) + 
  geom_segment(aes(x = locus, xend = locus, y= third_prop_real, yend = third_prop_error, color = third_higher_prop)) +
  scale_color_manual(values = c("red", "black")) +
  ylab("Frequency of alleles") +
  labs(title = "Frequency of 3 rarest alleles in the population with and without error") +
  theme_classic()



#Now subset the initial dataset to get wild + garden dfs?
perc_sampled <- .2
pops <- c(unique(genetic_data$Pop), "pop2")
genetic_data$location <- "wild"
garden_inds <- unlist(lapply(pops, function(x) {
  inds <- genetic_data$Ind[genetic_data$Pop == x]
  garden <- sample(inds, size = perc_sampled * length(inds), replace = F)
  return(garden)
}))

genetic_data[genetic_data$Ind %in% garden_inds, c("location")] <- "garden"

# Making the error dataset have all the info present in the real dataset 
genetic_data_error <- cbind(genetic_data$Ind, genetic_data$Pop, genetic_data_sim_error, genetic_data$location)
colnames(genetic_data_error) <- colnames(genetic_data)


# Now compare the wild data and the garden data (allele count, allele frequency in wild, heterozygosity at each allele in population)


# Get the lists that have all of the info I need for graphing
w_g_real_alleles_and_freqs <- get_alleles_and_freqs_w_g(genetic_data)
w_g_error_alleles_and_freqs <- get_alleles_and_freqs_w_g(genetic_data_error)

# Make the full dataframes that holds all summary info across all alleles/loci
loci <- seq(1, num_loci)
full_dfs_in_lists <- lapply(loci, compare_w_g, w_g_real_alleles_and_freqs, w_g_error_alleles_and_freqs)
allele_info_df <- do.call(rbind, lapply(full_dfs_in_lists, function(i){
  return(i[[1]])}))
loci_info_df <- do.call(rbind, lapply(full_dfs_in_lists, function(i){
  return(i[[2]])}))


#### Making some graphs####
allele_info_df %>%
  ggplot(alpha = .7) +
  geom_point(aes(x = all_alleles, y = allele_freq, color = dataset), shape = 1) +
  scale_color_manual(values = c("indianred", "gray",  "red", "black")) +
  facet_wrap(~ locus) +
  theme_classic()


allele_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot() +
  geom_boxplot(aes(y = allele_freq, x = dataset, color = dataset), outlier.shape = NA) +
  geom_point(aes(y = allele_freq, x = dataset, color = dataset), alpha = .3) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  #facet_wrap(~ locus) +
  theme_classic()

# Note an overly important graph, just shows differences between expected and obs het for each data set split by each locus
loci_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot() +
  geom_point(aes(y = exp_het, x = dataset, color = dataset)) +
  geom_point(aes(y = obs_het, x = dataset, color = dataset)) +
  geom_segment(aes(y = obs_het, yend = exp_het, x = dataset, xend= dataset), color = "black") +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_wrap(~ locus) +
  theme_classic()


# Comparison of observed heterozygosity of all loci in each dataset 
loci_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = dataset, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  #facet_wrap(~ locus) +
  theme_classic()


# Comparison of expected heterozygosity of all loci in each dataset (error is lower bc more alleles at a lower frequency)
loci_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot() +
  geom_boxplot(aes(y = exp_het, x = dataset, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  #facet_wrap(~ locus) +
  theme_classic()


# Comparison of number of across all loci in each dataset (indeed we see more alleles in the error dataset)
loci_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot() +
  geom_boxplot(aes(y = num_alleles, x = dataset, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  #facet_wrap(~ locus) +
  theme_classic()


# Across all loci, proportion of wild rare alleles captured in garden samples in error vs real dataset (rare loci are less likely to be represented in the garden sample in the error dataset)
loci_info_df %>%
  filter(dataset %in% c("wild_error", "wild_real")) %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot() +
  geom_boxplot(aes(y = prop_rare_wild_alleles_captured, x = dataset, color = dataset), outlier.shape = NA) +
  geom_jitter(aes(y = prop_rare_wild_alleles_captured, x = dataset, color = dataset), alpha = .3) +
  scale_color_manual(values = c("black", "red")) +
  #facet_wrap(~ locus) +
  theme_classic()


loci_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot(aes(x = rare_allele_count)) +
  geom_bar(aes(fill = dataset), position=position_dodge2(preserve = "single")) +
  scale_fill_manual(values = c("black","gray", "red",  "indianred")) +
  ylab("Number of loci with given amount of rare alleles") +
  xlab("Number of rare alleles") +
  theme_classic()
  #facet_wrap(~ dataset)


loci_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error")) %>%
  ggplot(aes(y = rare_allele_count)) +
  geom_boxplot(aes(x = dataset, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  ylab("Number of rare alleles/locus") +
  xlab("dataset") +
  theme_classic()
#facet_wrap(~ dataset)


#by locus: get unique alleles from wild pop, use an %in% to ask if wild alleles are in garden samples (maybe count the number of yesses) --> generate summary table per locus that will give if the wild allele is in the garden, how many times it's in the garden and frequency of that allele in the wild population 

#can have a three column table (allele name, count in garden), put that in an 3D array such that each plan of the array reps a locus

#should have some alleles that exist in the wild and not in the garden

# do above in errored and non-errored dataset and compare via boxplots


