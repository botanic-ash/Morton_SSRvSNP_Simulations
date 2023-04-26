# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%% SUBSET AND RESAMPLE MSAT SIMULATED DATASETS %%%
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
library(ggnewscale)

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
simulating_error <- function(data_for_editing, type1_rate, type2_rate, msat_length, replicate_number){
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
  
  # Adding replicate number to the output dfs
  locus_data_new <-  locus_data_new %>%
    mutate(rep = replicate_number,
           type1_rate = type1_rate, 
           type2_rate = type2_rate)
  diffs_in_genos <- diffs_in_genos %>%
    mutate(rep = replicate_number,
           type1_rate = type1_rate,
           type2_rate = type2_rate)
  
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

loci <- seq(1, num_loci)

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

# Make current genetic data 
genetic_data <- read.csv(paste(arq.wd, "/", genalex_list[1], sep=""), header=FALSE)[-c(1,2),]
#giving the data frame columns new names
names(genetic_data) = c("Ind", "Pop", loci_names)
genetic_data = as.data.frame(genetic_data[-1,]) #removing the first row and converting to a dataframe


#### All things below are done to a single entry in the genetic_data list####

#plan is to make the simulating error code output a list of dataframes so that everything else can be wrapped into an lapply call and then at the end I can rbind my data together, will just need to make sure I have a column which specifies which replicate it is and can then just group_by replicate and get medians   

###Adding error to the data
type1_rate = c(.01, .05, .1)
type2_rate = c(.01, .05, .1)

rate_list <- expand.grid(type1_rate, type2_rate) %>%
  rename(type1_rate = Var1, 
         type2_rate = Var2)
rate_combo_count <- seq(1, nrow(rate_list), 1)

# This is pretty clunky and slow but it get's the job done, could be parralelized fairly easily though I think
# output is a list of lists of lists, 1st list set reps all diff rate value combos, 2nd list reps each replicate error dataset with those error rate values, 3rd list contains the error dataset and the tracked error dataset
genetic_data_error <- lapply(rate_combo_count, function(x)
  lapply(reps, function(i)
       simulating_error(data_for_editing = genetic_data, type1_rate = rate_list$type1_rate[x], type2_rate = rate_list$type2_rate[x], msat_length = 1, replicate_number = i)))



genetic_data_track_error <- do.call(rbind, lapply(genetic_data_error, function(x)
  do.call(rbind, lapply(x, function(i) 
  colSums(i[[2]])/(nrow(i[[2]])*2)))))  # Calc total error rate at each locus

as.data.frame(genetic_data_track_error) %>%
  mutate(rep = as.factor(rep *2), 
         type1_rate = as.factor(type1_rate *2), 
         type2_rate = as.factor(type2_rate *2)) %>%
  pivot_longer(cols = starts_with("locus"), names_to = "locus", values_to = "error_rate") %>%
  ggplot() +
  geom_hline(aes(yintercept = median(error_rate)), alpha = .3) +
  geom_boxplot(aes(x = locus, y = error_rate), outlier.shape = NA) + 
  #geom_point(aes(x = locus, y = error_rate, color = rep)) +
  facet_wrap(~ type1_rate * type2_rate) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

as.data.frame(genetic_data_track_error) %>%
  mutate(rep = as.factor(rep *2), 
         type1_rate = as.factor(type1_rate *2), 
         type2_rate = as.factor(type2_rate *2)) %>%
  pivot_longer(cols = starts_with("locus"), names_to = "locus", values_to = "error_rate") %>%
  ggplot() +
  geom_hline(aes(yintercept = median(error_rate)), alpha = .3) +
  geom_boxplot(aes(y = error_rate, group = interaction(type1_rate, type2_rate), color = interaction(type1_rate, type2_rate))) +
  ylab("real error rate") +
  theme_classic() +
  #facet_wrap(~ type1_rate * type2_rate) +
  theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())



# Getting the frequencies of each allele from each locus from the "real" dataset
allele_freqs_real <- get_alleles_and_freqs(input_data = genetic_data)

# Get the number of alleles for each locus into a data frame
allele_counts_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_real, function(x) lengths(x)[1])))


# Pull out just the errored reads data frames
genetic_data_sim_error <- lapply(rate_combo_count, function(x)
  lapply(reps, function(i) genetic_data_error[[x]][[i]][[1]]))

# Getting the frequency of each allele from each locus from the dataset with error added to it
allele_freqs_error <- lapply(genetic_data_sim_error, function(x)
  lapply(x, get_alleles_and_freqs))


# Get the number of alleles for each locus for each rep into a data frame

allele_counts_error_df <- as.data.frame(do.call(rbind, lapply(rate_combo_count, function(x)
  as.data.frame(do.call(rbind, lapply(reps, function(i)
    as.data.frame(do.call(rbind, lapply(loci, function(z)
    return(lengths(allele_freqs_error[[x]][[i]][[z]])[1])))) %>%
      mutate(rep = i, 
             type1_rate = rate_list$type1_rate[x], 
             type2_rate = rate_list$type2_rate[x]))))))) %>%
  rename(allele_count = V1) %>%
  mutate(locus = rep(seq(1, 20, by = 1), nrow(rate_list)*length(reps))) 


comparison_allele_count <- do.call(rbind, replicate(n= nrow(rate_list)*length(reps), expr = allele_counts_df, simplify = FALSE))

prop_loci_diff <- allele_counts_error_df %>%
  mutate(equal_to_orig = allele_count == comparison_allele_count) %>%
  group_by(across(all_of(c("rep", "type1_rate", "type2_rate")))) %>%
  summarise(prop_diff = 1 - (sum(equal_to_orig)/20))

prop_loci_diff %>%
  ggplot(aes(group = interaction(type1_rate, type2_rate))) +
  geom_boxplot(aes(y = prop_diff, x = interaction(type1_rate, type2_rate), color = interaction(type1_rate, type2_rate))) +
  #scale_color_viridis(discrete = TRUE) +
  xlab("type 1 rate x type 2 rate") +
  ylab("proportion of loci where error cause change in allele #") +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
  #facet_wrap(~ type1_rate * type2_rate)


# Get the proportion of the 3 least common alleles and identity of those alleles for all loci using custom get_df_of_sm_alleles function
rarest_allele_prop_df <- as.data.frame(do.call(rbind, lapply(allele_freqs_real, get_df_of_sm_alleles))) %>%
    mutate(locus = loci, 
           rep = "no error", 
           type1_rate = 0, 
           type2_rate = 0)

rarest_allele_prop_df_error <- do.call(rbind, lapply(rate_combo_count, function(i) {
  single_rate_rarest_alleles <- do.call(rbind, lapply(reps, function(x) {
  # Make a data frame for each rep that consist of the info for the rarest alleles at each locus and note rep and locus
    single_rep_rarest_alleles <- as.data.frame(do.call(rbind, lapply(allele_freqs_error[[i]][[x]], get_df_of_sm_alleles))) %>%
      mutate(locus = loci, 
             rep = x, 
             type1_rate = rate_list$type1_rate[i], 
             type2_rate = rate_list$type2_rate[i])
    return(single_rep_rarest_alleles)
  }))
return(single_rate_rarest_alleles)
}))


# Merge the real and error dfs for visualization purposes
rarest_allele_prop_df %>% 
  rbind(rarest_allele_prop_df_error) %>%
  mutate(rep = as.factor(rep)) %>%
  # Make the graph of the freq of the 3 rarest alleles w/ and without error
  ggplot(aes(x = locus)) +
  geom_point(aes(y = smallest_prop, color = rep), alpha = .4) + 
  geom_point(aes(y = sec_smallest_prop, color = rep), alpha = .4) + 
  geom_point(aes(y = third_smallest_prop, color = rep), alpha = .4) + 
  #scale_color_manual(values = c(rep("red", length(reps)), "black")) +
  facet_wrap(~ type1_rate * type2_rate) +
  ylab("Frequency of alleles") +
  labs(title = "Frequency of 3 rarest alleles in the population with and without error") +
  theme_classic()

rarest_allele_prop_df_error %>% 
  group_by(locus, type1_rate, type2_rate) %>%
  summarise_at(vars("smallest_prop", "sec_smallest_prop", "third_smallest_prop"), median) %>%
  rbind(subset(rarest_allele_prop_df, select = -c(rep, smallest_allele_num, sec_smallest_allele_num, third_smallest_allele_num))) %>%
  # Make the graph of the freq of the 3 rarest alleles w/ and without error
  ggplot(aes(x =  interaction(type1_rate, type2_rate))) +
  geom_point(aes(y = smallest_prop, color =  interaction(type1_rate, type2_rate))) + 
  geom_point(aes(y = sec_smallest_prop, color =  interaction(type1_rate, type2_rate))) + 
  geom_point(aes(y = third_smallest_prop, color =  interaction(type1_rate, type2_rate))) + 
  facet_wrap(~locus) +
  theme_classic() +
  scale_color_manual(values = c("red", rep("black", length(rate_combo_count)))) +
  ylab("Frequency of alleles") +
  labs(title = "Frequency of 3 rarest alleles in the population with and without error by locus") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


#### Want to make a graph that shows the number of reps that have changed the identity of the 1st 2nd and 3rd least frequent allele !####



###Subsetting the initial dataset to get wild + garden dfs

perc_sampled <- .2
pops <- c(unique(genetic_data$Pop)) # Get the number of populations in the data_frame
genetic_data$location <- "wild" # Add a location column to the genetic data df
garden_inds <- unlist(lapply(pops, function(x) {
  inds <- genetic_data$Ind[genetic_data$Pop == x]
  garden <- sample(inds, size = perc_sampled * length(inds), replace = F) # Samples the specified percent of each pop
  return(garden)
}))

#I could return the above as a list that matches the index of the perc_sampled run.... and that way I can attach the garden inds within a for loop or apply type thing instead of continuously remaking the entire genetic dataset? --> it doesn't seem overly slow to edit these dfs with lapplys as done below

genetic_data[genetic_data$Ind %in% garden_inds, c("location")] <- "garden"

# Making the error datasets have all the info present in the real dataset 
genetic_data_error <- lapply(genetic_data_sim_error, function(i) {
  all_genetic_data_reps <- lapply(i, function(x) {
  single_genetic_data_error_rep <- x %>%
  mutate(location = genetic_data$location, 
         Ind = genetic_data$Ind, 
         Pop = genetic_data$Pop)
  return(single_genetic_data_error_rep)})
  return(all_genetic_data_reps)})


# Now compare the wild data and the garden data (allele count, allele frequency in wild, heterozygosity at each allele in population)


# Get the lists that have all of the info I need for graphing
w_g_real_alleles_and_freqs <- get_alleles_and_freqs_w_g(genetic_data)
w_g_error_alleles_and_freqs <- lapply(genetic_data_error, function(x)
  lapply(x, get_alleles_and_freqs_w_g))


# Make the full data frames that holds all summary info across all alleles/loci
full_dfs_in_lists <- lapply(loci, compare_w_g, w_g_real_alleles_and_freqs, w_g_error_alleles_and_freqs)

full_dfs_in_lists <- lapply(w_g_error_alleles_and_freqs, function(x) {
  lapply(x, function(i) {
  lapply(loci, compare_w_g, w_g_real_alleles_and_freqs, i)})})


# Data frame with all of the info at allele level
allele_info_df <- do.call(rbind, lapply(rate_combo_count, function(z){
  all_allele_info_per_rate_combo <- do.call(rbind, lapply(reps, function(i){
    all_allele_info_per_rep <- do.call(rbind, lapply(loci, function(x) {
    allele_info_per_rep <- full_dfs_in_lists[[z]][[i]][[x]][[1]] %>%
      mutate(rep = i, 
             type1_rate = rate_list$type1_rate[z], 
             type2_rate = rate_list$type2_rate[z]) 
    return(allele_info_per_rep)
    }))
  return(all_allele_info_per_rep)
  }))
return(all_allele_info_per_rate_combo)
}))


# Data frame with all of the info at locus level
loci_info_df <- do.call(rbind, lapply(rate_combo_count, function(z){
  all_loci_info_per_rate_combo <- do.call(rbind, lapply(reps, function(i){
    all_loci_info_per_rep <- do.call(rbind, lapply(loci, function(x) {
      loci_info_per_rep <- full_dfs_in_lists[[z]][[i]][[x]][[2]] %>%
        mutate(rep = i, 
               type1_rate = rate_list$type1_rate[z], 
               type2_rate = rate_list$type2_rate[z])
      return(loci_info_per_rep)
    }))
    return(all_loci_info_per_rep)
  }))
  return(all_loci_info_per_rate_combo)
}))

loci_info_df_for_graphing <- loci_info_df%>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error"), 
         type1_rate = ifelse(dataset %in% c("wild_real", "garden_real"), 0, type1_rate), 
         type2_rate = ifelse(dataset %in% c("wild_real", "garden_real"), 0, type2_rate),
         dataset_rates = ifelse(dataset %notin% c("wild_error" , "garden_error"), 
                                paste0(dataset),
                                paste0(dataset, "_", type1_rate, "_", type2_rate)), 
         w_g = ifelse(dataset %in% c("garden_error", "garden_real"), "g", "w"), 
         rep= as.factor(rep)) %>%
  mutate(dataset_rates = fct_relevel(dataset_rates, "garden_real", "wild_real", after = 0)) %>%
  distinct(.keep_all = TRUE)


allele_info_df_for_graphing <- allele_info_df %>%
  mutate(dataset = fct_relevel(dataset, "wild_real", "garden_real", "wild_error", "garden_error"), 
         dataset_rates = ifelse(dataset %notin% c("wild_error" , "garden_error"), 
                              paste0(dataset),
                              paste0(dataset, "_", type1_rate, "_", type2_rate)), 
         w_g = ifelse(dataset %in% c("garden_error", "garden_real"), "g", "w"), 
         rep= as.factor(rep)) %>%
  mutate(dataset_rates = fct_relevel(dataset_rates, "garden_real", "wild_real", after = 0))

#### Making some graphs####

#need to adjust these graphs to be able to look across rate combos!


allele_info_df_for_graphing  %>%
  group_by(all_alleles, locus, type1_rate, type2_rate, dataset) %>%
  summarise_at(vars("allele_freq"), median) %>%
#need to modify above summarise_at command to get median values of allele freqs across reps for each rate combo --> not sure this will work..... might need to think of another route....
  filter(locus == 1) %>%
  ggplot(alpha = .7) +
  geom_point(aes(x = all_alleles, y = allele_freq, color = dataset), shape = 1) +
  scale_color_manual(values = c("indianred", "gray",  "red", "black")) +
  facet_wrap(~ type1_rate * type2_rate) +
  theme_classic()


allele_info_df_for_graphing  %>%
  ggplot() +
  geom_boxplot(aes(y = allele_freq, x = dataset, color = dataset), outlier.shape = NA) +
  geom_point(aes(y = allele_freq, x = dataset, color = dataset), alpha = .3) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_wrap(~ type1_rate * type2_rate) +
  theme_classic()

# Note an overly important graph, just shows differences between expected and obs het for each data set split by each locus, worse with many rates
loci_info_df_for_graphing  %>%
  group_by(locus, type1_rate, type2_rate, dataset, dataset_rates) %>%
  summarise_at(vars("exp_het", "obs_het"), median) %>%
  ggplot() +
  geom_point(aes(y = exp_het, x = dataset_rates, color = dataset_rates)) +
  geom_point(aes(y = obs_het, x = dataset_rates, color = dataset_rates)) +
  geom_segment(aes(y = obs_het, yend = exp_het, x = dataset_rates, xend= dataset_rates), color = "black") +
  scale_color_manual(values = c("gray", "black", rep("indianred", length(rate_combo_count)), rep("red", length(rate_combo_count)))) +
  facet_wrap(~ locus) +
  theme_classic() +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank())


# Comparison of observed heterozygosity of all loci in each dataset (at higher type 2 error rates, this seems to get higher.... makes sense bc literally making more heterozygotes)
loci_info_df_for_graphing  %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = as.factor(type1_rate), color = as.factor(type2_rate))) +
  #scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_wrap(~ w_g) +
  theme_classic()


# Comparison of expected heterozygosity of all loci in each dataset (error is lower bc more alleles at a lower frequency)
loci_info_df  %>%
  ggplot() +
  geom_boxplot(aes(y = exp_het, x = dataset, color = dataset)) +
  scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  facet_wrap(~ type1_rate * type2_rate) +
  theme_classic()

# Plot expected het from the wild datasets of all combos of error rates against next to the single real
loci_info_df  %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = exp_het, x = interaction(type2_rate, type1_rate), color = type1_rate)) +
  #scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  theme_classic()

# Plot observed het from the wild datasets of all combos of error rates against next to the single real --> shows that very little change with increased type 1 rate but small increase as inc type 2 rate 
loci_info_df_for_graphing  %>%
  filter(dataset %notin% c("garden_error", "garden_real")) %>%
  ggplot() +
  geom_boxplot(aes(y = obs_het, x = interaction(type1_rate, type2_rate), color = type2_rate)) +
  #scale_color_manual(values = c("black","gray", "red",  "indianred")) +
  theme_classic()


# Comparison of number of across all loci in each dataset (as type 1 rate increases --> exp het decreases, as type 2 rate increases --> obs het increases)
loci_info_df_for_graphing  %>%
  ggplot() +
  geom_boxplot(aes(y = num_alleles, x = interaction(type1_rate, type2_rate), color = type2_rate), outlier.shape = NA) +
  geom_jitter(aes(y = num_alleles, x = interaction(type1_rate, type2_rate), color = type2_rate), alpha = .3) +
  facet_wrap(~ w_g) +
  theme_classic()


# Across all loci, proportion of wild rare alleles captured in garden samples in error vs real dataset (rare loci are less likely to be represented in the garden sample in the error dataset)
loci_info_df_for_graphing %>%
  ggplot() +
  geom_boxplot(aes(y = prop_rare_wild_alleles_captured, x = interaction(type1_rate, type2_rate), color = interaction(type1_rate, type2_rate)), outlier.shape = NA) +
  geom_jitter(aes(y = prop_rare_wild_alleles_captured, x = interaction(type1_rate, type2_rate), color = interaction(type1_rate, type2_rate)), alpha = .3) +
  theme_classic()


loci_info_df_for_graphing %>%
  filter(dataset %in% c("wild_real", "wild_error")) %>%
  group_by(locus, type1_rate, type2_rate, dataset, dataset_rates) %>%
  summarise_at(vars("rare_allele_count"), median) %>%
  ggplot(aes(x = rare_allele_count)) +
  geom_bar(aes(fill = dataset), position=position_dodge2(preserve = "single")) +
  #scale_fill_manual(values = c("black","gray", "red",  "indianred")) +
  scale_fill_manual(values = c("black","red")) +
  ylab("Number of loci with given amount of rare alleles") +
  xlab("Number of rare alleles") +
  theme_classic() +
  facet_wrap(~ type1_rate * type2_rate)

loci_info_df_for_graphing %>%
  filter(dataset %in% c("garden_real", "garden_error")) %>%
  group_by(locus, type1_rate, type2_rate, dataset, dataset_rates) %>%
  summarise_at(vars("rare_allele_count"), median) %>%
  ggplot(aes(x = rare_allele_count)) +
  geom_bar(aes(fill = dataset), position=position_dodge2(preserve = "single")) +
  #scale_fill_manual(values = c("black","gray", "red",  "indianred")) +
  scale_fill_manual(values = c("gray","indianred")) +
  ylab("Number of loci with given amount of rare alleles") +
  xlab("Number of rare alleles") +
  theme_classic() +
  facet_wrap(~ type1_rate * type2_rate)


loci_info_df_for_graphing %>%
  ggplot(aes(y = rare_allele_count)) +
  geom_boxplot(aes(x = interaction(type1_rate, type2_rate), color = as.factor(type2_rate))) +
  geom_jitter(aes(x = interaction(type1_rate, type2_rate), color = as.factor(type2_rate)), alpha = .3) +
  ylab("Number of rare alleles/locus") +
  xlab("Type 1 and Type 2 error rates") +
  theme_classic() +
  facet_wrap(~ w_g)


