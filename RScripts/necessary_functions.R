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


####Creating functions for importing initial MSAT data####

# Create a function that is the opposite of %in%
`%notin%` <- Negate(`%in%`)

#A function that converts an arlequin simulation file to a genepop file, edited from Austin's code bc I am not converting the genind files to editable data
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
simulating_error <- function(data_for_editing, syst_err_rate, stoch_err_rate, msat_length, replicate_number){
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
    e_1 = syst_err_rate / (1 + syst_err_rate)
    #e_2 = stoch_err_rate / (k - 1) # Old error rate, when type 2 error = change to any other already known geno
    e_2 = stoch_err_rate / 2 # Dividing by 2 bc in the new version of type 2 error the real geno can become either of 2 genos
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
           syst_err_rate = syst_err_rate, 
           stoch_err_rate = stoch_err_rate)
  diffs_in_genos <- diffs_in_genos %>%
    mutate(rep = replicate_number,
           syst_err_rate = syst_err_rate,
           stoch_err_rate = stoch_err_rate)
  
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
    # Make columns with observed heterozygosity
    mutate(wild_real = locus_real_data[[3]],  # Third item in locus_data list is a the number observed heterozygotes in the wild
           garden_real = locus_real_data[[6]], # Sixth item in locus_data list is a the number observed heterozygotes in the garden
           wild_error = locus_error_data[[3]],  # Third item in locus_data list is a the number observed heterozygotes in the wild
           garden_error = locus_error_data[[6]]) %>% # Sixth item in locus_data list is a the number observed heterozygotes in the garden 
    # Make above tidy
    pivot_longer(cols = c(wild_real, garden_real, wild_error, garden_error), names_to = "dataset", values_to = "obs_het") %>%
    
    # Calc expected heterozygosity with 1 - (exp frequencies of all homozygotes based on frequencies of each allele)
    mutate(exp_het = ifelse(dataset == "wild_real", 1-sum(na.omit(subset(alleles_df, dataset == "wild_real", select = c(allele_freq)))^2),
                            ifelse(dataset == "garden_real", 1-sum(na.omit(subset(alleles_df, dataset == "garden_real", select = c(allele_freq)))^2),
                                   ifelse(dataset == "wild_error", 1-sum(na.omit(subset(alleles_df, dataset == "wild_error", select = c(allele_freq)))^2), 1-sum(na.omit(subset(alleles_df, dataset == "garden_error", select = c(allele_freq)))^2))))) %>%
    
    # Calc F (inbreeding coefficient) with F = (Exp Het - Obs Het) / Exp Het 
    mutate(inbreeding_coeff = (exp_het - obs_het) / exp_het) %>%
    
    # Add column with number of alleles per dataset
    mutate(num_alleles = ifelse(dataset == "wild_real", length(wild_real_alleles), 
                                ifelse(dataset == "garden_real", length(garden_real_alleles),  
                                       ifelse(dataset == "wild_error", length(wild_error_alleles), length(garden_error_alleles)))),
           #Need to fix code below such that if there are 0 of a type of allele in the wild, representation will be 1?
           prop_rare_wild_alleles_captured = ifelseifelse(dataset == "wild_real", sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Rare" & alleles_df$in_garden_real == T)/sum(grepl("wild_real", alleles_df$dataset) & alleles_df$allele_freq_cat == "Rare"),
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

