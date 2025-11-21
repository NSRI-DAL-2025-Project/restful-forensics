# Author: Ambrocio Melvin A. Matias
# Contact Detail: ambrocio.matias@uq.net.au / dinmatias2006@gmail.com
# Date Edited: 27/07/2016

# Description: This code computes genetic distances and use these
# to determine barcoding gaps. Other summary statistics are also computed
# The code uses several packages such as 
# ape for the genetic distance computation, 
# phangorn in loading the file, though this can be replaced by ape functions
# tidyr for restructuring the format of distance data
# ggplot2 for the plot
# Final note, the major contribution of the author here is the code
# for the restructuring the distance data, which is tedious to perform
# manually.

# NOTE: I deemed it easier to use this one with R Studio as platform of R

# Requirements: First, you need to install R (https://www.r-project.org/)
# and, a higly recommended interface, R Studio (https://www.rstudio.com/).

# After installation, this script can be viewed in R Studio
# First, highlight the commands below starting with
# install.packages and then hit Ctrl + Enter (You need internet for this)

# install.packages('phangorn')
# install.packages('ape')
# install.packages('ggplot2')
# install.packages('tidyr')
# install.packages('magrittr')

# This will install the packages need for the script

library('phangorn')
library('ape')
library('ggplot2')
library('tidyr')
library('magrittr')
library('pegas')
library('adegenet')


# Clears the R environment before running the script
rm(list=ls())


#####################
#### INPUTS START ####
######################

## Revision 26/07/2016
## Rather than having input here,
## why not just (1) create a project
## (2) have a data folder where
## the fasta file can be placed
## (3) source every file in that one
## (4) have a text file input for
## the model and filename of LABELS
## Revision 26/07/2016

# Place the address of the directory where 
# the aligned fasta file is located. 
#file_dir <- '/home/aamatias/OneDrive/1 Matias Data/1 Half Beaks/'

## Instead of working file_dir use: 
## list.files("Data/")[grep(".fas", list.files("Data/"))]
## list.files("Data/")[grep(".csv", list.files("Data/"))]

# Place the name of the fasta file.
file_fas <-  list.files("Data/")[grep(".fas", list.files("Data/"))]


file("Model/model.txt") %>% 
  readLines(1) %>%
  strsplit(':') %>%
  unlist() -> model_in

model_in <- model_in[2]

# Choices for models:
# "raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", 
# "F84", "BH87", "T92", "TN93", "GG95", "logdet"
# "paralin", "indel", or "indelblock".

# Place the address of the csv containing the label (species name)
# of each sequence
label_in <- list.files("Data/")[grep(".csv", list.files("Data/"))]

# IMPORTANT: 
# The order of the entry in this fiel is important.
# It should be the same as the order of the sequences
# of specimens.
# Another note, there must only be genus and species
# name in the file. If you want to indicate a subspecies
# or a variety, place an underscore in the subspecies name
# e.g., Zenarchopterus philippinus var magatensis should be
# written as Zenarchopterus philippinus_var_magatensis

######################
#### INPUTS END   ####
######################

# file_name = complete directory and filename of fasta file
file_name <- paste("Data/", file_fas, sep ='')

# Load the aligned file using read.phyDat then convert it to DNAbin
# data.align = variable containing the loaded file
# This is the one that should be aligned
read.phyDat(file_name, format = 'fas', type = 'DNA') %>% 
  as.DNAbin() -> data_align

# 2) Compute Distance
# Compute pairwise distances between the samples

# Then compute the distance
data_dist_matrix <- dist.dna(x = data_align, model = model_in)
# Variance can be included here

# Note: The succeeding section is actually the tricky
# part done in excel. Here, it will simplified by using
# functions that 'melt' the distance matrix into a
# dataframe with long format
data_final <- data.frame(t(combn(attr(data_align, "dimnames")[[1]], 2)), 
                         as.numeric(data_dist_matrix), 
                         stringsAsFactors = FALSE)

colnames(data_final) <- c('taxa_1', 'taxa_2', 'distance')


paste("Data/", label_in, sep = "") %>% 
  read.csv2(stringsAsFactors = FALSE) -> data_labels

data_labels$Species <- as.character(data_labels$Species)



label_final <- data.frame(t(combn(data_labels$Species, 2)), 
                          stringsAsFactors = FALSE)

colnames(label_final) <- c("taxon_1", "taxon_2")


data_final <- data.frame(taxon_1_UI = label_final$taxon_1,
                         taxon_1_LF = data_final$taxa_1,
                         taxon_2_UI = label_final$taxon_2,
                         taxon_2_LF = data_final$taxa_2,
                         Distance = data_final$distance,
                         stringsAsFactors = FALSE)

# http://stackoverflow.com/questions/4350440/split-a-column-of-a-data-frame-to-multiple-columns  
sci_names_1 <- do.call(rbind, strsplit(data_final$taxon_1_UI, split=" "))

sci_names_2 <- do.call(rbind, strsplit(data_final$taxon_2_UI, split=" "))

sci_names <- data.frame(sci_names_1, 
                        sci_names_2,
                        stringsAsFactors = FALSE)


# Replace column names for clarity
colnames(sci_names) <- c('Genus_1', 'species_1', 'Genus_2', 'species_2')

# Concatenates original data with the genus and species
# names
data_final <- cbind(data_final, sci_names)

# This is a sweet function.
# The function basicall compares invidual rows of the two vectors
# specified below.
# This line basically determines if the distance is an interspecific
# or intraspecific comparison
# Intraspecific: Genus_1 == Genus_2 & Species_1 == Species_2
# Interspecies: Genus_1 == Genus_2 & Species_1 != Species_2
# Intergenre: Genus_1 != Genus_2
# marker = a variable to indicate if the distance is inter or intraspecific
marker <- ifelse(data_final$Genus_1 == data_final$Genus_2 & 
                   data_final$species_1 == data_final$species_2, 
                 'Intraspecific', ifelse(data_final$Genus_1 == data_final$Genus_2 & 
                                           data_final$species_1 != data_final$species_2, 
                                         'Interspecific', 'Intergenre'))

data_final <- cbind(data_final, marker)
#### Revisions 30/07/2016

data_final$marker <- factor(data_final$marker, levels(data_final$marker)[c(3,2,1)])

# Output Directory
# out_dir <- paste(file_dir, 'Output_Barcoding', sep ='')
# Replaced with just the Output_Barcoding beacause of
# the R project
dir.create('Output_Barcoding', 
           showWarnings = TRUE, 
           recursive = FALSE, mode = "0777")

### Need to revise
# Need to subset the data here
results_bargap <- data_final[data_final$marker == 'Interspecific' | data_final$marker == 'Intraspecific', ]
### Need to revise

# Plot for comparing Distance from different comparisons
results_groupplot <- ggplot(data_final, aes(x = marker, y = Distance)) + geom_boxplot() + 
  xlab('Type of Comparison') + ylab('Distance') + theme_bw() + 
  ggtitle('Distance vs Type of Comparison') + 
  guides(fill = guide_legend(title = NULL))

ggsave(filename = 'distances_group.pdf',
       path = 'Output_Barcoding/',
       plot = results_groupplot, device = pdf, height = 8.27, width = 11.69)

# Plots the barcoding gap
results_gapplot <- ggplot(results_bargap, aes(x = Distance, fill = as.factor(marker))) + 
  geom_histogram(alpha = 0.5, position = 'identity', aes(y=(..count..)/sum(..count..))) +
  xlab('Pairwise Genetic Distance') + ylab('Frequency') + theme_bw() + 
  ggtitle('IntraSpecific vs InterSpecific Distances') +
  guides(fill = guide_legend(title = NULL))

ggsave(filename = 'barcoding_gap.pdf',
       path = 'Output_Barcoding/',
       plot = results_gapplot, device = pdf, height = 8.27, width = 11.69)


# Note: This parameter aes(y=(..count..)/sum(..count..)) can be removed to revert
# back to count rather than frequency

# Computes the minimum, maximum and mean of the inter 
# and intraspecific distances
data_inter_summary <- c(min(results_bargap[results_bargap$marker == 'Interspecific', 5]), 
                        max(results_bargap[results_bargap$marker == 'Interspecific', 5]), 
                        mean(results_bargap[results_bargap$marker == 'Interspecific', 5]))


data_intra_summary <- c(min(results_bargap[results_bargap$marker == 'Intraspecific', 5]), 
                        max(results_bargap[results_bargap$marker == 'Intraspecific', 5]), 
                        mean(results_bargap[results_bargap$marker == 'Intraspecific', 5]))

# Combines the summary of the inter and intra
results_summary <- cbind(c('minimum', 'maximum', 'mean'), data_inter_summary, data_intra_summary)
colnames(results_summary) <- c('Measure', 'InterSpecific', 'Intraspecific')

# Relabel the rows of the final data table
# row.names(results_summary) <- (c('minimum', 'maximum', 'mean'))



# Series of commands to save the output as csv file
write.table(data_final[,c(1,2,3,4,5,10)], 
            file = 'Output_Barcoding/Distances.Melted.txt',
            dec = '.', sep = '\t', row.names = FALSE)
write.table(data_dist_matrix, 
            file = 'Output_Barcoding/Distances.Raw.txt',
            dec ='.', sep = '\t', row.names = FALSE)
write.table(results_summary, 
            file = 'Output_Barcoding/Distance.Summary.txt', 
            dec ='.', sep = '\t', row.names = FALSE)


### EDIT END 23/06/2016
#######################
## Revision 26/07/2016:
## MAJOR CHANGE:
## Replaced . with _ of each variable/object name
## Revision 31/07/2016:
## Added boxplot of distances vs category of comparison

