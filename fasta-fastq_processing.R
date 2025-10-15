# processing of fasta files

BiocManager::install("Biostrings") # to read the fasta files
library(Biostrings)

# read fasta
dna_seqs = readDNAStringSet("sample.fasta")
#QUESTION: MULTIPLE SEQUENCES??
## Other options, use subseq() to subset a certain length
## Use reverseComplement() to generate the complement strand
## use translate() to translate DNA to proteins


# NOTE: Install latex system for msaPrettyPrint() in Dockerfile
BiocManager::install("ggtree")
BiocManager::install("pwalign")
BiocManager::install("msa") # to perform msa using standard alignment algorithms
BiocManager::install() # to automatically install all dependencies

library("pwalign")
library(msa)

####################
# Installing latex #
####################



aligned_sequences = msa(dna_seqs, "Muscle") # options are "ClustalOmega", "Muscle", and "ClustalW"
# add show="complete" when printing results

alignment_scores = msaConservationScore(aligned_sequences, substitutionMatrix = BLOSUM62) # substitution matrix can be changed

# download aligned sequences
sink("aligned_seqs.fasta") #redirects all outputs to the specified file, sink() also creates a new file
print(aligned_sequences, show = "complete")
sink()

# Print a color-coded MSA
msaPrettyPrint(aligned_sequences, output="pdf", showNames="none", showLogo = "none")

# for trimming, use a linux-based program
# compilation instructions: https://trimal.readthedocs.io/en/latest/installation.html
#documentation: https://trimal.readthedocs.io/en/latest/usage.html

# the tree
# use ape for building the tree
# tidy tree to calculate based on different methods


####################
# ACTUAL FUNCTIONS #
####################

read_fasta <- function(zipped, directory){
   utils::unzip(zipped, 
                files = NULL, 
                list = FALSE, 
                overwrite = TRUE, 
                exdir = file.path(directory, "fasta_files"))
   
   data_path <- file.path(directory, "fasta_files")
   fasta_patterns <- paste("\\.fasta$", "\\.fa$", "\\.fna$", "\\fas$", sep = "|")
   fasta_files <- list.files(path = data_path, pattern = fasta_patterns, full.names = TRUE)

   dna_sequences <- Biostrings::readDNAStringSet(fasta_files)
   
   return(dna_sequences)
}

# requirements Latex and pdflatex
# Option 1: install latex using: sudo apt-get install textlive-full


msa_results <- function(files, algorithm, directory){
   
   # Creating Substitution Matrix
   personal_matrix <- pwalign::nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = TRUE, type = "DNA")
   gap_penalty <- -2
   dna_matrix_wgaps <- rbind(personal_matrix, gap_penalty)
   dna_matrix_wgaps <- cbind(dna_matrix_wgaps, gap_penalty)
   rownames(dna_matrix_wgaps)[nrow(dna_matrix_wgaps)] <- "-"
   colnames(dna_matrix_wgaps)[ncol(dna_matrix_wgaps)] <- "-"
   colnames(dna_matrix_wgaps) <- c("A", "C", "G", "T", "-")
   rownames(dna_matrix_wgaps) <-  c("A", "C", "G", "T", "-")
   dna_matrix_wgaps <- as.matrix(dna_matrix_wgaps)
   
   filename1 <- paste0(directory, "/aligned_seqs.txt")
   filename2 <- paste0(directory, "/aligned_seqs_wscores.txt")
   
   # perform msa
   aligned_sequences <- msa::msa(files,substitutionMatrix = dna_matrix_wgaps, method = algorithm) # ClustalW, ClustalOmega, MUSCLE
   output_aligned <- utils::capture.output(print(aligned_sequences, show = "complete"))
   writeLines(output_aligned, filename1)
   
   # calculate alignment score
   alignment_scores <- msa::msaConservationScore(aligned_sequences, substitutionMatrix = dna_matrix_wgaps)
   output_scores <- utils::capture.output(print(alignment_scores, show = "complete"))
   writeLines(output_scores, filename2)
   
   
   # saving the aligned sequences
   #sink(filename1)
   #print(aligned_sequences, show = "complete")
   #sink()
   
   #saving alignment with scores
   #sink(filename2)
   #print(alignment_scores, show = "complete")
   #sink()
   
   filename3 <- paste0(directory, "/aligned_seqs.pdf")
   # saving a pdf file
   msa::msaPrettyPrint(aligned_sequences, output="pdf", file = filename3, showNames= "none", showLogo = "none")
   
   # double check directory where this is saved
   return(list(
      fasta_file = filename1,
      msa_pdf = filename2,
      aligned = aligned_sequences
   ))
}


################
# VCF TO FASTA #
################

vcf_to_fasta <- function(vcf_file, reference, bcftools_path, directory){
   output_file <- paste0(directory, "consensus.fa")
   
   bcftools_path(stringr::str_c(
      "consensus -f ", reference, " ", vcf_file, " -o ", output_file
   ))
}

#####################
# Building the tree #
#####################

building_distance <- function(alignment, model){
   bins <- ape::as.DNAbin(alignment)
   distance <- ape::dist.dna(bins, model = model)
   # ape: nj or upgma
   # phangorn: maximum likelihood or parsimony
   nj_tree <- ape::nj(distance)
   upgma_tree <- ape::upgma(distance)
   
   # sample_tree <- ape::rcoal()
   ml_tree <- phangorn::pml() # what about other options?
                              # accepts alignment?? docu shows it accepts a tree, class phylo
}

# model for dist.dna
#"raw", "N", "TS", "TV", "JC69", "K80" (the default), "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", or "indelblock".

build_nj_tree <- function(alignment, outgroup = NULL, seed = 123, model = model){
   library(ape)
   library(ggtree)
   bins <- ape::as.DNAbin(alignment)
   distance <- ape::dist.dna(bins, model = model)
   nj_tree <- ape::nj(distance)
   
   # Rooting
   if (!is.null(outgroup) && outgroup %in% nj_tree$tip.label){
      nj_tree <- ape::root(nj_tree, outgroup = outgroup)
   }
   
   nj_tree <- ape::ladderize(nj_tree)
   
   num_sites <- ncol(bins)
   if (num_sites < 10){
      warning("Alignment has fewer than 10 sites. Skipping bootstrap.")
      
      tree_plot <- ggtree(nj_tree, branch.length = "none") +
         theme_tree2() +
         geom_tiplab() +
         ggtitle("NJ Tree")
      
      return(tree_plot)
   } 
   
   # Bootstrapping
   set.seed(seed)
   boots <- ape::boot.phylo(nj_tree, bins, 
                            FUN = function(x){
                               tree <- ape::nj(ape::dist.dna(x, model = model))
                               if (!is.null(outgroup) && outgroup %in% tree$tip.label){
                                  tree <- ape::root(tree, outgroup = outgroup)
                               }
                               ape::ladderize(tree)
                            }, rooted = TRUE
                            )
   
   boots[is.na(boots)] <- 0
   nj_tree$node.label <- as.character(boots)
   tree_plot <- ggtree(nj_tree, branch.length = "none") +
      theme_tree2() +
      geom_tiplab() +
      geom_text2(aes(subset = !isTip, label = label), hjust = -0.3) +
      ggtitle("NJ Tree") +
      xlim(0, 20)
   
   return(tree_plot)
}

build_upgma_tree <- function(alignment, outgroup = NULL, seed =123, model = model, ){
   library(ape)
   
   bins <- ape::as.DNAbin(alignment)
   distance <- ape::dist.dna(bins, model = model)
   upgma_tree <- upgma(distance)
   
   if (!is.null(outgroup) && outgroup %in% upgma_tree$tip.label){
      upgma_tree <- ape::root(upgma_tree, outgroup = outgroup)
   }
   
   upgma_tree <- ape::ladderize(upgma_tree)
   
   num_sites <- ncol(bins)
   if (num_sites < 10){
      warning("Alignment has fewer than 10 sites. Skipping bootstrap.")
      
      tree_plot <- ggtree(upgma_tree, branch.length = "none") +
         theme_tree2() +
         geom_tiplab() +
         ggtitle("UPGMA Tree")
      
      return(tree_plot)
   } # end of num_sites check
   
   
   set.seed(seed)
   boots <- ape::boot.phylo(upgma_tree, bins, 
                            FUN = function(x){
                               tree <- ape::nj(ape::dist.dna(x, model = model))
                               if (!is.null(outgroup) && outgroup %in% tree$tip.label){
                                  tree <- ape::root(tree, outgroup = outgroup)
                               }
                               ape::ladderize(tree)
                            }, rooted = TRUE, 
   )
   
   boots[is.na(boots)] <- 0
   upgma_tree$node.label <- as.character(boots)
   tree_plot <- ggtree(upgma_tree, branch.length = "none") +
      theme_tree2() +
      geom_tiplab() +
      geom_text2(aes(subset = !isTip, label = label), hjust = -0.3) +
      ggtitle("UPGMA Tree") +
   xlim(0, 20)
   
   return(tree_plot)
}


build_max_parsimony <- function(alignment, outgroup = NULL, use_midpoint = NULL, seed = 123, directory){
   library(phangorn)
   library(ape)
   
   bins <- ape::as.DNAbin(alignment)
   phy <- phangorn::phyDat(bins, type = "DNA")
   dm <- dist.ml(phy)
   start_tree <- NJ(dm)
   parsimony_tree <- optim.parsimony(start_tree, phy)
   
   # Rooting
   if (!is.null(outgroup) && outgroup %in% parsimony_tree$tip.label){
      parsimony_tree <- root(parsimony_tree, outgroup = outgroup, resolve.root = TRUE)
   } else if (use_midpoint){
      parsimony_tree <- midpoint(parsimony_tree)
   }
   
   # boostrapping
   set.seed(seed)
   bs_pars <- bootstrap.phyDat(phy, \(x) optim.parsimony(NJ(dist.ml(x)), x))
   
   # plot
   filename <- paste(directory, "parsimony_tree.png")
   png(filename, width = 800, height = 600)
   plotBS(parsimony_tree, bs_pars, main = "Parsimony Tree")
   dev.off()
   
   return(filename)
}

build_ml_tree <- function(alignment, 
                          outgroup = NULL, 
                          use_midpoint = FALSE, 
                          seed = 123, 
                          bs_reps = 100,
                          directory){
   library(phangorn)
   library(ape)
   
   bins <- ape::as.DNAbin(alignment)
   phy <- phyDat(bins, type = "DNA")
   
   dm <- dist.ml(phy)
   start_tree <- NJ(dm)
   
   fit <- pml(start_tree, data = phy)
   # find best-fit model
   model_test <- modelTest(phy, tree = start_tree)
   best_model <- model_test$Model[which.min(model_test$BIC)]
   fit_opt <- optim.pml(fit, model = best_model, optGamma = TRUE, optInv = TRUE, rearrangement = "stochastic")
   
   tree <- fit_opt$tree
   if (!is.null(outgroup) && outgroup %in% tree$tip.label){
      tree <- root(tree, outgroup = outgroup, resolve.root=TRUE)
   } else if (use_midpoint){
      tree <- midpoint(tree)
   }
   
   # bootstrapping
   set.seed(seed)
   bs <- bootstrap.pml(fit_opt, bs = bs_reps, optNni = TRUE)
   
   filename <-  paste(directory, "ml_tree.png")
   png(filename, width = 800, height = 600)
   plotBS(tree, bs, main = paste("ML Tree (", best_model, ")"))
   dev.off()
   return(list(
      best_model = best_model,
      filename = filename
   ))
   
}