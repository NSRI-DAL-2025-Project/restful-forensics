# TO DO #1: (DONE) ADDITIONAL ARGUMENTS FOR PLINK - OPTIONAL FILTERING PROCEDURES
# TO DO #2: ADD DOWNLOADABLE SAMPLE FILES - can do, create sample files

###========================FILE CONVERSION (SNIPPER)===============###

tosnipper <- function(input, references, target.pop = TRUE, population.name = NULL, markers = snps){
   if (!dir.exists(dir))  dir.create(dir, recursive = TRUE)
   
   if(!require("pacman")){
      install.packages("pacman")
   }
   
   pacman::p_load(pacman, tools, readr, readxl, vcfR, dplyr, purrr, plyr, openxlsx, install = TRUE)
   
   # load in input file
   if (tools::file_ext(input) == "csv"){
      input.file <- readr::read_csv(input)
      input.file <- input.file %>%
         rename(Ind = 1)
      
   } else if (tools::file_ext(input) == "xlsx"){
      input.file <- readxl::read_excel(input)
      input.file <- input.file %>%
         rename(Sample = 1)
   } else if (tools::file_ext(input) == "vcf"){
      vcf_file <- vcfR::read.vcfR(input)
      vcf_gt <- vcfR::extract.gt(vcf_file, return.alleles = TRUE)
      vcf_cols <- vcfR::getFIX(vcf_file)
      vcf_to_df <- data.frame(vcf_cols, vcf_gt)   
      names(vcf_to_df) <- sub('^X', '', names(vcf_to_df))
      
      drops  <- c("CHROM", "POS", "REF", "ALT", "QUAL", "FILTER") # remove unnecessary columns
      vcf_to_df <- vcf_to_df[ ,!(names(vcf_to_df) %in% drops)] 
      
      # transposed (rsID as column headers)
      vcf_to_excel <- data.frame(t(vcf_to_df))
      vcf_to_excel <- vcf_to_excel[-1,]
      Sample <- rownames(vcf_to_excel)
      input.file <- data.frame(Sample, vcf_to_excel)
   } else {
      stop("Not an xlsx, csv, or vcf file.")
   } # TO ADD: VCF FILE AS INPUT
   
   ## correct format
   tosnipper <- lapply(
      input.file, 
      function(x){
         gsub(pattern = "/", replacement = "", x = x, fixed = TRUE)
      }
   )
   
   tosnipper <- as.data.frame(tosnipper)
   
   # to check Sample column for compatibility
   if(class(tosnipper$Sample) != "character"){
      tosnipper$Sample <-  as.character(tosnipper$Sample)
   }
   
   
   # load in the reference
   if(tools::file_ext(references) == "xlsx"){
      reference <- readxl::read_excel(references) 
      reference <- reference %>%
         rename(Sample = 1)
   } else if(tools::file_ext(references) == "csv"){
      reference <- readr::read_csv(references)
      reference <- reference %>%
         rename(Sample = 1)
   } else {
      stop("Not an xlsx or csv file")
   }
   
   if("Sample" %in% names(reference)){
      library(dplyr) #in case there's an error loading the library
      
      if(class(reference$Sample) != "character"){
         reference$Sample <-  as.character(reference$Sample)
      }
      
      matched <- tosnipper %>% dplyr::left_join(reference, by = "Sample")
      last.col <- as.integer(ncol(matched)) # refer to the total number of columns
      sec.last <- last.col -1
      Superpop <- as.data.frame(matched[,last.col])
      Population <- as.data.frame(matched[,sec.last])
      
      Sample <- as.data.frame(matched$Sample)
      
      #this removes the new col added at the end of the df
      data <- as.data.frame(matched[,2:ncol(matched)-1])
      
      #to ensure that the Sample column is not included in the data df
      drops <- "Sample"
      data <- data[,!(names(data) %in% drops)]
      
      #merge in order
      to_excel <- dplyr::bind_cols(Population, Superpop, Sample, data)
      names(to_excel)[names(to_excel) == "matched[, last.col]"] <- "Superpop"
      names(to_excel)[names(to_excel) == "matched[, sec.last]"] <- "Population"
      names(to_excel)[names(to_excel) == "matched$Sample"] <- "Sample"
      
   } else {
      stop("No 'Sample' header in reference file")
   }
   
   
   # split per population
   tosnpr_split <- split(to_excel, to_excel$Population)
   
   library(purrr)
   tosnpr_split <- tosnpr_split %>% map(`rownames<-`, NULL)
   
   # add the no. to the label
   tosnpr_split <- lapply(
      tosnpr_split, 
      function(x){
         x$No <- rownames(x)
         as.data.frame(x)
         data <- as.data.frame(x[,2:ncol(x)-2])
         Sample <- x[,1]
         x <- dplyr::bind_cols(x$No, Sample, data)
      }
   )
   
   # merges all to a list
   merged <- plyr::ldply(tosnpr_split, data.frame)
   merged <- merged[,-c(1,3)]
   
   if(target.pop == TRUE){
      
      # Subset the target pop
      target <- merged[merged$Population == population.name,]
      target$snpr <- "0"
      
      non.target <- merged[merged$Population != population.name,]
      non.target$snpr <- "1"
      
      merged2 <- dplyr::bind_rows(target, non.target) 
   } else if(target.pop == FALSE) {
      
      merged2 <- merged
      merged2$snpr <- "1"
      
   } else{
      stop("Parameter target.pop is not stated.")
   }
   
   # add value
   # count number of samples
   sample.count <- as.integer(nrow(merged2))
   
   # get the total number of pops
   pop.only <- as.data.frame(merged2$Superpop)
   pop.only <- pop.only[!duplicated(pop.only), ]
   pop.only <- as.data.frame(pop.only)
   pop.count <- as.integer(nrow(pop.only))
   
   # to add the columns
   merged3 <- merged2[,-2]
   merged3 <- as.data.frame(merged3)
   names(merged3)[names(merged3) == "...1"] <- sample.count
   names(merged3)[names(merged3) == "Superpop"] <- markers
   names(merged3)[names(merged3) == "Sample"] <- pop.count
   names(merged3)[names(merged3) == "snpr"] <- ""
   
   # merged3[nrow(merged3)+4,] <- NA #or
   merged3 <- rbind(NA, merged3)
   merged3 <- rbind(NA, merged3)
   merged3 <- rbind(NA, merged3)
   merged3 <- rbind(NA, merged3)
   
   #file_path <- file.path(dir, "snipper.xlsx") # revised to add tmp dir
   #openxlsx::write.xlsx(merged3, file_path)
   return(merged3)
   
}

### =======================FILE CONVERSION (to CSV)================###
# Issue (06 August 2025): Pop added at the end of the file
vcftocsv <- function(vcf, ref = NULL, dir = tempdir()) { # set dir to tmp
   library(dplyr)
   
   # Parse .vcf input
   if (tools::file_ext(vcf) == "vcf") {
      vcf_object <- vcfR::read.vcfR(vcf, verbose = FALSE)
      genotypes <- vcfR::extract.gt(vcf_object, return.alleles = TRUE)
      columns <- as.data.frame(vcfR::getFIX(vcf_object))
      ID <- columns$ID
      raw_df <- data.frame(ID, genotypes)
      
      final_df <- data.frame(t(raw_df)) %>%
         janitor::row_to_names(row_number = 1) %>%
         tibble::rownames_to_column(var = "Sample")

      
   } else if (tools::file_ext(vcf) == ".gz") {
      directory <- paste0(dir, "./vcf_files")
      dir.create(directory, recursive = TRUE, showWarnings = FALSE)
      utils::untar(vcf, exdir = directory) # do i need to specify file.path?
      
      wb <- list.files(path = file.path(directory), pattern = ".vcf$", full.names = TRUE)
      dflist <- lapply(wb, function(i) {
         vcf_obj <- vcfR::read.vcfR(i, verbose = FALSE)
         genotypes <- vcfR::extract.gt(vcf_obj, return.alleles = TRUE)
         columns <- as.data.frame(vcfR::getFIX(vcf_obj))
         ID <- columns$ID
         temp_df <- data.frame(ID, genotypes)
         
         data.frame(t(temp_df)) %>%
            janitor::row_to_names(row_number = 1) %>%
            tibble::rownames_to_column(var = "Sample")
      })
      
      final_df <- dplyr::bind_rows(dflist, .id = "source", .fill = TRUE)

      
   } else {
      stop("Unsupported file type for VCF input.")
   }
   
   # Ensure reference is provided
   if (is.null(ref) || ref == "") stop("Reference should be provided")
   
   # Merge or assign population
   if (file.exists(ref)) {
      # If ref is a valid filepath, treat it as a reference file
      ref_data <- read.csv(ref)
      final_df <- left_join(final_df, ref_data, by = "Sample")
   } else if (is.character(ref)) {
      # If ref is a string, assign it directly to the Population column
      final_df$Population <- ref
   } else {
      stop("Invalid reference input: must be a file path or string.")
   }
   
   # solution to the pop and sample organization (06 August 2025)
   sample <- final_df[,1]
   nosample <- final_df[,-1]
   data_cols <- ncol(nosample) - 1 #number of data cols
   total_cols <- as.integer(ncol(final_df)) #number of total cols
   data_gt <- final_df[, 2:data_cols] #subset data
   pop <- final_df[,total_cols] # subset pop
   final_df <- bind_cols(sample, pop, data_gt)
   final_df <- final_df %>%
      rename(Sample = 1, Population = 2)
   
   return(final_df)
}

### =======================FILE CONVERSION (UAS2CSV)===============###
uas2csv <- function(files = files, population = pop_file, reference = FALSE, dir = dir){
   
   if(!require("pacman")) {
      install.packages("pacman")
   }
   
   pacman::p_load(tools, utils, readxl, stats, dplyr, purrr, tidyr, readr, install = TRUE)
   
   if(!file.exists(files)){
      stop("File does not exist in the working directory")
   } else {
      
      if(tools::file_ext(files) == "zip"){
         utils::unzip(files, 
                      files = NULL, 
                      list = FALSE, 
                      overwrite = TRUE, 
                      exdir = file.path(dir, "UAS_genotypes"))
         
         data_path <- file.path(dir, "UAS_genotypes")
         data.files <- list.files(path = data_path, full.names = TRUE)
         data.list <- list.files(path = data_path, pattern = ".xlsx", full.names = TRUE)
         all.list <- list()
         
         for (i in data.list) {
            all.list[[i]] = readxl::read_excel(i, sheet = 1, col_names = TRUE, row.names(data.files))}
         
         # correct names of dataframes
         new_colnames <- c("Sample", "ID", "Allele")
         dflist_new <- lapply(all.list, setNames, new_colnames)
         
         dflist_corrected <- lapply(
            dflist_new, 
            function(x){ 
               stats::aggregate(Allele ~ Sample + ID, x, paste, collapse = "/")
            })
         
         # pivot wider
         df_list <- purrr::map(dflist_corrected, ~(tidyr::pivot_wider(.x, names_from= Sample, values_from = Allele)))
         
         library(dplyr) # usually needs to be specified
         library(purrr)
         
         merged <- df_list %>% purrr::reduce(full_join, by= "ID")
         id <- merged$ID 
         
         #transpose
         correct_alleles <- merged
         correct_alleles <- data.frame(t(correct_alleles))
         names(correct_alleles) <- correct_alleles[1,] 
         correct_alleles <- correct_alleles[-1,] 
         
         samples <- data.frame(colnames(merged[,-1]))
         corrected <- correct_alleles
         
         Samples <- rownames(corrected)
         corrected <- data.frame(Samples, corrected)
         
         
      } else if(tools::file_ext(files) == "tar"){ 
         
         untar(files, files = NULL, list = FALSE, exdir = file.path(dir, "UAS_genotypes"))
         
         new <- list.files(path = file.path(dir, "UAS_genotypes"), recursive = TRUE, full.names = TRUE)
         file.copy(from = new, to = file.path(dir, "UAS_genotypes"), overwrite = TRUE)
         
         data_path <- file.path(dir, "UAS_genotypes")
         data.files <- list.files(path = data_path, full.names = TRUE)
         data.list <- list.files(path = data_path, pattern = ".xlsx$", full.names = TRUE)
         all.list <- list()
         
         for (i in data.list) {
            all.list[[i]] = readxl::read_excel(i, sheet = 1, col_names = TRUE, row.names(data.files))}
         
         # correct names of dataframes
         new_colnames <- c("Sample", "ID", "Allele")
         dflist_new <- lapply(all.list, setNames, new_colnames)
         
         dflist_corrected <- lapply(
            dflist_new, 
            function(x){ 
               stats::aggregate(Allele ~ Sample + ID, x, paste, collapse = "/")
            })
         
         # pivot wider
         df_list <- purrr::map(dflist_corrected, ~(tidyr::pivot_wider(.x, names_from= Sample, values_from = Allele)))
         
         library(dplyr) # usually needs to be specified
         library(purrr)
         
         merged <- df_list %>% purrr::reduce(full_join, by= "ID")
         id <- merged$ID 
         
         #transpose
         correct_alleles <- merged
         correct_alleles <- data.frame(t(correct_alleles))
         names(correct_alleles) <- correct_alleles[1,] # removes letters are columns, but duplicates ID row+column
         correct_alleles <- correct_alleles[-1,] #will remove the duplicate ID row
         
         samples <- data.frame(colnames(merged[,-1]))
         corrected <- correct_alleles
         
         Samples <- rownames(corrected)
         corrected <- data.frame(Samples, corrected)
         
      } else {
         stop("Not a zipped file. Accepted are zipped and tar files")
      }
      
      if(reference == FALSE){
         
         #readr::write_csv(corrected, file = "01_merged_typed_data.csv", path = dir)
         return(corrected)
      } else {
         
         if(tools::file_ext(population) == "csv"){
            pop_data <- readr::read_csv(population)
         } else if(tools::file_ext(population) == "xlsx"){
            pop_data <- readxl::read_excel(population)
         } else {
            stop("File is not a csv or xlsx file")
         }
         
         matched <- corrected %>% dplyr::left_join(pop_data, by = "Sample")
         data_length <- as.integer(ncol(corrected) - 1) # count the number of markers
         data_matched <- matched[,2:data_length] # subset the markers
         meta_begin <- as.integer(ncol(corrected) + 1) # estimate the beginning of the metadata
         meta <- matched[,meta_begin:ncol(matched)] # subset the metadata 
         
         final_df <- dplyr::bind_cols(matched$Sample, meta, data_matched)
         names(final_df)[names(final_df) == "matched$Sample"] <- "Sample"
         
         return(final_df)
         #readr::write_csv(final_df, path = dir, file = "01_merged_typed_data.csv")
      }
   }
} 

### ======================= MARKER EXTRACTION =======================###
extract_markers <- function(input.file, 
                            snps.list = NULL, 
                            pos.list = NULL, 
                            bed.file = NULL, 
                            bim.file = NULL, 
                            fam.file = NULL, 
                            plink_args = NULL,
                            output.dir = output.dir, 
                            merged.file = "final_merged.vcf",
                            plink_path = plink_path) {
   
   file_type <- detect_file_type(input.file, bed.file, bim.file, fam.file)
   
   commands <- extraction(file_type, 
                          input.file, 
                          snps.list, 
                          pos.list, 
                          bed.file, 
                          bim.file, 
                          fam.file, 
                          plink_args,
                          output.dir, 
                          merged.file,
                          plink_path)
   
   return(list(
      extracted_files = list.files(output.dir, full.names = FALSE)
   ))
}

# TO DO #1: ADDITIONAL ARGUMENTS FOR PLINK - OPTIONAL FILTERING PROCEDURES
detect_file_type <- function(input.file, bed.file = NULL, bim.file = NULL, fam.file = NULL) {
   file_extension <- tools::file_ext(input.file)
   
   if (!is.null(bed.file) && !is.null(bim.file) && !is.null(fam.file)) {
      return("PLINK")
   } else if (grepl("\\.vcf\\.gz$", input.file)) {
      return("VCF_GZ")
   } else if (file_extension == "vcf") {
      return("VCF")
   } else if (file_extension == "bcf") {
      return("BCF")
   } else {
      stop("Unsupported file type. Please provide a VCF, VCF.GZ, BCF, or PLINK file.")
   }
}

merge_vcf_files <- function(output.dir, merged.file) {
   vcf_files <- list.files(output.dir, pattern = "*.vcf$", full.names = TRUE)
   
   if (length(vcf_files) == 0) {
      stop("No extracted VCF files found for merging.")
   }
   
   # July 30 note: replaced str_c with system2 to source out bcftools
   merge_command <- system2("bcftools",
                            args = c("concat", "-o", merged.file),
                            stdout = paste(vcf_files, collapse = " "))
   
   system(merge_command)
   print(paste("Merged VCF file created:", merged.file))
}

extraction <- function(file_type, 
                       input.file, 
                       snps.list = NULL, 
                       pos.list = NULL, 
                       bed.file = NULL, 
                       bim.file = NULL, 
                       fam.file = NULL, 
                       plink_args = NULL,
                       output.dir, 
                       merged.file,
                       plink_path) {
   
   if (!dir.exists(output.dir)) {
      dir.create(output.dir, recursive = TRUE)
   }
   
   if (!is.null(snps.list)) {
      args <- if (!is.null(plink_args)) paste(plink_args, collapse = " ") else ""
      
      # Extract markers using rsID
      command <- switch(file_type,
                        "VCF" = stringr::str_c(plink_path, 
                                               " --vcf ", 
                                               input.file, 
                                               " --const-fid 0 --cow --extract ", 
                                               snps.list, 
                                               " ", args,
                                               " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                               file.path(output.dir, "rsid_extracted")),
                        "VCF_GZ" = stringr::str_c(plink_path, 
                                                  " --vcf ", input.file, 
                                                  " --const-fid 0 --cow --extract ", snps.list, 
                                                  " ", args,
                                                  " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                                  file.path(output.dir, "rsid_extracted")),
                        "BCF" = stringr::str_c(plink_path, 
                                               " --bcf ", 
                                               input.file, 
                                               " --const-fid 0 --cow --extract ", 
                                               snps.list, 
                                               " ", args,
                                               " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                               file.path(output.dir, "rsid_extracted")),
                        "PLINK" = stringr::str_c(plink_path, 
                                                 " --bfile ", 
                                                 sprintf("--bed %s --bim %s --fam %s", bed.file, bim.file, fam.file), 
                                                 " --const-fid 0 --cow --extract ", 
                                                 snps.list, 
                                                 " ", args,
                                                 " --keep-allele-order --allow-no-sex --allow-extra-chr --recode vcf --out ", 
                                                 file.path(output.dir, "rsid_extracted"))
      )
      system(command)  # Execute PLINK
   } else if (!is.null(pos.list)) {
      args <- if (!is.null(plink_args)) paste(plink_args, collapse = " ") else ""
      
      # Extract markers using positions (creates multiple output files)
      command_list <- lapply(seq_len(nrow(pos.list)), function(i) {
         chr_num <- pos.list[i, 1]
         start_bp <- pos.list[i, 2]
         end_bp <- pos.list[i, 3]
         
         # Subtract 1 from start position for filename
         filename_base <- paste0("chr", chr_num, "_", start_bp - 1)
         output_file <- file.path(output.dir, paste0(filename_base, ".vcf"))
         
         cmd <- stringr::str_c(
            plink_path, " --", tolower(file_type), " ", input.file,
            " --cow --chr ", chr_num,
            " --from-bp ", start_bp,
            " --to-bp ", end_bp,
            " ", args,
            " --recode vcf --keep-allele-order --out ", tools::file_path_sans_ext(output_file)
         )
         
         system(cmd)  # Run immediately
         return(cmd)
      })
      
      merge_vcf_files(output.dir, merged.file)  # Merge extracted files
      
      return(command_list)
   } else {
      stop("Either `snps.list` or `pos.list` must be provided.")
   }
   
   return(command)
}


### ========================== CONCORDANCE =======================###

concordance <- function(file1, file2, haplotypes = FALSE){
   
   if (!require("pacman")){
      install.packages("pacman")
   }
   
   pacman::p_load(readr, readxl, tools, dplyr, janitor, purrr, tibble, tidyselect, QurvE, tidyr, ggplot2, forcats, install = TRUE)
   
   # check if file exists
   # read in the file
   if(!file.exists(file1)){
      stop("First file does not exist in the working directory")
   } else {
      if(tools::file_ext(file1) == "csv"){
         file1 <- readr::read_csv(file1)
      } else if(tools::file_ext(file1) == "xlsx"){
         file1 <- readxl::read_excel(file1)
      }
   }
   
   if(!file.exists(file2)){
      stop("Second file does not exist in the working directory")
   } else {
      if(tools::file_ext(file2) == "csv"){
         file2 <- readr::read_csv(file2)
      } else if(tools::file_ext(file2) == "xlsx"){
         file2 <- readxl::read_excel(file2)
      }
   }
   
   file1 <- file1 %>%
      rename(Ind = 1)
   file2 <- file2 %>%
      rename(Ind = 1)
   
   file_list <- list(file1, file2)
   
   # check intersecting values
   # samples are in the first column
   overlaps <- as.list(intersect(file1$Ind, file2$Ind))
   
   
   #subset
   file_list2 <- lapply(
      file_list,
      function(x){
         x[x$Ind %in% overlaps, ]
      }
   )
   
   # transpose
   file_list3 <- lapply(
      file_list2,
      function(x){
         library(dplyr)
         data.frame(t(x)) %>%
            janitor::row_to_names(row_number = 1) %>%
            tibble::rownames_to_column(., var = "markers") 
      }
   )
   
   # rearrange column order
   overlaps <- as.character(overlaps)
   markers1 = file_list3[[1]]$markers
   markers2 = file_list3[[2]]$markers
   
   file_list4 <- lapply(
      file_list3,
      function(x){
         #markers = x$markers
         relocate(x, any_of(overlaps)) 
      }
   )
   
   
   #re-add marker column
   file_list4[[1]]$markers = markers1
   file_list4[[2]]$markers = markers2
   
   
   # merge 
   merged <- file_list4 %>% purrr::reduce(full_join, by= "markers")
   ID <- merged$markers
   
   merged <- lapply(
      merged, 
      function(x){
         gsub(pattern = "|", replacement = "/", x = x, fixed = TRUE)}
   )
   
   merged <- as.data.frame(merged)
   
   if(haplotypes == TRUE){
      
      print("Assuming the data are haplotypes.")
      
   } else if(haplotypes == FALSE){
      merged <- merged %>% mutate(across(tidyselect::everything(), ~ case_when(
         . == "A" ~ "A/A",
         . == "T" ~ "T/T",
         . == "C" ~ "C/C",
         . == "G" ~ "G/G",
         . == "T/C" ~ "C/T",
         . == "T/G" ~ "G/T",
         . == "T/A" ~ "A/T",
         . == "G/A" ~ "A/G",
         . == "C/G" ~ "G/C",
         . == "C/A" ~ "A/C",
         TRUE ~ .x)))
   } else {
      stop("Parameter haplotype is required.")
   }
   
   
   # subset based on those with .x and .y
   # rearrange based on order
   overlap1 <- merged %>% select(dplyr::ends_with(".x"))
   overlap2 <- merged %>% select(dplyr::ends_with(".y"))
   
   for_conc <- QurvE::zipFastener(overlap1, overlap2, along = 2)
   for_conc2 <- data.frame(ID, for_conc)
   for_conc2[is.na(for_conc2)] <- "N"
   names(for_conc2) <-  sub('^X', '', names(for_conc2))
   
   concordance <- bind_cols(for_conc2 %>%
                               tidyr::gather(var, val, -matches("(.x$|ID)")) %>%
                               select(ID,val), for_conc2 %>%
                               tidyr::gather(var2, val2, -matches("(*.y$|ID)")) %>%
                               select(val2)) %>%
      add_count(ID) %>%
      group_by(ID) %>%
      summarise(
         Total = paste((ncol(for_conc2) - 1)/2),
         Incomparable = paste(sum(val == "N") + sum(val2 == "N")),
         Concordant = paste(sum(val == val2)),
         Discordant = paste(((ncol(for_conc2) - 1)/2)- sum(val == val2) - (sum(val == "N") + sum(val2 == "N")))
      ) %>%
      left_join(for_conc2, by = c("ID" = "ID")) 
   
   
   readr::write_csv(concordance, file = "concordance.csv")
   
   pivot <- concordance[,1]
   pivot2 <- concordance[,3:5]
   pivot <- data.frame(pivot, pivot2)
   pivot <- pivot %>%
      tidyr::pivot_longer(!ID,
                          names_to = 'Condition',
                          values_to = 'Count'
      )
   
   # select 
   Count <- as.integer(pivot$Count)
   rsID <- pivot$ID
   Condition <- pivot$Condition
   
   visual <- data.frame(rsID, Count, Condition)
   
   library(ggplot2)
   
   plot_conc <- visual %>%
      arrange(Count) %>%
      arrange(Condition) %>%
      mutate(rsID = forcats::fct_inorder(rsID)) %>%
      ggplot(aes(fill = Condition, x= rsID, y=Count)) +
      geom_bar(position = "stack", stat = "identity") +
      theme(
         axis.text.x = element_text(
            angle = 90,
            vjust = .3), 
         panel.background = element_blank()) +
      scale_fill_manual(values= c("#1ca7ec",
                                  "#fb7a8e",
                                  "#1f2f98"
      ))
   
   return(list(
      results = concordance,
      plot = plot_conc 
   ))   
}


### ================= POP STAT ====================###
load_input_file <- function(input) {
   if (tools::file_ext(input) == "csv") {
      return(readr::read_csv(input))
   } else if (tools::file_ext(input) == "xlsx") {
      return(readxl::read_excel(input))
   } else {
      stop("Input file should be in csv or xlsx format.")
   }
}


clean_input_data <- function(file) {
   library(dplyr)
   
   file <- lapply(file, function(x) gsub("|", "/", x, fixed = TRUE))
   file <- as.data.frame(file)
   file[is.na(file)] <- "N"
   
   file <- file %>%
      mutate(across(everything(), ~ case_when(. == "N/A" ~ "N", . == "NA" ~ "N", TRUE ~ .x))) %>%
      rename(Ind = 1, Pop = 2)
   
   file <- as.data.frame(file)
   
   return(file)
}


convert_to_genind <- function(file) {
   ind <- as.character(file$Ind)
   pop <- as.character(file$Pop)
   fsnps_geno <- file[, 3:ncol(file)]
   
   fsnps_gen <- adegenet::df2genind(fsnps_geno, 
                                    ind.names = ind, 
                                    pop = pop, 
                                    sep = "/", 
                                    NA.char = "N", 
                                    ploidy = 2, 
                                    type = "codom")
   
   devtools::source_url("https://raw.githubusercontent.com/Tom-Jenkins/utility_scripts/master/TJ_genind2genepop_function.R")
   fsnps_gen@pop <- as.factor(file$Pop)
   
   return(fsnps_gen)
}

# Need to correct, genind2hierfstat is outputting an error 
compute_population_stats <- function(fsnps_gen) {
   mar_matrix <- allelic.richness(genind2hierfstat(fsnps_gen))$Ar %>%
      apply(MARGIN = 2, FUN = mean) %>%
      round(digits = 3)
   mar_list <- as.list(mar_matrix)
   
   basic_fsnps <- hierfstat::basic.stats(fsnps_gen, diploid = TRUE)
   Ho <- apply(basic_fsnps$Ho, 2, mean, na.rm = TRUE) %>% round(2)
   He <- apply(basic_fsnps$Hs, 2, mean, na.rm = TRUE) %>% round(2)
   heterozygosity_df <- data.frame(Population = names(Ho), Ho, He) %>%
      tidyr::pivot_longer(cols = c("Ho", "He"), names_to = "Variable", values_to = "Value")
   
   # T-Test between Heterozygosities
   ttest_df <- data.frame(
      Locus = rownames(basic_fsnps$perloc),
      basic_fsnps$perloc,
      row.names = NULL
   )
   
   # Inbreeding Coefficient
   fis_values <- apply(basic_fsnps$Fis, 2, mean, na.rm = TRUE) %>%
      round(3)
   fis_df <- data.frame(Population = names(fis_values), Fis = fis_values)
   
   # Allele frequencies
   devtools::source_url("https://raw.githubusercontent.com/Tom-Jenkins/utility_scripts/master/TJ_genind2genepop_function.R")
   fsnps_gpop <- genind2genpop(fsnps_gen)
   allele_freqs <- t(makefreq(fsnps_gpop, quiet = FALSE, missing = NA)) %>%
      as.data.frame()
   
   return(list(
      mar_list = mar_list,
      heterozygosity = heterozygosity_df,
      ttest = ttest_df,
      inbreeding_coeff = fis_df,
      allele_frequencies = allele_freqs
   ))
}

compute_hardy_weinberg <- function(fsnps_gen) {
   # Hardy-Weinberg Equilibrium (List for export)
   fsnps_hwe <- as.numeric(round(pegas::hw.test(fsnps_gen, B = 1000), 6)) 
   
   # Chi-square test (Matrix for export, Data Frame for ggplot)
   fsnps_hwe_test <- data.frame(sapply(adegenet::seppop(fsnps_gen), 
                                       function(ls) pegas::hw.test(ls, B=0)[,3]))
   
   fsnps_hwe_chisq_matrix <- as.matrix(fsnps_hwe_test)  # Convert to matrix
   fsnps_hwe_chisq_df <- as.data.frame(t(fsnps_hwe_chisq_matrix)) %>% tibble::rownames_to_column("Population")
   
   return(list(
      hw_summary = fsnps_hwe,
      hw_matrix = fsnps_hwe_test,
      hw_dataframe = fsnps_hwe_chisq_df
   ))
}


compute_fst <- function(fsnps_gen) {
   # Pairwise Fst matrix using Weir & Cockerham 1984 method
   fst_matrix_raw <- hierfstat::genet.dist(fsnps_gen, method = "WC84") %>%
      round(3)
   
   # Safeguard against empty results
   fst_list <- if (length(fst_matrix_raw) == 0) {
      list(message = "No Fst values calculated")
   } else {
      as.list(as.matrix(fst_matrix_raw))
   }
   
   # Tidy format for ggplot or tabular view
   fst_df <- as.data.frame(as.matrix(fst_matrix_raw)) %>%
      tibble::rownames_to_column("Site1") %>%
      tidyr::pivot_longer(cols = -Site1, names_to = "Site2", values_to = "Fst")
   
   return(list(
      fst_matrix = fst_list,   
      fst_dataframe = fst_df   
   ))
}


plot_heterozygosity <- function(Het_fsnps_df, out_dir) {
   library(ggplot2)
   out_path <- file.path(out_dir, "heterozygosity_plot.png")
   
   Het_fsnps_df$Variable <- as.factor(Het_fsnps_df$Variable)
   
   p <- ggplot(Het_fsnps_df, aes(x = Population, y = Value, fill = Variable)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.6), colour = "black") +
      scale_y_continuous(expand = c(0, 0), limits = c(0, 0.5)) +
      scale_fill_manual(
         values = c("royalblue", "#bdbdbd"),
         labels = c(expression(italic("H")[o]), expression(italic("H")[e]))
      ) +
      labs(y = "Heterozygosity") +
      theme(axis.text.x = element_text(size = 10, angle = 90, vjust = 0.5, face = "bold"))
   
   #ggsave(out_path, plot = p, width = 9, dpi = 300)
   return(p)
}

plot_fst_heatmap_interactive <- function(fst_df) {
   library(ggplot2)
   p <- ggplot(fst_df, aes(x = Site1, y = Site2, fill = Fst, label = round(Fst, 3))) +
      geom_tile(color = "black") +
      geom_text(aes(label = round(Fst, 3)), size = 3, color = "black") +
      scale_fill_gradient2(
         low = "blue", mid = "pink", high = "red",
         midpoint = max(fst_df$Fst, na.rm = TRUE) / 2
      ) +
      labs(x = "Site 1", y = "Site 2", fill = "Fst") +
      theme_minimal(base_size = 11) +
      theme(
         axis.text = element_text(face = "bold"),
         axis.text.x = element_text(angle = 45, hjust = 1)
      )
   
   #plotly::ggplotly(p, tooltip = c("x", "y", "label", "fill"))
   return(p)
}


export_results <- function(stats_matrix, hw_matrix, fst_matrix, dir = tempdir()) {
   timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
   out_file <- file.path(dir, paste0("population-statistics-results_", timestamp, ".xlsx"))
   
   datasets <- list(
      "Mean Allelic Richness"     = as.data.frame(stats_matrix$mar_list),
      "Heterozygosities"          = as.data.frame(stats_matrix$heterozygosity),
      "T-test per Locus"          = as.data.frame(stats_matrix$ttest),
      "Inbreeding Coefficient"    = as.data.frame(stats_matrix$inbreeding_coeff),
      "Allele Frequencies"        = as.data.frame(stats_matrix$allele_frequencies),
      "Hardy-Weinberg Equilibrium"= as.data.frame(hw_matrix$hw_summary),
      "Chi-square Test Results"   = as.data.frame(hw_matrix$hwe_matrix),
      "Pairwise Fst Matrix"       = as.data.frame(fst_matrix$fst_matrix)
   )
   
   #openxlsx::write.xlsx(datasets, file = out_file)
   return(out_file)
}



### ====================PCA======================###

compute_pca <- function(fsnps_gen) {
   library(ade4)
   library(adegenet)
   library(stats)
   x <- tab(fsnps_gen, NA.method = "mean")
   set.seed(9999)
   pca1 <- dudi.pca(x, scannf = FALSE, scale = FALSE, nf = min(7, ncol(x)))
   
   percent <- pca1$eig / sum(pca1$eig) * 100
   ind_coords <- as.data.frame(pca1$li)  # Convert matrix to data frame
   colnames(ind_coords) <- paste0("PC", seq_len(ncol(ind_coords)))
   
   ind_coords$Ind <- indNames(fsnps_gen)
   ind_coords$Site <- fsnps_gen@pop
   
   centroid <- aggregate(ind_coords[, -c(ncol(ind_coords), ncol(ind_coords)-1)], 
                         by = list(ind_coords$Site), 
                         FUN = mean)
   colnames(centroid)[1] <- "Site"
   centroid <- as.data.frame(centroid)  # Convert matrix to data frame
   
   return(list(pca1 = pca1, percent = percent, ind_coords = ind_coords, centroid = centroid))
}


get_colors_labels <- function(fsnps_gen, use_default, input_labels = NULL, input_colors = NULL) {
   library(ade4)
   library(adegenet)
   if (use_default) {
      labels <- levels(as.factor(fsnps_gen@pop))
      colors <- RColorBrewer::brewer.pal(n = length(labels), name = "Set1")
   } else {
      labels <- input_labels
      colors <- input_colors
      
      if (is.null(labels) || is.null(colors)) {
         stop("Custom labels and colors must be provided if use_default is FALSE.")
      }
      
      if (length(labels) != length(colors)) {
         stop("Number of labels and colors must match.")
      }
   }
   
   return(list(labels = labels, colors = colors))
}

plot_pca <- function(ind_coords, centroid, percent, labels_colors, filename, width = 8, height = 8, pc_x = 1, pc_y = 2) {
   library(ggplot2)
   library(stats)
   # Convert matrices to data frames if needed
   if (!is.data.frame(ind_coords)) ind_coords <- as.data.frame(ind_coords)
   if (!is.data.frame(centroid)) centroid <- as.data.frame(centroid)
   
   if (!is.list(labels_colors)) {
      stop("labels_colors must be a named list containing 'labels' and 'colors'.")
   }
   
   # Ensure Site column exists and matches colors
   ind_coords$Site <- factor(ind_coords$Site, levels = labels_colors$labels)
   colors_named <- setNames(labels_colors$colors, labels_colors$labels)
   
   # Labels for axes
   xlab <- paste("PC", pc_x, " (", format(round(percent[pc_x], 1), nsmall = 1), "%)", sep = "")
   ylab <- paste("PC", pc_y, " (", format(round(percent[pc_y], 1), nsmall = 1), "%)", sep = "")
   
   # Create plot
   plot <- ggplot(data = ind_coords, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y))) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(aes(fill = Site), shape = 21, size = 4, show.legend = FALSE) +
      geom_label_repel(data = centroid, aes_string(x = paste0("PC", pc_x), y = paste0("PC", pc_y), label = "Site"), 
                       size = 4, show.legend = FALSE, max.overlaps = Inf) +
      scale_fill_manual(values = colors_named) +
      scale_colour_manual(values = colors_named) +
      labs(x = xlab, y = ylab)
   
   # Save plot
   ggsave(filename = filename, plot = plot, width = width, height = height, dpi = 600)
}


explore.pca <- function(input, default.colors.labels = TRUE, pca.labels = NULL, color.palette = NULL, set.size = FALSE, width = NULL, height = NULL, add.pc = FALSE, add.pc.x = NULL, add.pc.y = NULL) {
   file <- load_input_file(input)
   file <- clean_input_data(file)
   
   fsnps_gen <- convert_to_genind(file)
   
   pca_results <- compute_pca(fsnps_gen)
   
   labels_colors <- get_colors_labels(fsnps_gen, default.colors.labels, pca.labels, color.palette)
   
   if (set.size) {
      if (is.null(width) || is.null(height)) {
         stop("If set.size is TRUE, both width and height must be provided.")
      }
   } else {
      width <- 8
      height <- 8
   }
   
   plot_pca(pca_results$ind_coords, pca_results$centroid, pca_results$percent, labels_colors, filename = "pca.png", width = width, height = height, pc_x = 1, pc_y = 2)
   
   if (add.pc) {
      plot_pca(pca_results$ind_coords, pca_results$centroid, pca_results$percent, labels_colors, filename = "pca2.png", width = width, height = height, pc_x = add.pc.x, pc_y = add.pc.y)
   }
}

### ===================== STRUCTURE ANALYSIS ===============###
run_structure <- function(
      file,
      k.range = 1:5,
      num.k.rep = 3,
      burnin = 1000,
      numreps = 1000,
      noadmix = FALSE,
      phased = FALSE,
      ploidy = 2,
      linkage = FALSE,
      structure_path = str_path,
      dir = dir_input #, # changed to structure_files
      #plot_dir = file.path(output_base_dir, "evanno_plots"),
      #clumpp_plots = file.path(output_base_dir, "str_plots")
      #plot.out = TRUE,
      #delete.files = TRUE
){
   dir.create(dir)
   #dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
   #dir.create(clumpp_plots, recursive = TRUE, showWarnings = FALSE)
   
   ext <- tools::file_ext(file)
   
   raw <- if (ext == "csv") {
      readr::read_csv(file)
   } else if (ext == "xlsx") {
      readxl::read_excel(file)
   } else {
      stop("Unsupported file type. Please provide a CSV or XLSX.")
   }
   
   library(dplyr) # to make sure it works
   
   # Homogenize file format
   # check allele format
   raw <- as.data.frame(lapply(raw, function(x) gsub("|", "/", x, fixed = TRUE)))
   raw[is.na(raw)] <- "N"
   raw <- dplyr::mutate(raw, across(everything(), ~ dplyr::case_when(. == "N/A" ~ "N", . == "NA" ~ "N", TRUE ~ .x))) |> 
      dplyr::rename(Ind = 1, Pop = 2)
   
   # For Plotting
   populations_df <- raw[,1:2]
   colnames(populations_df) <- c("Label", "Population")
   populations_df$Label <- rownames(populations_df)
   
   ### Change pop to numeric - for STRUCTURE
   pop_df <- as.data.frame(raw$Pop) %>%
      dplyr::rename(pops = 1)
   
   # Get total no. of pops
   pop_df_unique <- as.data.frame(pop_df[!duplicated(pop_df), ]) %>%
      dplyr::rename(pops = 1)
   # add row names as numbers
   pop_df_unique$num <- rownames(pop_df_unique)
   write.csv(pop_df_unique, file = "population_order.csv")
   
   # replace the pops in the original df (pop_df) with the numbers
   pop_df_corr <- left_join(pop_df, pop_df_unique, by = "pops") 
   pops <- pop_df_corr$num
   
   ### Change Ind to numeric
   ind_only <- as.data.frame(raw$Ind)
   ind_only$num <- rownames(ind_only)
   
   ind <- as.character(ind_only$num)
   pop <- as.character(pops)
   geno <- raw[, 3:ncol(raw)]
   
   genind_obj <- adegenet::df2genind(geno, ind.names = ind, pop = pop, sep = "/", NA.char = "N", ploidy = 2, type = "codom")
   genind_obj@pop <- as.factor(pop)
   
   str.file <- to_structure(genind_obj, file = "structure_input.str", dir = dir)
   
   result <- running_structure(
      input_file = str.file,
      k.range = k.range,
      num.k.rep = num.k.rep,
      burnin = burnin,
      numreps = numreps,
      noadmix = noadmix,
      output_dir = dir,
      structure_path = structure_path
   )
   
   # download each item on the list of the q matrices
   ############ NEED TO VERIFY NAME
   
   qmatrices_list <- q_matrices()
   for(matrices in qmatrices_list){
      con <- file.path(paste0(dir, names(matrices), "matrices"))
      writeLines(as.character(item), con = dir)
   }
   
   #devtools::install_github("sa-lee/starmie")
   str.dir <- dir
   str.files <- list.files(str.dir, pattern = "\\_f$", full.names = TRUE)
   str.data <- lapply(str.files, starmie::loadStructure)
   
   for (i in str.data){
      #path_plot <- paste(clumpp_plots, i$K)
      file_name <- paste0(dir, "/", i$K, "plot.png")
      plotQ(i, populations_df, outfile = file_name)
   }
   
}


## Adapted from the dartR package
to_structure <- function(genind_obj, file = "structure_input.str", include_pop = TRUE, dir = output_base_dir) {
   out_path <- file.path(dir, file)
   # Get basic info
   ind <- adegenet::indNames(genind_obj)
   pop <- if (include_pop) as.character(genind_obj@pop) else rep(1, length(ind))
   ploidy <- max(genind_obj@ploidy)
   loci <- adegenet::locNames(genind_obj)
   n_loc <- length(loci)
   
   # Convert allele matrix to integer format
   allele_matrix <- matrix(-9, nrow = length(ind), ncol = n_loc * ploidy)
   colnames(allele_matrix) <- paste(rep(loci, each = ploidy), paste0(".a", 1:ploidy), sep = "")
   
   for (i in seq_along(ind)) {
      for (j in seq_along(loci)) {
         allele_values <- genind_obj@tab[i, grep(paste0("^", loci[j], "\\."), colnames(genind_obj@tab))]
         allele_values[is.na(allele_values)] <- -9
         allele_matrix[i, ((j - 1) * ploidy + 1):(j * ploidy)] <- allele_values
      }
   }
   
   final_data <- data.frame(ID = ind, POP = pop, allele_matrix, stringsAsFactors = FALSE)
   
   
   write.table(final_data, file = out_path, quote = FALSE, sep = " ", row.names = FALSE, col.names = FALSE)
   return(out_path)
}

## Directly from the dartR package, revised sections with comments
.structureParseQmat <- function (q.mat.txt, 
                                 pops) {
   
   q.mat.txt <- sub("[*]+", "", q.mat.txt)
   q.mat.txt <- sub("[(]", "", q.mat.txt)
   q.mat.txt <- sub("[)]", "", q.mat.txt)
   q.mat.txt <- sub("[|][ ]+$", "", q.mat.txt)
   cols1to4 <- c("row", "id", "pct.miss", "orig.pop")
   
   strsplit(q.mat.txt, " ") %>% 
      purrr::map(function(q) {
         q <- q[!q %in% c("", " ", ":")] %>% 
            as.character() %>% 
            rbind() %>% 
            as.data.frame(stringsAsFactors = FALSE)
         
         stats::setNames(q, 
                         c(cols1to4,
                           paste("Group", 1:(ncol(q) - 4), sep = ".")))
         
      }) %>% 
      dplyr::bind_rows() %>% 
      dplyr::mutate_at(dplyr::vars("row",
                                   "pct.miss",
                                   "orig.pop",
                                   dplyr::starts_with("Group.")), as.numeric) %>% 
      dplyr::mutate(orig.pop = if (!is.null(pops)) {
         pops[.data$orig.pop]
      } else {
         .data$orig.pop
      })
}


## Directly from the dartR package, revised sections with comments
structureRead <- function(file,
                          pops = NULL) {
   
   if (!file.exists(file)) {
      stop(error(paste("the file '", file, "' can't be found.", 
                       sep = "")))
   }
   
   result <- scan(file, "character", quiet = TRUE)
   loc <- grep("Estimated", result, ignore.case = FALSE, 
               value = FALSE)
   est.ln.prob <- as.numeric(result[loc[1] + 6])
   loc <- grep("likelihood", result, ignore.case = FALSE, 
               value = FALSE)
   mean.lnL <- as.numeric(result[loc[1] + 2])
   var.lnL <- as.numeric(result[loc[2] + 2])
   loc <- grep("MAXPOPS", result, value = F)
   maxpops <- result[loc]
   maxpops <- sub("MAXPOPS=", "", maxpops)
   maxpops <- as.integer(sub(",", "", maxpops))
   loc <- grep("GENSBACK", result, value = F)
   gensback <- result[loc]
   gensback <- sub("GENSBACK=", "", gensback)
   gensback <- as.integer(sub(",", "", gensback))
   smry <- c(k = maxpops, est.ln.prob = est.ln.prob, mean.lnL = mean.lnL, 
             var.lnL = var.lnL)
   result <- scan(file, "character", sep = "\n", 
                  quiet = TRUE)
   first <- grep("(%Miss)", result, value = FALSE) + 1
   last <- grep("Estimated Allele", result, value = FALSE) - 
      1
   tbl.txt <- result[first:last]
   tbl.txt <- sub("[*]+", "", tbl.txt)
   tbl.txt <- sub("[(]", "", tbl.txt)
   tbl.txt <- sub("[)]", "", tbl.txt)
   tbl.txt <- sub("[|][ ]+$", "", tbl.txt)
   prior.lines <- grep("[|]", tbl.txt)
   
   no.prior <- if (length(prior.lines) < length(tbl.txt)) {
      no.prior.q.txt <- if (length(prior.lines) == 0){ 
         tbl.txt
      }else{
         tbl.txt[-prior.lines]
      }
      .structureParseQmat(no.prior.q.txt, pops)
   } else { 
      NULL
   }
   
   # Initial section's error: no function to return from, jumping to top level
   # Resolved by wrapping in a function 'maxpops_check'
   maxpops_check <- function(maxpops, smry, no.prior){   
      if (maxpops == 1) {
         no.prior$row <- NULL
         return(list(summary = smry, q.mat = no.prior, prior.anc = NULL))
      }}
   
   has.prior <- if (length(prior.lines) > 0) {
      prior.txt <- strsplit(tbl.txt[prior.lines], "[|]")
      prior.q.txt <- unlist(lapply(prior.txt, function(x) x[1]))
      df <- .structureParseQmat(prior.q.txt, pops)
      
      prior.anc <- purrr::map(prior.txt, function(x) {
         anc.mat <- matrix(NA, nrow = maxpops, ncol = gensback + 1)
         rownames(anc.mat) <- paste("Pop", 1:nrow(anc.mat),  sep = ".")
         colnames(anc.mat) <- paste("Gen", 0:gensback, sep = ".")
         
         x <- sapply(strsplit(x[-1], "\\s|[:]"), function(y) {
            y <- y[y != ""]
            y[-1]
         })
         
         for (i in 1:ncol(x)) {
            pop <- as.numeric(x[1, i])
            anc.mat[pop, ] <- as.numeric(x[-1, i])
         }
         anc.mat
      }) %>% stats::setNames(df$id)
      prob.mat <- t(sapply(1:nrow(df), function(i) {
         pop.probs <- rowSums(prior.anc[[i]])
         pop.probs[is.na(pop.probs)] <- df$Group.1[i]
         pop.probs
      }))
      colnames(prob.mat) <- paste("Group", 1:ncol(prob.mat), 
                                  sep = ".")
      df$Group.1 <- NULL
      df <- cbind(df, prob.mat)
      list(df = df, prior.anc = prior.anc)
   }else{
      NULL
   }
   
   has.prior.df <- if(is.null(has.prior)){ 
      NULL
   }else{ 
      has.prior$df
   }
   
   q.mat <- rbind(no.prior, has.prior.df)
   q.mat <- q.mat[order(q.mat$id), ] # removed q.mat$row since row is nonexistent
   #q.mat$row <- NULL
   rownames(q.mat) <- NULL
   
   # Initial single-line code did not work, revised as ff::
   subset <- q.mat[,-(1:3)]
   subset <- as.data.frame(lapply(subset, as.numeric))
   normalized <- subset/rowSums(subset)
   q.mat <- cbind(q.mat[,1:3], normalized)
   
   #q.mat[, -c(1:3)] <- t(apply(q.mat[, -c(1:3)], 1, function(i) i/sum(i)))
   
   prior.anc <- if (is.null(has.prior)){ 
      NULL
   }else {
      has.prior$prior.anc
   }
   
   list(summary = smry, q.mat = q.mat, prior.anc = prior.anc)
}


## Directly from the dartR package, revised sections with comments
utils.structure.evanno <- function (sr, plot = TRUE) 
{
   if (!"structure.result" %in% class(sr)) {
      stop(error("'sr' is not a result from 'structure.run'."))
   }
   k.tbl <- table(sapply(sr, function(x) x$summary["k"]))
   
   
   #f (length(k.tbl) < 3) 
   # stop(error("must have at least two values of k."))
   sr.smry <- t(sapply(sr, function(x) x$summary))
   ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[,"k"], mean)
   sd.ln.k <- tapply(sr.smry[, "est.ln.prob"], sr.smry[,"k"], sd)
   
   ln.pk <- diff(ln.k)
   ln.ppk <- abs(diff(ln.pk))
   delta.k <- sapply(2:(length(ln.k) - 1), function(i) {
      abs(ln.k[i + 1] - (2 * ln.k[i]) + ln.k[i - 1])/sd.ln.k[i]
   })
   df <- data.frame(k = as.numeric(names(ln.k)), 
                    reps = as.numeric(table(sr.smry[, "k"])),
                    mean.ln.k = as.numeric(ln.k), 
                    sd.ln.k = as.numeric(sd.ln.k), 
                    ln.pk = c(NA, ln.pk), 
                    ln.ppk = c(NA, ln.ppk, NA), 
                    delta.k = c(NA, delta.k, NA))
   
   rownames(df) <- NULL
   df$sd.min <- df$mean.ln.k - df$sd.ln.k
   df$sd.max <- df$mean.ln.k + df$sd.ln.k
   plot.list <- list(mean.ln.k = 
                        ggplot2::ggplot(df, 
                                        ggplot2::aes_string(x = "k", 
                                                            y = "mean.ln.k")) + 
                        ggplot2::ylab("mean LnP(K)") + 
                        ggplot2::geom_segment(ggplot2::aes_string(x = "k", 
                                                                  xend = "k", 
                                                                  y = "sd.min", 
                                                                  yend = "sd.max")), 
                     ln.pk = ggplot2::ggplot(df[!is.na(df$ln.pk), ], 
                                             ggplot2::aes_string(x = "k", 
                                                                 y = "ln.pk")) +
                        ggplot2::ylab("LnP'(K)"), 
                     ln.ppk = ggplot2::ggplot(df[!is.na(df$ln.ppk), ],
                                              ggplot2::aes_string(x = "k", 
                                                                  y = "ln.ppk")) +
                        ggplot2::ylab("LnP''(K)"))
   
   if (!all(is.na(df$delta.k))) {
      plot.list$delta.k <- ggplot2::ggplot(df[!is.na(df$delta.k), 
      ], ggplot2::aes_string(x = "k", y = "delta.k")) + 
         ggplot2::ylab(expression(Delta(K)))
   }
   for (i in 1:length(plot.list)) {
      plot.list[[i]] <- plot.list[[i]] + ggplot2::geom_line() + 
         ggplot2::geom_point(fill = "white", shape = 21, 
                             size = 3) + ggplot2::xlim(c(1, max(df$k))) + 
         ggplot2::theme(axis.title.x = ggplot2::element_blank())
   }
   if (plot) {
      p <- plot.list %>% purrr::map(function(x) {
         ggplot2::ggplot_gtable(ggplot2::ggplot_build(x))
      })
      maxWidth <- do.call(grid::unit.pmax, purrr::map(p, function(x) x$widths[2:3]))
      for (i in 1:length(p)) p[[i]]$widths[2:3] <- maxWidth
      p$bottom <- "K"
      p$ncol <- 2
      do.call(gridExtra::grid.arrange, p)
   }
   df$sd.min <- df$sd.max <- NULL
   invisible(list(df = df, plots = plot.list))
}


### Adapted from the dartR package

running_structure <- function(input_file,
                              k.range,
                              num.k.rep,
                              burnin,
                              numreps,
                              noadmix = FALSE,
                              phased = FALSE,
                              ploidy = 2,
                              linkage = FALSE,
                              structure_path = structure_path,
                              output_dir = dir #,
                              #plot_dir = file.path(output_dir, "evanno_plots")
){
   
   
   # Validate K range
   if (length(k.range) < 2) stop("Provide at least two K values for Evanno analysis.")
   #if (!dir.exists(plot_dir)) dir.create(plot_dir, recursive = TRUE)
   
   numinds <- length(readLines(input_file))
   numloci <- (ncol(read.table(input_file, header = FALSE, sep = " ")) - 2)/2
   
   base_label <- tools::file_path_sans_ext(basename(input_file))
   
   rep.df <- expand.grid(rep = 1:num.k.rep, k = k.range)
   rep.df$run <- paste0(base_label, ".k", rep.df$k, ".r", rep.df$rep)
   
   
   out_files <- lapply(1:nrow(rep.df), function(run_label) {
      #k <- rep.df[i, "k"]
      #####################
      out_path <- file.path(output_dir, rep.df[run_label, "run"])
      
      mainparams <- paste0(output_dir, "/mainparams")
      extraparams <- paste0(output_dir, "/extraparams")
      
      #mainparams <- paste0(dir, "mainparams")
      #extraparams <- paste0(dir, "extraparams")
      
      if(is.null(ploidy)){
         ploidy = "2"
      } else {
         ploidy = ploidy
      }
      
      writeLines(paste("#define", c(
         paste("MAXPOPS", rep.df[run_label, "k"]),
         paste("BURNIN", burnin),
         paste("NUMREPS", numreps),
         paste("INFILE", input_file),
         #paste("OUTFILE", out_path),
         paste("NUMINDS", numinds),
         paste("NUMLOCI", numloci),
         paste("PLOIDY", ploidy),
         "MISSING -9", 
         "ONEROWPERIND 1",
         "LABEL 1", 
         "POPDATA 1",
         "POPFLAG 0", 
         "LOCDATA 0", 
         "PHENOTYPE 0",
         "EXTRACOLS 0", 
         "MARKERNAMES 0",
         "RECESSIVEALLELES 0",
         "MAPDISTANCES 0",
         paste("PHASED", ifelse(phased, 1, 0)),
         "MARKOVPHASE 0",
         "NOTAMBIGUOUS -999"
      )), con = mainparams)
      
      # change alpha value divide 1 with max K value #######################################################
      alpha = 1/max(k.range)
      
      writeLines(paste("#define", c(
         paste("NOADMIX", ifelse(noadmix, 1, 0)), 
         paste("LINKAGE", ifelse(linkage, 1, 0)),
         "USEPOPINFO 0",
         "FREQSCORR 1",
         "ONEFST 0",
         "INFERALPHA 1",
         "POPALPHAS 0",
         paste("ALPHA", alpha),
         "INFERLAMBDA 0",
         "POPSPECIFICLAMBDA 0",
         "LAMBDA 1.0",
         "FPRIORMEAN 0.01",
         "FPRIORSD 0.05",
         "UNIFPRIORALPHA 1",
         "ALPHAMAX 10.0",
         "ALPHAPRIORA 1.0",
         "ALPHAPRIORB 2.0",
         "LOG10RMIN -4.0",
         "LOG10RMAX 1.0",
         "LOG10RSTAT -2.0",
         "GENSBACK 2",
         "MIGRPRIOR 0.01",
         "PFROMPOPFLAGONLY 0",
         "LOCISPOP 1",
         "LOCPRIORINIT 1.0",
         "MAXLOCPRIOR 20.0",
         "PRINTNET 1",
         "PRINTLAMBDA 1",
         "PRINTQSUM 1",
         "SITEBYSITE 0",
         "PRINTQHAT 0",
         "UPDATEFREQ 100",
         "PRINTLIKES 0",
         "INTERMEDSAVE 0",
         "ECHODATA 1",
         "ANCESTDIST 0",
         "NUMBOXES 1000",
         "ANCESTPINT 0.90",
         "COMPUTEPROB 1",
         "ADMBURNIN 500",
         "ALPHAPROPSD 0.025",
         "STARTATPOPINFO 0",
         "RANDOMIZE 0",
         "SEED 2245",
         "METROFREQ 10",
         "REPORTHITRATE 0"
      )), con = extraparams)
      
      #cmd <- paste(structure_path,
      #             "-i", input_file,
      #             "-m", mainparams,
      #             "-e", extraparams,
      #             "-o", out_path)
      
      # for(i in 1:nrow(rep.df)){
      cmd <- paste(structure_path,
                   "-i", input_file,
                   "-K", rep.df[run_label, "k"],
                   "-m", mainparams,
                   "-e", extraparams,
                   "-o", out_path)
      
      
      message("Running STRUCTURE: ", cmd)
      log <- tryCatch(system(cmd, intern = TRUE), error = function(e) e$message)
      writeLines(log, paste0(out_path, "_log.txt"))
      
      final_out <- paste0(output_dir, "_f")
      if (!file.exists(final_out) && file.exists(output_dir)) final_out <- output_dir
      if (!file.exists(final_out)) {
         warning("Missing output file for run: ", run_label)
         return(NULL)
      }
      
      result <- tryCatch(structureRead(final_out), error = function(e) {
         warning("Failed to parse STRUCTURE output: ", run_label)
         return(NULL)
      })
      result$label <- run_label
      result
      
   })
   
   run.result <- Filter(Negate(is.null), out_files)
   message("Successful STRUCTURE runs: ", length(run.result))
   names(run.result) <- sapply(run.result, `[[`, "label")
   class(run.result) <- c("structure.result", class(run.result))
   
   if (length(run.result) < 3) stop("Evanno analysis needs at least three successful runs.")
   
   ev <- utils.structure.evanno(run.result, plot = FALSE)
   
   # 
   lapply(names(ev$plots), function(pname) {
      plot_obj <- ev$plots[[pname]]
      png_path <- file.path(output_dir, paste0("Evanno_", pname, ".png"))
      grDevices::png(filename = png_path, width = 1600, height = 1200, res = 200)
      print(plot_obj)
      grDevices::dev.off()
      message("Saved PNG plot: ", png_path)
   })
   
   #if (plot.out && "delta.k" %in% names(ev$plots)) {
   suppressMessages(print(ev$plots$delta.k))
   #}
   
   #if (delete.files) unlink(output_dir, recursive = TRUE)
   
   invisible(list(
      results = run.result,
      evanno = ev,
      plot.paths = list.files(output_dir, full.names = TRUE)
   ))
}

# Plot STRUCTURE
plotQ <- function(qmat, populations_df, outfile = outfile) {
   
   # Revised to be compatible with large list of matrices
   facet = FALSE
   K <- qmat$K
   Label <- seq_len(nrow(qmat$ancest_df))
   
   clusters <- qmat$ancest_df %>%
      as.data.frame() %>%
      select(contains("Cluster"))
   colnames(clusters) <- gsub(" ", ".", colnames(clusters))
   
   df <- data.frame(ID = as.data.frame(Label), clusters)
   
   if (is.null(populations_df)){
      #Generate a plot without any family information
      #Q_melt <- reshape2::melt(Q, variable.name="Cluster")
      Q_melt <- stats::reshape(
         df, 
         varying = list(names(clusters)),
         v.names = "Value",
         timevar = "Cluster",
         times = paste0("K", seq_len(K)),
         direction = "long"
      )
      
      colnames(Q_melt) <- c("Label", "Cluster", "Value")
   } else{
      if (length(populations_df$Label) != length(rownames(df)))
         stop("Unequal sample size.")
      colnames(populations_df) <- c("Label", "Population")
      Q_merge <- merge(populations_df, df, by="Label")
      #Q_melt <- melt(Q_merge, id.vars=c("Label", "Population"), variable.name="Cluster")
      
      Q_melt <- stats::reshape(
         Q_merge, 
         varying = list(names(clusters)),
         v.names = "Value",
         timevar = "Cluster",
         times = paste0("K", seq_len(K)),
         direction = "long"
      )
   }
   
   #Generate plot
   # Updated to proper ggplot syntax
   Q_melt <- Q_melt[order(Q_melt$Cluster),]
   Q_melt$Label <- factor(Q_melt$Label)
   
   gg <- ggplot(Q_melt, aes(x=Label, y=Value, fill=Cluster))
   if (!is.null(populations_df)){
      if (facet){
         gg <- gg + facet_grid( Cluster ~ Population, scales = "free_x", space = "free_x")
      } else{
         gg <- gg + facet_grid( . ~ Population, scales = "free_x", space = "free_x")
      }
   } else{
      if (facet){
         gg <- gg + facet_grid( Cluster ~ ., scales = "free_x", space = "free_x")
      }
   }
   gg <- gg + geom_bar(stat = "identity", width=1)
   gg <- gg + scale_y_continuous(expand=c(0,0), breaks=c(0.25,0.75))
   gg <- gg + coord_cartesian(ylim=c(0,1))
   gg <- gg + xlab("Sample ID") + ylab("Proportion of cluster") 
   gg <- gg + theme_bw()
   gg <- gg + guides(fill=guide_legend(title="Cluster"))
   gg <- gg + theme(axis.text.x = element_text(angle = 90))
   
   ggplot2::ggsave(outfile, plot = gg, width = 12, height = 10, dpi = 600)
}

q_matrices <- function(dir = dir){
   output <- list.files(path = dir, pattern = "\\_f$", full.names = TRUE)
   output_list <- lapply(output, function(filepath){
      lines <- readLines(filepath)
      
      start <- grep("Inferred ancestry of individuals", lines) + 2
      end <- grep("^Estimated Allele Frequencies in each cluster", lines)[1] - 1
      
      qmatrices <- lines[start:end]
      
      qmatrices_data <- do.call(rbind, lapply(qmatrices, function(line){
         section <- unlist(strsplit(line, ":"))
         
         if(length(section) == 2) {
            props <- as.numeric(strsplit(trimws(section[2]), "\\s+")[[1]])
            return(props)
         } else {
            return(NULL)
         }
      }))
      return(qmatrices_data)
   })
   names(output_list) <- sapply(output_list, function(n) as.character(output$k))
   return(output_list) ########### DOUBLE CHECK NAMES
}


