library(shiny)
library(shinyjs)
library(ggplot2)
library(bslib)
source("functions.R", local = TRUE)
useShinyjs()

# load plink
plink_path <- Sys.which("plink")
if (plink_path == "") plink_path <- "/usr/local/bin/plink"

#structure_path <- Sys.which("structure")
#if (structure_path == "") structure_path <- "/usr/local/bin/console/structure"
#structure_path <- "./structure.exe"