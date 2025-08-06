library(shiny)
library(shinyjs)
library(ggplot2)
source("functions.R")

# load plink
plink_path <- Sys.which("plink")
if (plink_path == "") plink_path <- "/usr/local/bin/plink"

structure_path <- Sys.which("structure")
if (structure_path == "") plink_path <- "/usr/local/bin/structure"
