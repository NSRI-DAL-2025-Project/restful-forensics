# problem 1: docker logo not showing
# proble 2: copyright not showing properly, add at the bottom of the application
# problem 3: previously uploaded file showing upon reload, should be wiped clean
# proble 4: wont run in docker container
library(bslib)
library(shinyjs)
source("functions.R", local = TRUE)
#shiny::addResourcePath('www', '/srv/shiny-server/restful-forensics/www')
useShinyjs()

ui <- navbarPage(
   title = div(
      tags$img(src = "www/logo.png", height = "30px", style = "display: inline-block; vertical-align: center;"),
      tags$span("RESTful Forensics",
                style = "font-family: Carme, sans-serif; font-size: 26px; color: #92b2e4; vertical-align: middle; padding-left: 0px;") #,
      #tags$div(
      #   style = "position: fixed; bottom: 0, width: 100%; background-color: transparent; padding: 8px; text-align: center; font-size: 10px; color: #666;",
      #   HTML("&copy; 2025 DNA Analysis Laboratory, Natural Sciences Research Institute, University of the Philippines Diliman. All rights reserved.")
      #)
   ),
   
   tags$head(
      tags$link(rel = "stylesheet", href = "https://fonts.googleapis.com/css2?family=Carme&display=swap"),
      tags$style(HTML("
            body {font-family: 'Carme';}
            .navbar-nav > li:nth-child(1) > a { background-color: transparent !important; }
            .navbar-nav > li:nth-child(2) > a { background-color: transparent !important; }
            .navbar-nav > li:nth-child(3) > a { background-color: transparent !important; }
            .navbar-nav > li:nth-child(4) > a { background-color: #75a2bf !important; }
            .navbar-nav > li:nth-child(5) > a { background-color: #5e8cad !important; }
            .navbar-nav > li:nth-child(6) > a { background-color: #46769b !important; }
            .navbar-nav > li:nth-child(7) > a { background-color: #2f5f8a !important; }
            .navbar-nav > li:nth-child(8) > a { background-color: #174978 !important; }
            .navbar-nav > li:nth-child(9) > a { background-color: #003366 !important; }
      
            .clickable-card {
                border: 1px solid #ccc;
                border-radius: 6px;
                padding: 15px;
                margin-bottom: 10px;
                box-shadow: 2px 2px 6px rgba(0,0,0,0.1);
                cursor: pointer;
                background-color: #f9f9f9;
              }
              .card-header {
                font-weight: bold;
                font-size: 16px;
                color: #1c4e80;
              }
              .card-body {
                display: none;
                margin-top: 10px;
                font-size: 14px;
                color: #333;
              }
              
              .inner-card {
                background-color: #eef5fb;
                border: 1px solid #bdd3eb;
                border-radius: 6px;
                padding: 12px;
                margin-top: 10px;
                margin-bottom: 10px;
              }
              .inner-card h5 {
                margin-top: 0;
                font-size: 15px;
                color: #1c4e80;
              }
      
          "))
   ), # end of tags$head
   tags$script(HTML("
            $(document).on('click', '.clickable-card', function() {
              $(this).find('.card-body').slideToggle('fast');
            });
          ")),
   
   
   tabPanel(title = HTML("<span style = 'color:#   ;'>Homepage</span>"),
               p("This application is a compilation of the work on ancestry informative markers by the DNA Analysis Laboratory with an ongoing effort to expand to other marker types.")
                  
               
   ), # end of tab panel for homepage
   
   ## 1. Instructions Tab ----
   tabPanel(title = HTML("<span style = 'color:#ffffff;'>Instructions</span>"),
            fluidPage(
               
               div(class = "clickable-card",
                   div(class = "card-header", "🔄 File Conversion"),
                   div(class = "card-body",
                       p("Convert various files to commonly used input files."),
                       
                       div(class = "inner-card",
                           h5("A. Convert files to CSV and add population info"),
                           p("Input file: VCF, BCF, or PLINK (.bed, .bim, .fam) files."),
                           p("Expected output file: CSV file.")
                       ),
                       div(class = "inner-card",
                           h5("B. Convert ForenSeq UAS outputs to wide format"),
                           p("Input file: Compressed folder (.zip or .tar) of XLSX files."),
                           p("Expected output file: Single CSV file (merged XLSX files).")
                       ),
                       div(class = "inner-card",
                           h5("C. Convert CSV or XLSX files to a SNIPPER-compatible file"),
                           p("Input file: CSV or XLSX file."),
                           p("Expected output file: XLSX file.")
                       )
                       
                   )
               ),
               ## extraction
               div(class = "clickable-card",
                   div(class = "card-header", "🧬 SNP Extraction"),
                   div(class = "card-body",
                       p("Extract markers from sequenced data and perform concordance analysis."),
                       
                       div(class = "inner-card",
                           h5("A. Extract SNPs based on rsID or GRCh37/GRCh38 position"),
                           p("Input file:"),
                           p("(1) VCF, BCF, or PLINK (.bed, .bim, .fam) files."),
                           p("(2) Markers/position list — you may type rsIDs manually, upload a list, or use a POS txt file."),
                           p("The position list (txt file) should include:"),
                           tags$ul(
                              tags$li("[1] Chromosome number (integer)"),
                              tags$li("[2] Starting base-pair position"),
                              tags$li("[3] Final base-pair position")
                           )
                       ),
                       
                       div(class = "inner-card",
                           h5("B. Concordance analysis between files with the same samples"),
                           p("Input files: Two CSV or XLSX files."),
                           p("Expected output:"),
                           tags$ul(
                              tags$li("Concordance table"),
                              tags$li("Concordance plot")
                           )
                       )
                       
                   )
               ), # end of div for tab2
               
               ### POP Stat
               div(class = "clickable-card",
                   div(class = "card-header", "📝 Population Statistics"),
                   div(class = "card-body",
                       p("Calculate private alleles, heterozygosity, inbreeding coefficients, allele frequencies, and other basic populations statistics"),
                       p("Input file: CSV or XLSX file"),
                       p("Expected output files:"),
                       tags$ul(
                          tags$li("XLSX file with all table results"),
                          tags$li("Heterozygosity Plot"),
                          tags$li("Fst Plot")
                       )
                   )
               ), #end of div for popstat
               
               div(class = "clickable-card",
                   div(class = "card-header", "🔍 Exploratory Analysis"),
                   div(class = "card-body",
                       p("Run principal component analysis"),
                       p("Input file: CSV or XLSX file"),
                       p("Expected output file: PNG plots"))
               ), # end of div for pca
               
               div(class = "clickable-card",
                   div(class = "card-header", "📊 STRUCTURE Analysis"),
                   div(class = "card-body",
                       p("Generate STRUCTURE input files and pong compatible files. Visualize the possible results"),
                       p("Input file: CSV or XLSX file"),
                       p("Expected output file: Zipped files and PNG plots")))
            ) # end of fluidpage
            
   ), # end of tab panel 
   
   ## FILE CONVERSION
   tabPanel(
      title = HTML("<span style = 'color:#ffffff;'>File Conversion</span>"),
      tabsetPanel(
         
         # Subtab 1: Convert to CSV
         tabPanel("Convert files to CSV",
                  useShinyjs(),
                  sidebarLayout(
                     sidebarPanel(
                        radioButtons("inputType", "Choose Input File Type",
                                     choices = c("VCF file" = "vcf", "BCF file" = "bcf", "PLINK files (.bed/.bim/.fam)" = "plink")),
                        
                        conditionalPanel(
                           condition = "input.inputType == 'vcf'",
                           fileInput("vcfFile", "Upload VCF File")
                        ),
                        # added 4 Aug (missed)
                        conditionalPanel(
                           condition = "input.inputType == 'bcf'",
                           fileInput("bcfFile", "Upload BCF File")
                        ),
                        
                        conditionalPanel(
                           condition = "input.inputType == 'plink'",
                           fileInput("bedFile", "Upload BED File"),
                           fileInput("bimFile", "Upload BIM File"),
                           fileInput("famFile", "Upload FAM File")
                        ),
                        
                        radioButtons("poptype", "Do samples come from a single population?",
                                     choices = c("Yes" = "single", "No" = "multiplepop")),
                        
                        conditionalPanel(
                           condition = "input.poptype == 'multiplepop'",
                           fileInput("multiplepop", "Input reference file with sample ID and population"),
                           helpText("*Accepts XLSX and CSV files")
                        ),
                        
                        conditionalPanel(
                           condition = "input.poptype == 'single'",
                           textAreaInput("typePop", "Enter population", rows = 1)
                        ),
                        
                        actionButton("convertCSV", "Convert File to CSV", icon = icon("arrow-up-right-from-square"))
                        
                     ),
                     mainPanel(
                        tableOutput("previewTable"),
                        downloadButton("downloadConvertedCSV", "Download Converted CSV")
                     )
                  )
         ),
         
         tabPanel("Widen ForenSeq UAS files",
                  useShinyjs(),
                  sidebarLayout(
                     sidebarPanel(
                        fileInput("uas_zip", "Upload ZIP or TAR file",
                                  accept = c(".zip", ".tar")),
                        helpText("*Accepts compressed files containing XLSX files."),
                        fileInput("ref_file", "Optional Reference File (CSV or XLSX)",
                                  accept = c(".csv", ".xlsx")),
                        actionButton("run_uas2csv", "Run Conversion"),
                        br(), br(),
                        downloadButton("downloadUAScsv", "Download Converted CSV")
                     ),
                     mainPanel(
                        tableOutput("previewTableUAS"),
                        h4("Example Input XLSX Format"),
                        fluidRow(
                           column(6,
                                  h5("All alleles of available SNPs per sample are listed in a long format."),
                                  tableOutput("exampleXLSX")
                           )
                        )
                        #possibly add an option to view the first few lines of the result,
                     )
                  )
         ), #end of tabpanel
         
         
         tabPanel("Convert to SNIPPER-analysis ready file",
                  useShinyjs(),
                  sidebarLayout(
                     sidebarPanel(
                        fileInput("convertFile", "Upload File"),
                        helpText("*Accepts VCF, XLSX, and CSV files"),
                        fileInput("refFile", "Upload Reference File"),
                        helpText("*Accepts XLSX and CSV files"),
                        
                        checkboxInput("targetPop", "Subset Target Population?", value = FALSE),
                        textInput("targetPopName", "Target Population Name"),
                        actionButton("convertBtn", "Convert Format", icon = icon("arrow-up-right-from-square"))
                     ),
                     mainPanel(
                        fluidRow(
                           column(6,
                                  h5("Sample input file"),
                                  tableOutput("exampleTableSnipper")
                           ),
                           column(6,
                                  h5("Sample reference file"),
                                  tableOutput("exampleRefSnipper"))
                        ), # end of fluid row
                        tableOutput("previewTableSNIPPER"),
                        downloadButton("downloadConverted", "Download Converted File")
                     )
                  )
         )
         
      )
   ), # end of tabpanel
   
   ## MARKER EXTRACTION
   tabPanel(title = HTML("<span style = 'color:#ffffff;'>SNP Extraction</span>"),
            tabsetPanel(
               tabPanel("SNP Extraction",
                        useShinyjs(),
                        sidebarLayout(
                           sidebarPanel(
                              tabPanel("Marker Extraction",
                                       fluidPage(
                                          fileInput("markerFile", "Upload Genotype (VCF, BCF or PLINK) File"),
                                          
                                          radioButtons("markerType", "Choose Marker Type",
                                                       choices = c("rsid", "pos"), inline = TRUE),
                                          
                                          conditionalPanel(
                                             condition = "input.markerType == 'rsid'",
                                             radioButtons("rsidInputType", "RSID Input",
                                                          choices = c("manual", "upload")),
                                             conditionalPanel(
                                                condition = "input.rsidInputType == 'manual'",
                                                textAreaInput("typedRSIDs", "Enter RSIDs (one per line)", rows = 5)
                                             ),
                                             conditionalPanel(
                                                condition = "input.rsidInputType == 'upload'",
                                                fileInput("markerList1", "Upload RSID List File")
                                             )
                                          ),
                                          
                                          conditionalPanel(
                                             condition = "input.markerType == 'pos'",
                                             fileInput("markerList2", "Upload POS List (.csv, .xlsx)")
                                          ),
                                          
                                          fileInput("bedFile", "PLINK BED file (optional)"),
                                          fileInput("bimFile", "PLINK BIM file (optional)"),
                                          fileInput("famFile", "PLINK FAM file (optional)"),
                                          textAreaInput("plink_args",
                                                        label = "Additional PLINK Arguments",
                                                        placeholder = "--maf 0.05 --geno 0.1",
                                                        rows = 3,
                                                        width = "100%"),
                                          helpText("See https://www.cog-genomics.org/plink/ for options."),
                                          
                                          actionButton("extractBtn", "Run Marker Extraction", icon = icon("play")),
                                          
                                          downloadButton("downloadExtracted", "Download Merged VCF")
                                       )
                              )
                           ),
                           mainPanel(
                              h4("Example Input Formats"),
                              fluidRow(
                                 column(6,
                                        h5("rsID Format"),
                                        tableOutput("exampleRSID")
                                 ),
                                 column(6,
                                        h5("Position Format"),
                                        tableOutput("examplePOS")
                                 )
                              ),
                              hr(),
                              downloadButton("downloadExtracted", "Download Extracted Markers (VCF)"),
                              downloadButton("downloadConverted", "Download Extracter Markers (CSV)")
                           )
                        )
               ),
               tabPanel(HTML("<span style = 'color:#000000;'>Concordance Analysis</span>"),
                        useShinyjs(),
                        sidebarLayout(
                           sidebarPanel(
                              fileInput("concordanceFile1", "Upload File A"),
                              fileInput("concordanceFile2", "Upload File B"),
                              checkboxInput("isHaplotype", "Treat data as haplotypes", value = FALSE),
                              actionButton("compareBtn", "Run Concordance Analysis", icon = icon("play"))
                           ),
                           mainPanel(
                              h4("Example Input Formats"),
                              fluidRow(
                                 column(6,
                                        h5("File Format (for concordance)"),
                                        tableOutput("exampleTable")
                                 )
                              ),
                              hr(),
                              h4("Concordance Summary Table"),
                              tableOutput("concordanceResults"),
                              hr(),
                              h4("Concordance Plot"),
                              imageOutput("concordancePlot"),
                              hr(),
                              downloadButton("downloadConcordance", "Download Concordance Results"),
                              downloadButton("downloadConcordancePlot", "Download Plot")
                           )
                        )
               )
            )
   ), # end of tabpanel
   
   ## POP STAT
   tabPanel(HTML("<span style = 'color:#ffffff;'>Population Statistics</span>"),
            tabsetPanel(
               tabPanel("Perform Analysis",
                        useShinyjs(),
                        fileInput("popStatsFile", "Upload CSV or XLSX Dataset"),
                        
                        actionButton("runPopStats", "Analyze", icon = icon("magnifying-glass-chart")),
                        downloadButton("downloadStats", "Download Results (Excel)"),
                        
                        hr(),
                        h4("Example: Population File Format"),
                        tableOutput("examplePop")
               ), 
               tabPanel("1 Private Alleles",
                        h4("Private Alleles Summary"),
                        uiOutput("privateAllelePlot"),
               ),
               tabPanel("2 Heterozygosity",
                        h4("Observed vs Expected Heterozygosity"),
                        DT::dataTableOutput("heterozygosity_table"),
                        hr(),
                        h4("Heterozygosity Plot"),
                        imageOutput("heterozygosity_plot"),
                        downloadButton("downloadHeterozygosityPlot", "Download Plot")
               ),
               tabPanel("3 Inbreeding Coefficients",
                        h4("Inbreeding Coefficient by Population"),
                        DT::dataTableOutput("inbreeding_table")
               ),
               tabPanel("4 Allele Frequencies",
                        h4("Allele Frequency Table"),
                        DT::dataTableOutput("allele_freq_table")
               ),
               tabPanel("5 Hardy-Weinberg Equilibrium",
                        h4("HWE P-value Summary"),
                        uiOutput("hwe_summary"),
                        h4("Population-wise HWE Chi-Square Table"),
                        DT::dataTableOutput("hwe_chisq_table")
               ),
               tabPanel("6 Fst Values",
                        h4("Pairwise Fst Matrix"),
                        uiOutput("fstMatrixUI"),  
                        h4("Tidy Pairwise Fst Data"),
                        DT::dataTableOutput("fstDfTable"),
                        
                        hr(),
                        h4("Fst Heatmap"),
                        plotly::plotlyOutput("fst_heatmap_interactive", height = "600px")
                        #imageOutput("fst_heatmap"),
                        #downloadButton("downloadFstHeatmap", "Download Heatmap")
                        
               )
            )
   ), # end of tabpanel
   
   ## PCA
   tabPanel(HTML("<span style = 'color:#ffffff;'>Exploratory Analysis</span>"),
            sidebarLayout(
               sidebarPanel(
                  fileInput("pcaFile", "Upload SNP Data (in CSV or XLSX) for PCA", accept = c(".csv", ".txt")),
                  checkboxInput("useDefaultColors", "Use Default Colors and Labels", TRUE),
                  conditionalPanel(
                     condition = "!input.useDefaultColors",
                     fileInput("pcaLabels", "Upload PCA Labels"),
                     fileInput("colorPalette", "Upload Color Palette")
                  ),
                  numericInput("pcX", "PC Axis X", value = 1, min = 1),
                  numericInput("pcY", "PC Axis Y", value = 2, min = 1),
                  actionButton("runPCA", "Run PCA Analysis", icon = icon("play"))
               ),
               mainPanel(
                  plotOutput("pcaPlot"),
                  imageOutput("concordancePlot"),
                  downloadButton("downloadPCAPlot", "Download PCA Plot"),
                  
                  hr(),
                  h4("Example: PCA Input Format"),
                  tableOutput("examplePCA")
               )
            )
   ), # end of tabpanel
   
   ## STRUCTURE Analysis
   tabPanel(
      HTML("<span style='color:#ffffff;'>STRUCTURE Analysis</span>"),
      sidebarLayout(
         sidebarPanel(
            fileInput("structureFile", "Upload STRUCTURE Input"),
            numericInput("kMin", "Min K", value = 2, min = 1),
            numericInput("kMax", "Max K", value = 5, min = 1),
            numericInput("numKRep", "Replicates per K", value = 5, min = 1),
            numericInput("burnin", "Burn-in Period", value = 1000),
            numericInput("numreps", "MCMC Reps After Burn-in", value = 10000),
            checkboxInput("noadmix", "No Admixture Model", value = FALSE),
            checkboxInput("phased", "Phased Genotype", value = FALSE),
            numericInput("ploidy", "Ploidy Level", value = 2),
            checkboxInput("linkage", "Use Linkage Model", value = FALSE),
            actionButton("runStructure", "Run STRUCTURE", icon = icon("play"))
         ),
         mainPanel(
            h4("STRUCTURE Visualization"),
            imageOutput("structurePlot"),
            downloadButton("downloadAllResults", "Download All STRUCTURE Results (ZIP)"),
            downloadButton("downloadQmatrices", "Download Q-matrices"),
            downloadButton("downloadStructurePlot", "Download STRUCTURE Plot"),
            h4("Run Summary"),
            tableOutput("structureSummary"),
            h4("All STRUCTURE Plots"),
            uiOutput("structurePlots")
         )
      )
   ), #end of tabpanel
   p("© 2025 DNA Analysis Laboratory, Natural Sciences Research Institute, University of the Philippines Diliman. All rights reserved."),
   div(
      style = "position: fixed; bottom: 0, width: 100%; background-color: transparent; padding: 8px; text-align: center; font-size: 10px; color: #666;"
   )
   )

server <- function(input, output, session) {
      
      lastAction <- reactiveVal(Sys.time())
      
      # timer
      observe({
         invalidateLater(300000, session)  # 5 minutes
         
         # set cleanup time by 5 mins
         if (difftime(Sys.time(), lastAction(), units = "secs") > 300) {
            # resets
            concordanceResult(NULL)
            concordancePlotPath(NULL)
            fsnps_gen(NULL)
            population_stats(NULL)
            hardy_weinberg_stats(NULL)
            fst_stats(NULL)
            
            # this would reset the ui
            shinyjs::reset("formPanel")
            
            # to clean the files
            unlink(tempdir(), recursive = TRUE)
            
            showNotification("Session cleaned due to inactivity.", type = "message")
         }
      })
   
   
   ## Store plots for viewing and downloading
   plots <- reactiveValues(pca = NULL, het = NULL, fst = NULL)
   status <- reactiveVal("Waiting for input...")
   
   
   # ===================== FILE CONV
   convertedCSV <- reactiveVal(NULL)
   
   observe({
      file_ready <- FALSE
      pop_ready <- FALSE
      
      # Check file inputs
      if (input$inputType == "vcf") {
         file_ready <- !is.null(input$vcfFile)
      } else if (input$inputType == "plink") {
         file_ready <- !is.null(input$bedFile) &&
            !is.null(input$bimFile) &&
            !is.null(input$famFile)
      } else if (input$inputType == "bcf") {
         file_ready <- !is.null(input$bcfFile)
      }
      
      # Check population input
      if (input$poptype == "multiplepop") {
         pop_ready <- !is.null(input$multiplepop)
      } else if (input$poptype == "single") {
         pop_ready <- nchar(trimws(input$typePop)) > 0
      }
      
      # Enable or disable convert button
      if (file_ready && pop_ready) {
         shinyjs::enable("convertCSV")
      } else {
         shinyjs::disable("convertCSV")
      }
   })
   
   observeEvent(input$convertCSV, {
      
      lastAction(Sys.time())
      
      #csv_result <- vcftocsv(vcf = vcfPath, ref = refValue)
      #convertedCSV(csv_result)
      
      outputDir <- tempdir()
      outputName <- "converted_to_csv.csv"
      
      # Set reference input
      refValue <- if (input$poptype == "multiplepop") {
         input$multiplepop$datapath
      } else {
         input$typePop
      }
      
      if (input$inputType == "vcf") {
         ext <- tools::file_ext(input$vcfFile$name)
         if (ext %in% c("vcf", "gz")) {
            vcfPath <- input$vcfFile$datapath
            #vcftocsv(vcf = vcfPath, ref = refValue)
            
            csv_result <- vcftocsv(vcf = vcfPath, ref = refValue)
            convertedCSV(csv_result)
         } else {
            showNotification("Unsupported VCF format. Please upload a .vcf or .vcf.gz file.", type = "error")
            return()
         }
         
      } else if (input$inputType == "bcf") {
         temp_output <- file.path(outputDir, "bcftovcf")
         command <- stringr::str_c(
            plink_path, " --bcf ", input$bcfFile$datapath,
            " --const-fid 0 --cow --keep-allele-order --allow-no-sex --allow-extra-chr",
            " --recode vcf --out ", temp_output
         )
         system(command)
         vcfPath <- paste0(temp_output, ".vcf")
         #vcftocsv(vcf = vcfPath, ref = refValue)
         csv_result <- vcftocsv(vcf = vcfPath, ref = refValue)
         convertedCSV(csv_result)
         
      } else if (input$inputType == "plink") {
         bed <- input$bedFile$datapath
         bim <- input$bimFile$datapath
         fam <- input$famFile$datapath
         
         outputVCF <- file.path(outputDir, "plink2vcf")
         command <- stringr::str_c(
            plink_path, " --bed ", bed,
            " --bim ", bim, " --fam ", fam,
            " --const-fid 0 --cow --keep-allele-order --allow-no-sex --allow-extra-chr",
            " --recode vcf --out ", outputVCF
         )
         system(command)
         vcfPath <- paste0(outputVCF, ".vcf")
         #vcftocsv(vcf = vcfPath, ref = refValue)
         csv_result <- vcftocsv(vcf = vcfPath, ref = refValue)
         convertedCSV(csv_result)
      }
      
      output$downloadConvertedCSV <- downloadHandler(
         filename = function() { "converted_to_csv.csv" },
         content = function(file) {
            readr::write_csv(convertedCSV(), file)
         }
      )
      
      output$previewTable <- renderTable({
         req(convertedCSV())
         head(convertedCSV(), 10)  # Preview top 10 rows
      })
   }
   )
   
   ### For SNIPPER
   output$exampleTableSnipper <- renderTable({
      data.frame(
         Ind = c("sample1", "sample2", "sample3", "sample4"),
         rs101 = c("A/A", "A/T", "T/T", "A/T"),
         rs102 = c("G/C", "G/C", "G/G", "G/C"),
         rs103 = c("C/C", "C/G", "G/G", "G/G")
      )
   })
   
   output$exampleRefSnipper <- renderTable({
      data.frame(
         Sample.Name = c("sample1", "sample2", "sample3", "sample4"),
         Population = c("Malaysia", "Mexico", "Greece", "South Korea"),
         Superpopulation = c("Southeast Asia", "North and South America", "Europe", "East Asia")
      )
   })
   
   convertedSNIPPER <- reactiveVal(NULL)
   
   observe({
      hasFile <- !is.null(input$convertFile)
      hasRef <- !is.null(input$refFile)
      
      if (input$targetPop) {
         hasTarget <- nchar(trimws(input$targetPopName)) > 0
      } else {
         hasTarget <- TRUE  # Only required if checkbox is checked
      }
      
      ready <- hasFile && hasRef && hasTarget
      
      if (ready) {
         shinyjs::enable("convertBtn")
      } else {
         shinyjs::disable("convertBtn")
      }
   })
   
   observeEvent(input$convertBtn, {
      lastAction(Sys.time())
      
      inputPath <- input$convertFile$datapath
      refPath <- input$refFile$datapath
      targetSet <- input$targetPop
      targetName <- if (targetSet) input$targetPopName else NULL
      
      inputData <- read.csv(inputPath, header = TRUE)
      numMarkers <- ncol(inputData) - 2
      
      outputName <- "snipper.xlsx"
      
      snipper.file <- tosnipper(input = inputPath,
                                references = refPath,
                                target.pop = targetSet,
                                population.name = targetName,
                                markers = numMarkers)
      
      #xlsx_file <- file.path(paste(outputDir, outputName))
      convertedSNIPPER(snipper.file)
      
      output$downloadConverted <- downloadHandler(
         filename = function() { outputName },
         content = function(file) {
            file.copy(snipper.file, file)
         }
      )
      
      # try to preview the table
      output$previewTableSNIPPER <- renderTable({
         req(convertedSNIPPER())
         head(convertedSNIPPER(), 10)  # Preview top 10 rows
      })
      
   })
   
   ### UAS to CSV
   convertedUAS <- reactiveVal(NULL)
   
   output$exampleXLSX <- renderTable({
      data.frame(
         Sample.Name = c("sample1","sample1", "sample1", "sample1", "sample1", "sample1", "sample2", "sample3", "sample3", "sample3"),
         Locus = c("rs01", "rs01", "rs02", "rs02", "rs02", "rs03", "rs01", "rs01", "rs02", "rs03"),
         Allele = c("A", "T", "C", "A", "G", "T", "A", "T", "G", "A")
      )
   })
   
   observe({
      zip_ready <- !is.null(input$uas_zip)
      shinyjs::toggleState("run_uas2csv", condition = zip_ready)
   })
   
   observeEvent(input$run_uas2csv, {
      lastAction(Sys.time())
      
      req(input$uas_zip)
      
      temp_dir <- tempdir()
      input_path <- file.path(temp_dir, input$uas_zip$name)
      file.copy(input$uas_zip$datapath, input_path, overwrite = TRUE)
      
      # Handle optional reference input
      ref_value <- NULL
      use_reference <- FALSE
      
      if (!is.null(input$ref_file)) {
         ref_value <- input$ref_file$datapath
         use_reference <- TRUE
      }
      
      tryCatch({
         widened.file <- uas2csv(files = input_path,
                                 population = ref_value,
                                 reference = use_reference,
                                 dir = temp_dir)
         convertedUAS(widened.file)
         showNotification("Conversion complete!", type = "message")
      }, error = function(e) {
         showNotification(paste("Error:", e$message), type = "error")
      })
      
      outputName <- "01_merged_typed_data.csv"
      #csv_file <- file.path(paste(temp_dir, outputName))
      
      output$downloadUAScsv <- downloadHandler(
         filename = function() {
            outputName
         },
         content = function(file) {
            file.copy(widened.file, file)
         }
      )
      
      output$previewTableUAS <- renderTable({
         req(convertedUAS())
         head(convertedUAS(), 10)  # Preview top 10 rows
      })
      
   })
   
   ## Concordance Analysis
   observe({
      toggleState("compareBtn", !is.null(input$concordanceFile1) && !is.null(input$concordanceFile2))
   })
   
   concordanceResult <- reactiveVal(NULL)
   concordancePlotPath <- reactiveVal(NULL)
   
   observeEvent(input$compareBtn, {
      lastAction(Sys.time())
      
      tryCatch({
         
         req(input$concordanceFile1$datapath, input$concordanceFile2$datapath)
         
         haplo_flag <- input$isHaplotype
         file1_path <- input$concordanceFile1$datapath
         file2_path <- input$concordanceFile2$datapath
         
         result <- concordance(file1_path, file2_path, haplotypes = haplo_flag)
         
         concordanceResult(result$results)
         concordancePlotPath(result$plot)
         
         
      }, error = function(e) {
         showNotification(paste("Error during analysis:", e$message), type = "error", duration = 10)
      })
   })
   
   output$concordanceResults <- renderTable({
      req(concordanceResult())
      concordanceResult()
   })
   
   output$concordancePlot <- renderPlot({
      req(concordancePlotPath())
      concordancePlotPath()
   })
   
   output$downloadConcordance <- downloadHandler(
      filename = function() {"concordance.csv"},
      content = function(file) {
         file.copy("concordance.csv", file)
      }
   )
   
   output$concordancePlot <- renderImage({
      req(concordancePlotPath())
      
      list(
         src = concordancePlotPath(),
         contentType = "image/png",
         alt = "Concordance Plot"
      )
   }, deleteFile = FALSE)
   
   
   ## MARKER EXTRACTION
   output$exampleRSID <- renderTable({
      data.frame(
         rsID = c("rs101", "rs102", "rs103", "rs104")
      )
   })
   
   output$examplePOS <- renderTable({
      data.frame(
         Chromosome = c("1", "2"),
         Start_BP = c("104500", "205300"),
         End_BP = c("104700", "205700")
      )
   })
   
   output$exampleTable <- renderTable({
      data.frame(
         Ind = c("sample1", "sample2", "sample3"),
         rs101 = c("A/A", "A/T", "T/T"),
         rs102 = c("G/C", "G/C", "G/G"),
         rs103 = c("C/C", "C/G", "G/G")
      )
   })
   
   
   # START MARKER EXTRACTION
   observe({
      isFileUploaded <- !is.null(input$markerFile)
      
      isRSIDReady <- input$markerType == "rsid" && (
         (input$rsidInputType == "manual" && nzchar(input$typedRSIDs)) ||
            (input$rsidInputType == "upload" && !is.null(input$markerList1))
      )
      
      isPOSReady <- input$markerType == "pos" && !is.null(input$markerList2)
      
      toggleState("extractBtn", isFileUploaded && (isRSIDReady || isPOSReady))
   })
   
   withProgress(message = "Extracting markers...", value = 0, {
      observeEvent(input$extractBtn, {
         lastAction(Sys.time())
         
         req(input$markerFile)
         
         snps_list <- if (input$markerType == "rsid") {
            if (input$rsidInputType == "manual") {
               temp <- tempfile(fileext = ".txt")
               writeLines(strsplit(input$typedRSIDs, "\n")[[1]], temp)
               temp
            } else if (!is.null(input$markerList1)) {
               input$markerList1$datapath
            }
         } else NULL
         
         pos_list <- if (input$markerType == "pos" && !is.null(input$markerList2)) {
            ext <- tools::file_ext(input$markerList2$name)
            if (ext == "csv") read.csv(input$markerList2$datapath, header = FALSE)
            else if (ext %in% c("xlsx", "xls")) readxl::read_excel(input$markerList2$datapath, col_names = FALSE)
            else {
               showNotification("Invalid file type", type = "error")
               return(NULL)
            }
         } else NULL
         
         tryCatch({
            # load other parameters
            plink_args <- if (!is.null(input$plink_args) && nzchar(input$plink_args)) {
               strsplit(input$plink_args, "\\s+")[[1]]
            } else NULL
            

            temp_dir <- tempdir()
            extract_markers(
               input.file  = input$markerFile$datapath,
               snps.list = snps_list,
               pos.list = pos_list,
               bed.file = input$bedFile$datapath,
               bim.file = input$bimFile$datapath,
               fam.file = input$famFile$datapath,
               plink_args = plink_args,
               output.dir  = temp_dir,
               merged.file = "final_merged.vcf",
               plink_path  = plink_path
            )
            
            output$downloadExtracted <- downloadHandler(
               filename = function() { "final_merged.vcf" },
               content = function(file) {
                  source_path <- file.path(temp_dir, "final_merged.vcf")
                  file.copy(source_path, file)
               }
            )
         }, error = function(e) {
            showNotification(paste("Error:", e$message), type = "error")
         })
      })
   })
   
   
   # POP STAT
   output$examplePop <- renderTable({
      data.frame(
         Sample = c("Sample1", "Sample2", "Sample3", "Sample4"),
         Population = c("POP1", "POP2", "POP3", "POP4"),
         rs101 = c("A/A", "A/T", "A/A", "T/T"),
         rs102 = c("G/G", "C/C", "G/C", "G/G")
      )
   })
   
   observe({
      file_ready <- !is.null(input$popStatsFile)
      shinyjs::toggleState("runPopStats", condition = file_ready)
   })
   
   observeEvent(input$runPopStats, {
      lastAction(Sys.time())
      
      req(input$popStatsFile)
      
      fsnps_gen <- reactive({
         req(input$popStatsFile)
         df <- load_input_file(input$popStatsFile$datapath)
         cleaned <- clean_input_data(df)
         convert_to_genind(cleaned)
      })
      
      withProgress(message = "Running population analysis...", value = 0, {
         
         tryCatch({
            req(fsnps_gen())
            
            incProgress(0.4, detail = "Computing private alleles...")
            priv_alleles <- poppr::private_alleles(fsnps_gen())
            if (is.null(priv_alleles)) priv_alleles <- list(message = "No private alleles detected")
            
            incProgress(0.6, detail = "Computing population statistics...")
            population_stats <- reactive({
               req(fsnps_gen())
               compute_population_stats(fsnps_gen())
            })
            
            output$heterozygosity_table <- DT::renderDataTable({
               population_stats()$heterozygosity
            })
            
            output$inbreeding_table <- DT::renderDataTable({
               population_stats()$inbreeding_coeff
            })
            
            output$ttest_table <- DT::renderDataTable({
               population_stats()$ttest
            })
            
            output$allele_freq_table <- DT::renderDataTable({
               population_stats()$allele_frequencies
            }, options = list(scrollX = TRUE))
            
            output$privateAllelePlot <- renderUI({
               if (is.list(priv_alleles) && "message" %in% names(priv_alleles)) {
                  # Display styled message if no private alleles
                  tags$div(style = "color:gray; font-style:italic; margin-top:10px;",
                           priv_alleles$message)
               } else {
                  # Display as a nicely formatted table
                  DT::dataTableOutput("privateAlleleTable")
               }
            })
            
            output$privateAlleleTable <- DT::renderDataTable({
               as.data.frame(priv_alleles)
            }, options = list(pageLength = 10, scrollX = TRUE))
            
            ## Heterozygosity Plotting
            output$heterozygosity_plot <- renderImage({
               req(population_stats()$heterozygosity)
               
               plot_path <- plot_heterozygosity(
                  Het_fsnps_df = population_stats()$heterozygosity,
                  out_dir = tempdir()
               )
               
               list(
                  src = plot_path,
                  contentType = "image/png",
                  alt = "Heterozygosity Plot",
                  width = "100%"
               )
            }, deleteFile = TRUE)
            
            output$downloadHeterozygosityPlot <- downloadHandler(
               filename = function() {
                  "heterozygosity_plot.png"
               },
               content = function(file) {
                  # Recompute the plot to save as requested file
                  plot_path <- plot_heterozygosity(
                     Het_fsnps_df = population_stats()$heterozygosity,
                     out_dir = tempdir()
                  )
                  file.copy(plot_path, file)
               }
            )
            
            ## HWE
            incProgress(0.8, detail = "Running HWE and FST calculations...")
            hardy_weinberg_stats <- reactive({
               req(fsnps_gen())
               compute_hardy_weinberg(fsnps_gen())
            })
            
            output$hwe_summary <- renderUI({
               pvals <- hardy_weinberg_stats()$hw_summary
               tagList(
                  h4("HWE P-values across loci"),
                  verbatimTextOutput("hwe_summary_text")
               )
            })
            
            output$hwe_summary_text <- renderText({
               paste0("P-values: ", paste(round(hardy_weinberg_stats()$hw_summary, 4), collapse = ", "))
            })
            
            output$hwe_chisq_table <- DT::renderDataTable({
               hardy_weinberg_stats()$hw_dataframe
            }, options = list(scrollX = TRUE))
            
            ## FST
            incProgress(1.0, detail = "Still running HWE and FST calculations...")
            fst_stats <- reactive({
               req(fsnps_gen())
               compute_fst(fsnps_gen())
            })
            
            output$fstMatrixUI <- renderUI({
               fst <- fst_stats()$fsnps_fst_matrix
               if (is.list(fst) && "message" %in% names(fst)) {
                  tags$p(style = "color:gray;", fst$message)
               } else {
                  DT::dataTableOutput("fstMatrixTable")
               }
            })
            
            output$fstMatrixTable <- DT::renderDataTable({
               matrix_data <- matrix(unlist(fst_stats()$fst_matrix),
                                     nrow = sqrt(length(fst_stats()$fst_matrix)),
                                     byrow = TRUE)
               rownames(matrix_data) <- colnames(matrix_data) <- attr(fsnps_gen(), "pop.names")
               as.data.frame(matrix_data)
            }, options = list(scrollX = TRUE))
            
            output$fstDfTable <- DT::renderDataTable({
               fst_stats()$fst_dataframe
            })
            
            # plotting
            output$fst_heatmap_interactive <- plotly::renderPlotly({
               req(fst_stats()$fst_dataframe)
               plot_fst_heatmap_interactive(fst_stats()$fst_dataframe)
            })
            
            output$downloadFstHeatmap <- downloadHandler(
               filename = function() { "fst_heatmap.png" },
               content = function(file) {
                  plot_path <- plot_fst_heatmap(fst_stats()$fst_dataframe, out_dir = tempdir())
                  file.copy(plot_path, file)
               }
            )
            
            ## download all results
            stats_matrix <- compute_population_stats(fsnps_gen())
            hw_matrix <- compute_hardy_weinberg(fsnps_gen())
            fst_matrix <- compute_fst(fsnps_gen())
            
            output$downloadStatsXLSX <- downloadHandler(
               filename = function() {
                  timestamp <- format(Sys.time(), "%Y%m%d_%H%M")
                  paste0("population-statistics-results_", timestamp, ".xlsx")
               },
               content = function(file) {
                  req(stats_matrix, hw_matrix, fst_matrix)
                  path <- export_results(stats_matrix, hw_matrix, fst_matrix, dir = tempdir())
                  file.copy(path, file)
               }
            )
            
            incProgress(1, detail = "Finalizing output...")
            
         }, error = function(e) {
            showNotification(paste("Population stats error:", e$message), type = "error")
         })
      }) #end of withProgress
   }) 
   
   # PCA
   output$examplePCA <- renderTable({
      data.frame(
         Sample = c("Sample1", "Sample2", "Sample3", "Sample4"),
         Population = c("POP1", "POP2", "POP3", "POP4"),
         rs101 = c("A/A", "A/T", "A/A", "T/T"),
         rs102 = c("G/G", "C/C", "G/C", "G/G")
      )
   })
   
   observe({
      hasDataFile <- !is.null(input$pcaFile)
      
      hasLabelsFile <- !is.null(input$pcaLabels)
      hasColorFile <- !is.null(input$colorPalette)
      
      usingDefaults <- input$useDefaultColors
      
      readyForPCA <- hasDataFile && (usingDefaults || (hasLabelsFile && hasColorFile))
      
      toggleState("runPCA", readyForPCA)
   })
   
   observeEvent(input$runPCA, {
      lastAction(Sys.time())
      
      req(input$pcaFile)
      
      withProgress(message = "Running PCA...", {
         tryCatch({
            incProgress(0.2, detail = "Loading input file...")
            df <- load_input_file(input$pcaFile$datapath)
            cleaned <- clean_input_data(df)
            fsnps_gen <- convert_to_genind(cleaned)
            
            incProgress(0.4, detail = "Preparing color and label sets...")
            
            # Read custom labels/colors if needed
            labels <- NULL
            colors <- NULL
            
            if (!input$useDefaultColors) {
               req(input$pcaLabels, input$colorPalette)
               
               labels <- readLines(input$pcaLabels$datapath)
               colors <- readLines(input$colorPalette$datapath)
               
               # Optional: trim whitespace, ensure length match
               labels <- trimws(labels)
               colors <- trimws(colors)
               
               if (length(labels) != length(colors)) {
                  stop("Mismatch between number of labels and colors.")
               }
            }
            
            labels_colors <- get_colors_labels(
               fsnps_gen = fsnps_gen,
               use_default = input$useDefaultColors,
               input_labels = labels,
               input_colors = colors
            )
            
            incProgress(0.6, detail = "Computing PCA...")
            pca_results <- compute_pca(fsnps_gen)
            
            incProgress(0.8, detail = "Rendering PCA plot...")
            
            output$pcaPlot <- renderPlot({
               plot_pca(
                  ind_coords = pca_results$ind_coords,
                  centroid = pca_results$centroid,
                  percent = pca_results$percent,
                  labels_colors = labels_colors,
                  filename = NULL,  # rendering only
                  pc_x = input$pcX,
                  pc_y = input$pcY
               )
            })
            
            output$downloadPCAPlot <- downloadHandler(
               filename = function() {
                  paste0("pca_plot_", Sys.Date(), ".png")
               },
               content = function(file) {
                  plot_pca(
                     ind_coords = pca_results$ind_coords,
                     centroid = pca_results$centroid,
                     percent = pca_results$percent,
                     labels_colors = labels_colors,
                     filename = file,
                     pc_x = input$pcX,
                     pc_y = input$pcY
                  )
               }
            )
            
         }, error = function(e) {
            showNotification(paste("PCA Error:", e$message), type = "error")
         })
      })
   })
   
   ## STRUCTURE ANALYSIS
   observe({
      file_ready <- !is.null(input$structureFile)
      shinyjs::toggleState("runStructure", condition = file_ready)
   })
   
   
   observeEvent(input$runStructure, {
      lastAction(Sys.time())
      
      req(input$structureFile)
      
      temp_dir <- tempdir()
      
      result <- run_structure(
         file = input$structureFile$datapath,
         k.range = input$kMin:input$kMax,
         num.k.rep = input$numKRep,
         burnin = input$burnin,
         numreps = input$numreps,
         noadmix = input$noadmix,
         phased = input$phased,
         ploidy = input$ploidy,
         linkage = input$linkage,
         dir = temp_dir,
         structure_path = structure_path
      )
      
      zip_path <- file.path(result, "structure_results.zip")
      zip(zip_path, files = list.files(result$directory, full.names = TRUE))
      
      output$downloadAllResults <- downloadHandler(
         filename = function() {
            paste("STRUCTURE_results_", Sys.Date(), ".zip", sep = "")
         },
         content = function(file) {
            file.copy(zip_path, file)
         }
      )
      
      plots <- file.path(result, "str_plots.zip")
      zip(plots, files = list.files(result$directory, pattern = "\\.png$", full.names = TRUE))
      
      output$downloadStructurePlot <- downloadHandler(
         filename = function() {
            paste("STR_plots_", Sys.Date(), ".zip", sep = "")
         },
         content = function(file){
            file.copy(plots, file)
         }
      )
      
      matrices <- file.path(result, "matrices.zip")
      zip(matices, files = list.files(results$directory, pattern = "\\matrices", full.names = TRUE))
      output$downloadQmatrices <- downloadHandler(
         filename = function(){
            paste("QMatrices", Sys.Date(), ".zip", sep = "")
         },
         content = function(file){
            file.copy(matrices, file)
         }
      )
      
      output$structurePlots <- renderUI({
         plot_files <- list.files(result$directory, pattern = "\\.png$", full.names = TRUE)
         
         # Generate an imageOutput for each plot
         lapply(seq_along(plot_files), function(i) {
            local({
               f <- plot_files[i]
               id <- paste0("plot_", i)
               
               output[[id]] <- renderImage({
                  list(src = f, contentType = "image/png", width = "100%")
               }, deleteFile = FALSE)
               
               imageOutput(id, width = "100%", height = "auto")
            })
         })
      })
   })
   
}

shinyApp(ui, server)