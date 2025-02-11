qPCRWesternAnalysisUI <- function() {
  tabPanel("qPCR/Western Analysis",
           sidebarPanel(
             tags$h3("Input:"),
             fileInput("qPCR_western_files", 
                      "Upload formatted qPCR/Western data file (.csv)", 
                      multiple = FALSE, 
                      accept = ".csv"),
             textInput("Control_gene_input", "control gene to compare", placeholder = NULL, width = NULL, value = "b-actin"),
             textInput("Control_name_input", "control name used in wellmaps", placeholder = "mir+vehcon", width = NULL, value = NULL),
             #textInput("chemical_symbols", "input as list if chemical symbols were used ", placeholder = "phe,9p,bap", width = NULL, value = NULL),
             #textInput("chemical_name_full", "corresponding full chemical names to symbols ", placeholder = "Phenanthrene,9-ChloroPhenanthrene,BaP", width = NULL, value = NULL),
             radioButtons("qPCR_or_Western", "Select an Option:", 
                          choices = c("qPCR" = "qPCR")),
                          #, "Western" = "western")),
             radioButtons("RQ_or_DDCT", "Select an Option:", 
                          choices = c("RQ" = "RQ")),
                          #, "DDCT" = "DDCT")),
             actionButton("submit_qpcr", "Submit")
           ), # sidebarPanel
           mainPanel(
             uiOutput("tables2"),
             #downloadButton("downloadData2", "Download Processed Data")
           ) # mainPanel
           
  )
}