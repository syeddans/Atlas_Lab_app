qPCRWesternAnalysisUI <- function() {
  tabPanel("qPCR/Western Analysis",
           sidebarPanel(
             tags$h3("Input:"),
             fileInput("qPCR_western_files", "Structure files as stated and upload as .zip", multiple = FALSE, accept = NULL, width = NULL),
             textInput("Control_gene_input", "control gene to compare", placeholder = NULL, width = NULL, value = "b-actin"),
             textInput("Control_name_input", "control name used in wellmaps", placeholder = NULL, width = NULL, value = "mir+vehcon"),
             textInput("chemical_symbols", "input as list if chemical symbols were used ", placeholder = NULL, width = NULL, value = "phe,9p,bap"),
             textInput("chemical_name_full", "corresponding full chemical names to symbols ", placeholder = NULL, width = NULL, value = "Phenanthrene,9-ChloroPhenanthrene,BaP"),
             radioButtons("qPCR_or_Western", "Select an Option:", 
                          choices = c("qPCR" = "qPCR", "Western" = "western")),
             radioButtons("RQ_or_DDCT", "Select an Option:", 
                          choices = c("RQ" = "RQ", "DDCT" = "DDCT")),
             actionButton("submit_qpcr", "Submit")
           ), # sidebarPanel
           mainPanel(
             uiOutput("tables2"),
             #downloadButton("downloadData2", "Download Processed Data")
           ) # mainPanel
           
  )
}