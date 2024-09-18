bodipyColumbusAnalysisUI <- function() {
  tabPanel("BODIPY Columbus Analysis",
           sidebarPanel(
             tags$h3("Input:"),
             fileInput("Columbus_Data_Files", "input Columbus Output files as .zip", multiple = FALSE, accept = NULL, width = NULL),
             fileInput("Wellmap_Files", "input Columbus Output files as .zip", multiple = FALSE, accept = NULL, width = NULL),
             textInput("Control_Name_Input", "control name used in wellmaps", value = NULL, width = NULL, placeholder = "DMSO"),
             actionButton("submit", "Submit")
           ), # sidebarPanel
           mainPanel(
             uiOutput("tables"),
             downloadButton("downloadData", "Download Processed Data")
           ) # mainPanel
           
  )
}