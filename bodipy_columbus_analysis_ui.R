bodipyColumbusAnalysisUI <- function() {
  tabPanel("BODIPY Columbus Analysis",
           sidebarPanel(
             tags$h3("Input:"),
             textInput("Processed_file_name", "Analysis Name", 
                       value = paste("processed_data_", Sys.Date(),sep = "")),
             fileInput("Columbus_Data_Files", "input Columbus Output files as .zip", multiple = FALSE, accept = NULL, width = NULL),
             fileInput("Wellmap_Files", "input Columbus Output files as .zip", multiple = FALSE, accept = NULL, width = NULL),
             fluidRow(
               column(6, textInput("Control_Dose", "Dose of Control used in wellmaps", 
                                   value = NULL, placeholder = 0)),
               column(6, textInput("Control_Name_Input", "Control name used in wellmaps", 
                                   value = NULL, placeholder = "DMSO"))
             ),
             actionButton("submit", "Submit"),
             fileInput("temp_data_file", "input partial processed data", multiple = FALSE, accept = NULL, width = NULL),
             actionButton("temp_data_submit", "Submit Temp Data"),
           ), # sidebarPanel
           mainPanel(
             uiOutput("tables")
             
           ) # mainPanel
           
  )
}