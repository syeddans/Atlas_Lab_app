qpcr_file_formatter_ui <- function() {
  shiny::tabPanel(
    "qPCR File Formatter",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # File input selection
        shiny::radioButtons(
          "file_input_type",
          "Choose input type:",
          choices = c(
            "Single set of files" = "single",
            "Folder of folders" = "folder"
          )
        ),
        
        # Conditional UI based on input type
        shiny::conditionalPanel(
          condition = "input.file_input_type == 'single'",
          shiny::fileInput("wellmap_file", "Select wellmap file (CSV)"),
          shiny::fileInput("genemap_file", "Select genemap file (CSV)"),
          shiny::fileInput("ct_file", "Select CT file (CSV)"),
          shiny::fileInput("repmap_file", "Select replicate map file (CSV)")
        ),
        
        shiny::conditionalPanel(
          condition = "input.file_input_type == 'folder'",
          shiny::fileInput("qpcr_folder", "Select folder (ZIP)")
        ),
        
        shiny::actionButton("submit_format", "Process Files")
      ),
      shiny::mainPanel(
        shiny::uiOutput("download_section")
      )
    )
  )
}
