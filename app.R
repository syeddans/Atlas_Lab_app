# Load R packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(readr)
library(DT)
library("tools")
library(rmarkdown)

options(shiny.maxRequestSize = 800*1024^2)

#temp_dir = tempdir()
#unlink(temp_dir, recursive = TRUE)

  # Define UI
  ui <- fluidPage(theme = shinytheme("cerulean"),
                  
    navbarPage(
      "Atlas Lab",
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
               
      ), # Navbar 1, tabPanel
      #tabPanel("Navbar 2", "This panel is intentionally left blank"),
      #tabPanel("Navbar 3", "This panel is intentionally left blank")
  
    ) # navbarPage
  ) # fluidPage

  
  # Define server function  
  server <- function(input, output, session) {
    # Create a temporary directory to store the extracted files
    temp_dir <- tempdir()
    results_dir <<- paste0(temp_dir, "/output/")
    
    # Check if results directory exists, and create it if not
    if (!dir.exists(results_dir)) {
      dir.create(results_dir, recursive = TRUE)
    }
    
    biological_rep_folders <- NULL
    wellmap_folder <- NULL

    # Ensure temp directory is deleted when session ends
    session$onSessionEnded(function() {
      if (dir.exists(temp_dir)) {
        unlink(temp_dir, recursive = TRUE)
        print("temp dir deleted")
        # Reset temp_dir after deletion to ensure a fresh one is created
        temp_dir <<- tempdir()
      }
    })
    
    # Submit 1: Handles unzipping and rendering the initial biological replicate checkboxes to choose which is included in analysis
    observeEvent(input$submit, {
      req(input$Columbus_Data_Files)
      req(input$Wellmap_Files)
      req(input$Control_Name_Input)
      
      # Unzip files
      unzip(input$Columbus_Data_Files$datapath, exdir = temp_dir)
      unzip(input$Wellmap_Files$datapath, exdir = temp_dir)
      
      biological_rep_folders <<- list.dirs(
        paste0(temp_dir, "/", file_path_sans_ext(input$Columbus_Data_Files$name)), 
        recursive = FALSE, full.names = TRUE
      )
      wellmap_folder <<- paste0(temp_dir, "/", file_path_sans_ext(input$Wellmap_Files$name))

      insertUI(
        selector = "#submit",  # Place it after the initial submit button
        where = "afterEnd",
        ui = tagList(
          checkboxGroupInput("checkboxes_replicates", "Select Biological Replicates to Remove", 
                             choices = basename(biological_rep_folders)),
          actionButton("submit2", "Submit 2")
        )
      )
    })
    
    # Submit 2: shows general stats of each well to allow removal of technical replicates if needed 
    observeEvent(input$submit2, {
      # Clear the output UI before rendering new checkboxes for wells
      output$tables <- renderUI({})
      
      #notfication to show processing to user
      id <- showNotification("Processing, please wait...", duration = NULL, closeButton = FALSE)
      
      # Logic for excluding selected replicates
      for (folder in biological_rep_folders) {
        if (basename(folder) %in% input$checkboxes_replicates) {
          next  # Skip folders that should be excluded
        }
        system(paste("Rscript", shQuote("Columbus_Analysis/Columbus Number of Cells with Lipid pipeline.R"),
                     shQuote(folder), shQuote(results_dir), shQuote(wellmap_folder)))
      }
      removeNotification(id)
      
      # Render well checkboxes after processing files
      output$tables <- renderUI({
        ui_elements <- list()
        
        for (file in list.files(results_dir, full.names = TRUE)) {
          if (file_path_sans_ext(basename(file)) %in% file_path_sans_ext(input$checkboxes_replicates)) {
            next  # Skip folders that should be excluded
          }

          # Create custom labels for each well
          file_df <- read.table(file, sep = "\t", header = TRUE)
          checkbox_labels <- paste(
            file_df$WellName, ":", file_df$chemical, 
            "-", file_df$Cells.with.Lipid, "/", file_df$Number.of.Cells
          )
          
          # Add checkbox group for wells
          ui_elements[[length(ui_elements) + 1]] <- checkboxGroupInput(
            paste0("checkboxes_wells_", basename(file)),  # Unique input ID based on file name
            label = paste("Select Wells to exclude from", basename(file)),
            choices = checkbox_labels
          )
        }
        
        # Add a third submit button
        ui_elements[[length(ui_elements) + 1]] <- actionButton("submit3", "Submit 3")
        
        # Return the UI elements to display
        tagList(ui_elements)
      })
    })
    
    # Submit 3: Process the selected wells to generate shiny report
    observeEvent(input$submit3, {
      for (file in list.files(results_dir, full.names = TRUE)) { 
        file_df <- read.table(file, sep = "\t", header = TRUE)
        
        # Each group of checkboxes of the technical replicates is given a unique Id, this is to access it
        input_id <- paste0("checkboxes_wells_", basename(file))
        selected_values <- input[[input_id]]  # Get selected wells
        
        if (!is.null(selected_values)) {
          selected_well_names <- sapply(selected_values, function(x) {
            strsplit(x, ":")[[1]][1]  # Extract well name (e.g., "A1")
          })
          
          updated_file_df <- file_df[!file_df$WellName %in% selected_well_names, ]
          
          # Save updated file
          write.table(updated_file_df, file = paste0(results_dir, "/", basename(file)),
                      sep = "\t", row.names = FALSE)
        }
      }
      
      # Render the R Markdown report
      rmarkdown::render("Columbus_Analysis/Rshiny_file_output.Rmd",
                        params = list(results_dir = results_dir,
                                      control_name = input$Control_Name_Input),
                        output_dir = temp_dir)
                        #testing output_dir
                        #output_dir = "~/storage/Atlas_lab_app/Columbus_Analysis")
      id <- showNotification("Data Processed and ready to be downloaded", duration = NULL, closeButton = TRUE)
    })
    
    # Download processed file
    output$downloadData <- downloadHandler(
      # Define the file name for download 
      # TODO: maybe edit file name 
      filename = function() {
        paste("processed_data_", Sys.Date(), ".html", sep = "")
      },
      
      content = function(file) {
        # Copy the file from the temp directory to the download location
        file.copy(file.path(temp_dir, "Rshiny_file_output.html"), file)
      }
    )
    output$downloadInfo <- renderText({
      paste("Click the button to download the processed file.")
    })
    
    
    
  }
  
  # Create Shiny object
  shinyApp(ui = ui, server = server)
