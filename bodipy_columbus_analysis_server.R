bodipyColumbusAnalysisServer <- function(input, output, session) {
  # Create a temporary directory to store the extracted files
  temp_dir <- tempdir()
  results_dir <<- paste0(temp_dir, "/output/")
  temp_analysis_dir <- paste0(temp_dir, "/temp_analysis_data/")
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
  
  # Observe the Temp Data Submit button
  observeEvent(input$temp_data_submit, { 
    req(input$temp_data_file)  # Ensure the temp data file input is available
    
    unzip(input$temp_data_file$datapath, exdir = results_dir)
    
    # Define the path to the temp_analysis_data folder
    temp_data_folder <- file.path(results_dir, file_path_sans_ext(input$temp_data_file$name))
    print("temp folder")
    print(temp_data_folder)
    print("inside temp folder")
    print(list.files(temp_data_folder, full.names = TRUE))
    # Move files from temp_analysis_data to results_dir
    file.copy(list.files(temp_data_folder, full.names = TRUE), results_dir, overwrite = TRUE)
    
    # Optionally, remove the temp_analysis_data folder if no longer needed
    unlink(temp_data_folder, recursive = TRUE)
    
    
    
    
    
    # Unzip the uploaded temp data file
    
    #file.copy(from = input$temp_data_file$datapath, to = results_dir)
    # Directly render the checkbox wells after unzipping
    render_checkbox_wells()
  })
  
  # Submit 2: Process the selected wells to generate shiny report
  observeEvent(input$submit2, {
    req(biological_rep_folders)  # Ensure the folders are available
    
    # Clear the output UI before rendering new checkboxes for wells
    output$tables <- renderUI({})
    
    # Notification to show processing to user
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
    
    print(list.files(results_dir, full.names = TRUE))
    dir.create(temp_analysis_dir)
    file.copy(from = list.files(results_dir, full.names = TRUE), to = temp_analysis_dir, overwrite = TRUE)
    print(list.files(temp_analysis_dir, full.names = TRUE))
    # Render well checkboxes after processing files
    render_checkbox_wells()
  })
  
  # Function to render well checkboxes
  render_checkbox_wells <- function() {
    output$tables <- renderUI({
      ui_elements <- list()
      #print(list.dirs(temp_dir))
      print("results dir")
      print(list.files(results_dir))
      for (file in list.files(results_dir, full.names = TRUE)) {
        if (!is.null(input$checkboxes_replicates)){ 
          if (file_path_sans_ext(basename(file)) %in% file_path_sans_ext(input$checkboxes_replicates)) {
            next  # Skip folders that should be excluded
          }
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
      ui_elements[[length(ui_elements) + 1]] <- downloadButton("download_temp_data", "Download Temp Data")
      
      # Return the UI elements to display
      tagList(ui_elements)
    })
  }
  
  # Submit 3: Process the selected wells to generate shiny report
  observeEvent(input$submit3, {
    for (file in list.files(results_dir, full.names = TRUE)) { 
      file_df <- read.table(file, sep = "\t", header = TRUE)
      
      # Each group of checkboxes of the technical replicates is given a unique Id, this is to access it
      input_id <- paste0("checkboxes_wells_", basename(file))
      selected_values <- input[[input_id]]  # Get selected wells
    
      if (!is.null(selected_values)) {
        selected_well_names <- sapply(selected_values, function(x) {
          trimws(strsplit(x, ":")[[1]][1])
        })
        print(selected_well_names)
        print("---------- selected values")
        print(selected_values)
        print(file_df)
        print("------------")
        updated_file_df <- file_df[!file_df$WellName %in% selected_well_names, ]
        print(updated_file_df)
        # Save updated file
        #TODO: currently cannot reselect data to reinclude after removing it because the file is being replaced
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
  
  output$download_temp_data <- downloadHandler(
    # Define the file name for download
    filename = function() {
      paste("temp_data", Sys.Date(), ".zip", sep = "")
    },
    
    content = function(file) {
      zip(file, basename(list.files(temp_analysis_dir), path = temp_analysis_dir))
    }
  )
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
}

