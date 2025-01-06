bodipyColumbusAnalysisServer <- function(input, output, session) {
  # Create a temporary directory to store the extracted files
  temp_dir <- tempdir()
  results_dir <<- paste0(temp_dir, "/output/")
  temp_analysis_dir <- paste0(temp_dir, "/temp_analysis_data/")
  ui_elements <- list()
  analysis_choices <- c(
    "Ratio of +/- Lipid (BODIPY) cells" = "1", 
    "Avg Number of Lipid per cell" = "2", 
    "Avg Lipid Intensity per cell" = "3", 
    "Mean Lipid Size" = "4"
  )
  # Check if results directory exists, and create it if not
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  if (!dir.exists(temp_analysis_dir)) {
    dir.create(temp_analysis_dir, recursive = TRUE)
  }
  
  biological_rep_folders <- NULL
  wellmap_folder <- NULL
  
  
  # Ensure temp directory is deleted when session ends
  session$onSessionEnded(function() {
    if (dir.exists(temp_dir)) {
      file.copy(from = temp_dir, to="~/storage/Atlas_lab_app/Extra/temp_dir_dump", recursive= TRUE)
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
    req(input$Processed_file_name)
    
    id <- showNotification("One moment", duration = NULL, closeButton = FALSE)
    # Unzip files
    biological_rep_dir <- paste0(temp_dir, "/biological_rep_dir/")
    if (!dir.exists(biological_rep_dir)) {
      dir.create(biological_rep_dir, recursive = TRUE)
    }
    unzip(input$Columbus_Data_Files$datapath, exdir = biological_rep_dir)
    
    wellmap_dir <- paste0(temp_dir, "/wellmap_dir/")
    if (!dir.exists(wellmap_dir)) {
      dir.create(wellmap_dir, recursive = TRUE)
    }
    unzip(input$Wellmap_Files$datapath, exdir = wellmap_dir, junkpaths = TRUE)

    biological_rep_folders <<- list.dirs(paste0(biological_rep_dir, "/", file_path_sans_ext(input$Columbus_Data_Files$name)), 
      recursive = FALSE, full.names = TRUE
    )
    wellmap_folder <<- wellmap_dir
    
    #print("temp dirs")
    #print(list.files(temp_dir, full.names = TRUE, recursive =  TRUE))
    
    insertUI(
      selector = "#submit",  # Place it after the initial submit button
      where = "afterEnd",
      ui = tagList(
        checkboxGroupInput("checkboxes_replicates", "Select Biological Replicates to Remove", 
                           choices = basename(biological_rep_folders)),
        actionButton("submit2", "Submit 2")
      )
    )
    removeNotification(id)
  })
  
  # Observe the Temp Data Submit button
  observeEvent(input$temp_data_submit, { 
    req(input$temp_data_file)  # Ensure the temp data file input is available
    req(input$Control_Name_Input)
    
    unzip(input$temp_data_file$datapath, exdir = temp_analysis_dir, junkpaths = TRUE)
    output$tables <- renderUI({
      ui_elements[[length(ui_elements) + 1]] <-selectInput("data_analysis_type", "Select Type of Analysis", choices = analysis_choices)
      ui_elements[[length(ui_elements) + 1]] <-selectInput("plot_type", "Select Plot Type", choices = c("bar"="bar",
                                                                                                        "box" = "box"))
      
      ui_elements[[length(ui_elements) + 1]] <- actionButton("submit3", "Submit 3")
      tagList(ui_elements)
    })
    #render_checkbox_wells()
    print(input$plot_type)
  })
  
  # Submit 2: Process the selected wells to run analysis
  observeEvent(input$submit2, {
    req(biological_rep_folders)  # Ensure the folders are available
    
    # Clear the output UI before rendering new checkboxes for wells
    output$tables <- renderUI({})
    
    # Notification to show processing to user
    id <- showNotification("Processing, please wait...", duration = NULL, closeButton = FALSE)
    
    # Logic for excluding selected replicates
    for (folder in biological_rep_folders) {
      if ((basename(folder) %in% input$checkboxes_replicates) || (basename(folder) %in% file_path_sans_ext(list.files(temp_analysis_dir)))){
        next  # Skip folders that should be excluded whether they were already prcoessed or choosen not to be prcoessed 
      }
      system(paste("Rscript", 
                   shQuote("Columbus_Analysis/Columbus Number of Cells with Lipid pipeline.R"),
                   shQuote(folder), 
                   shQuote(temp_analysis_dir), 
                   shQuote(wellmap_folder), 
                   input$data_analysis_type))
    }
    removeNotification(id)
  
    #file.copy(from = list.files(results_dir, full.names = TRUE), to = temp_analysis_dir, overwrite = TRUE)

    # Render well checkboxes after processing files
    output$tables <- renderUI({
      ui_elements[[length(ui_elements) + 1]] <-selectInput("data_analysis_type", "Select Type of Analysis", choices = analysis_choices)
      ui_elements[[length(ui_elements) + 1]] <- actionButton("submit3", "Submit 3")
      tagList(ui_elements)
    })
  })
  observeEvent(input$submit3, {
    
    render_checkbox_wells()
  })
  # Function to render well checkboxes
  render_checkbox_wells <- function() {
    output$tables <- renderUI({
      
      #print(list.dirs(temp_dir))
      for (file in list.files(temp_analysis_dir, full.names = TRUE)) {
        if (!is.null(input$checkboxes_replicates)){ 
          if (file_path_sans_ext(basename(file)) %in% file_path_sans_ext(input$checkboxes_replicates)) {
            next  # Skip folders that should be excluded
          }
        }
        
        # Create custom labels for each well
        file_df <- read.table(file, sep = "\t", header = TRUE)
        checkbox_labels <- paste(
          file_df$WellName, ":", file_df$chemical, "-",
          if (input$data_analysis_type == "1") {
            paste(file_df$Cells_with_Lipid, "/", file_df$Number_of_Cells)
          } else if (input$data_analysis_type == "2") {
            paste(file_df$Number_of_Lipid, "/", file_df$Number_of_Cells)
          } else if (input$data_analysis_type == "3") {
            paste(file_df$Avg_Lipid_Intensity)
          } else if (input$data_analysis_type == "4") {
            paste(file_df$Avg_Lipid_Area)
          }
        )
        
        # Add checkbox group for wells
        ui_elements[[length(ui_elements) + 1]] <- checkboxGroupInput(
          paste0("checkboxes_wells_", basename(file)),  # Unique input ID based on file name
          label = paste("Select Wells to exclude from", basename(file)),
          choices = checkbox_labels
        )
      }
      
      # Add a third submit button
      ui_elements[[length(ui_elements) + 1]] <- actionButton("submit4", "Submit 4")
      ui_elements[[length(ui_elements) + 1]] <- downloadButton("download_temp_data", "Download Temp Data")
      ui_elements[[length(ui_elements) + 1]] <- downloadButton("downloadData", "Download Processed Data")
      
      # Return the UI elements to display
      tagList(ui_elements)
    })
  }
  
  # Submit 3: Process the selected wells to generate shiny report
  observeEvent(input$submit4, {
    file.copy(from = list.files(temp_analysis_dir, full.names = TRUE), to = results_dir, recursive = TRUE)
    
    for (file in list.files(temp_analysis_dir, full.names = TRUE)) { 
      file_df <- read.table(file, sep = "\t", header = TRUE)
      
      # Each group of checkboxes of the technical replicates is given a unique Id, this is to access it
      input_id <- paste0("checkboxes_wells_", basename(file))
      selected_values <- input[[input_id]]  # Get selected wells
    
      if (!is.null(selected_values)) {
        selected_well_names <- sapply(selected_values, function(x) {
          trimws(strsplit(x, ":")[[1]][1])
        })

        updated_file_df <- file_df[!file_df$WellName %in% selected_well_names, ]

        # Save updated file
        write.table(updated_file_df, file = paste0(results_dir, "/", basename(file)),
                    sep = "\t", row.names = FALSE)
      }
    }
    print(list.files(results_dir))
    # Render the R Markdown report
    rmarkdown::render("Columbus_Analysis/Rshiny_file_output.Rmd",
                      params = list(title = input$Processed_file_name,
                                    results_dir = results_dir,
                                    control_name = input$Control_Name_Input,
                                    control_dose = input$Control_Dose),
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
      zip(file, temp_analysis_dir)
    }
  )
  # Download processed file
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0(input$Processed_file_name,".html", sep = "")
    },
    
    content = function(file) {
      # Copy the file from the temp directory to the download location
      file.copy(file.path(temp_dir, "Rshiny_file_output.html"), file)
    }
  )
}

