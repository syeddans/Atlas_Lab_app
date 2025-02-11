qpcr_file_formatter_server <- function(input, output, session) {
  # Create a temporary directory to store the extracted files
  temp_dir <- tempdir()
  results_dir <- paste0(temp_dir, "/output/")
  formatted_data_dir <- paste0(
    "~/storage/Atlas_lab_app/qPCR_formatted_data/"
  )
  # Check if directories exist, create if not
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  if (!dir.exists(formatted_data_dir)) {
    dir.create(formatted_data_dir, recursive = TRUE)
  }

  # Ensure temp directory is deleted when session ends
  session$onSessionEnded(function() {
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
      print("temp dir deleted")
      temp_dir <- tempdir()
    }
  })

  # Function to process uploaded files and format data
  process_qpcr_files <- function(
    wellmap_file,
    genemap_file,
    ct_file,
    repmap_file
  ) {
    # Read the files - handle both direct file paths and uploaded file objects
    wellmap <- read.csv(if(is.character(wellmap_file)) wellmap_file else wellmap_file$datapath, header = FALSE)
    genemap <- read.csv(if(is.character(genemap_file)) genemap_file else genemap_file$datapath, header = FALSE)
    ctmap <- read.csv(if(is.character(ct_file)) ct_file else ct_file$datapath, header = FALSE)
    repmap <- read.csv(if(is.character(repmap_file)) repmap_file else repmap_file$datapath, header = FALSE)
    
    # Create grid of row and column indices
    rows <- 2:nrow(wellmap)
    cols <- 2:ncol(wellmap)
    grid <- expand.grid(row = rows, col = cols)
    
    # Validate dimensions and create formatted data
    tryCatch({
      # Check for out of bounds indices
      if (any(grid$row > nrow(wellmap)) || any(grid$col > ncol(wellmap))) {
        stop("Grid indices exceed wellmap dimensions")
      }
      if (any(grid$row > nrow(genemap)) || any(grid$col > ncol(genemap))) {
        stop("Grid indices exceed genemap dimensions")
      }
      if (any(grid$row > nrow(ctmap)) || any(grid$col > ncol(ctmap))) {
        stop("Grid indices exceed ctmap dimensions")
      }
      if (any(grid$row > nrow(repmap)) || any(grid$col > ncol(repmap))) {
        stop("Grid indices exceed repmap dimensions")
      }

      formatted_data <- dplyr::tibble(
        well = paste0(wellmap[grid$row, 1], wellmap[1, grid$col]),
        chemical = as.character(wellmap[cbind(grid$row, grid$col)]),
        gene = as.character(genemap[cbind(grid$row, grid$col)]),
        CTvalue = as.numeric(ctmap[cbind(grid$row, grid$col)]),
        replicate = as.numeric(repmap[cbind(grid$row, grid$col)])
      ) %>%
      dplyr::filter(!is.na(CTvalue)) %>%
      dplyr::mutate(chemical = tolower(gsub(" ", "", chemical)))
      
      return(formatted_data)
      
    }, error = function(e) {
      # Only print debugging info when error occurs
      cat("\nError in data formatting. Debugging information:\n")
      cat("Error message:", conditionMessage(e), "\n")
      cat("\nDimensions of maps:\n")
      cat("wellmap:", dim(wellmap), "\n")
      cat("genemap:", dim(genemap), "\n")
      cat("ctmap:", dim(ctmap), "\n")
      cat("repmap:", dim(repmap), "\n")
      cat("\nGrid dimensions:\n")
      cat("row range:", range(grid$row), "\n")
      cat("col range:", range(grid$col), "\n")
      
      return(NULL)
    })
  }

  # Function to process folder of folders using parallel processing
  process_qpcr_folder <- function(folder_path) {
    message("Starting folder processing...")
    message("Input folder path: ", folder_path)
    
    # Extract the main zip file
    qpcr_folder_dir <- paste0(temp_dir, "/qpcr_folder/")
    if (!dir.exists(qpcr_folder_dir)) {
      dir.create(qpcr_folder_dir, recursive = TRUE)
    }
    
    tryCatch({
      unzip(folder_path, exdir = qpcr_folder_dir)
      message("Main zip file extracted successfully")
    }, error = function(e) {
      message("Error extracting main zip: ", e$message)
      stop("Failed to extract main zip file: ", e$message)
    })
    
    # Find the main folder (excluding __MACOSX)
    main_folders <- list.dirs(qpcr_folder_dir, full.names = TRUE, recursive = FALSE)
    main_folders <- main_folders[!grepl("__MACOSX", main_folders)]
    main_folder <- main_folders[1]  # Take the first non-MACOSX folder
    # Process each rep folder
    all_results <- lapply(main_folder, function(main_folder) {
      message("\nProcessing folder: ", basename(main_folder))
      
      # Get all plate directories within this rep
      plate_dirs <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE)
      
      if (length(plate_dirs) == 0) {
        message("No plate directories found in folder ", basename(main_folder))
        return(NULL)
      }
      
      # Process each plate in this rep
      plate_results <- lapply(plate_dirs, function(plate_dir) {
        message("\nProcessing plate folder: ", basename(plate_dir))
        
        wellmap_file <- list.files(plate_dir, pattern = ".*wellmap\\.csv$", 
                                  full.names = TRUE, ignore.case = TRUE)[1]
        genemap_file <- list.files(plate_dir, pattern = ".*genemap\\.csv$", 
                                  full.names = TRUE, ignore.case = TRUE)[1]
        ct_file <- list.files(plate_dir, pattern = ".*CTmap\\.csv$", 
                             full.names = TRUE, ignore.case = TRUE)[1]
        repmap_file <- list.files(plate_dir, pattern = ".*(-repmap|-replicate)\\.csv$", 
                                 full.names = TRUE, ignore.case = TRUE)[1]  # New replicate map file
        
        if (is.na(wellmap_file) || is.na(genemap_file) || is.na(ct_file) || is.na(repmap_file)) {
          message("Skipping incomplete plate in: ", basename(plate_dir))
          return(NULL)
        }
        
        message("Processing files:")
        message("Wellmap: ", wellmap_file)
        message("Genemap: ", genemap_file)
        message("CT file: ", ct_file)
        message("Replicate map file: ", repmap_file)
        
        tryCatch({
          result <- process_qpcr_files(wellmap_file, genemap_file, ct_file, repmap_file)
          if (!is.data.frame(result)) {
            message("Warning: Invalid result format for plate ", basename(plate_dir))
            return(NULL)
          }
          return(result)
        }, error = function(e) {
          message("Error processing plate ", basename(plate_dir), ": ", e$message)
          return(NULL)
        })
      })
      
      # Remove NULL results and combine plates
      valid_results <- Filter(Negate(is.null), plate_results)
      if (length(valid_results) == 0) {
        message("No valid results for folder ", basename(main_folder))
        return(NULL)
      }
      
      dplyr::bind_rows(valid_results)
    })
    
    # Remove NULL results and combine all reps
    valid_results <- Filter(Negate(is.null), all_results)
    if (length(valid_results) == 0) {
      stop("No valid results found in any plate folder")
    }
    
    # Combine all results
    all_formatted_data <- dplyr::bind_rows(valid_results)
    
    if (!is.data.frame(all_formatted_data) || nrow(all_formatted_data) == 0) {
      stop("Failed to create valid formatted data")
    }
    
    # Clean up
    unlink(qpcr_folder_dir, recursive = TRUE)
    message("Processing complete. Found ", nrow(all_formatted_data), " rows of data")
    # Convert character columns to lowercase and remove whitespace
    all_formatted_data <- all_formatted_data %>%
      mutate(across(where(is.character), ~tolower(gsub("\\s+", "", .))))
    return(all_formatted_data)
  }

  # Handle file submission
  shiny::observeEvent(input$submit_format, {
    id <- shiny::showNotification("Processing files...", duration = NULL, type = "message")
    
    tryCatch({
      if (!is.null(input$qpcr_folder)) {
        formatted_data <- process_qpcr_folder(input$qpcr_folder$datapath)
      } else {
        formatted_data <- process_qpcr_files(
          input$wellmap_file,
          input$genemap_file,
          input$ct_file,
          input$repmap_file  # Keep repmap_file input for individual file processing
        )
      }

      # Save to temp directory
      output_file <- file.path(temp_dir, "formatted_qpcr_data.csv")
      write.csv(formatted_data, file = output_file, row.names = FALSE)
      # Find rows with any NA values
      NArows <- formatted_data[apply(formatted_data, 1, function(x) any(is.na(x))), ]
      
      if(nrow(NArows) > 0) {
        # Print to console
        message("Rows containing NA values:")
        print.data.frame(NArows, row.names = TRUE, max = NULL)
        
        # Create notification message
        na_message <- paste0(
          " ", nrow(NArows), " missing values found in the data inputted.\n",
          "Rows: ", paste(capture.output(print.data.frame(NArows, row.names = TRUE, max = NULL)), collapse="\n")
        )
        
        # Show notification
        shiny::showNotification(na_message, type = "warning", duration = NULL)
      } else {
        message("No NA values found in the data")
        shiny::showNotification("No NA values found in the data", type = "message")
      }

      output$download_section <- shiny::renderUI({
        shiny::downloadButton("download_formatted", "Download Formatted Data")
      })
      
      shiny::removeNotification(id)
      shiny::showNotification("Files processed successfully!", type = "message")
    }, error = function(e) {
      shiny::removeNotification(id)
      shiny::showNotification(
        paste("Error processing files:", e$message),
        type = "error"
      )
    })
  })

  # Simple download handler
  output$download_formatted <- shiny::downloadHandler(
    filename = function() {
      paste0("qPCR_formatted_data_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
    },
    content = function(file) {
      file.copy(file.path(temp_dir, "formatted_qpcr_data.csv"), file)
    }
  )
}
