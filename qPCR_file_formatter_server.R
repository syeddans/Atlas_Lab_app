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

  # Function to process single set of files
  process_qpcr_files <- function(wellmap_file, genemap_file, ct_file, repmap_file) {
    # Read the files
    wellmap <- read.csv(wellmap_file, header = FALSE)
    genemap <- read.csv(genemap_file, header = FALSE)
    ctmap <- read.csv(ct_file, header = FALSE)
    repmap <- read.csv(repmap_file, header = FALSE)  # New replicate map file
    
    # Create grid of row and column indices
    rows <- 2:nrow(wellmap)
    cols <- 2:ncol(wellmap)
    grid <- expand.grid(row = rows, col = cols)
    
    # Vectorized operations using dplyr
    formatted_data <- dplyr::tibble(
      well = paste0(wellmap[grid$row, 1], wellmap[1, grid$col]),
      chemical = as.character(wellmap[cbind(grid$row, grid$col)]),
      gene = as.character(genemap[cbind(grid$row, grid$col)]),
      CTvalue = as.numeric(ctmap[cbind(grid$row, grid$col)]),
      replicate = as.numeric(repmap[cbind(grid$row, grid$col)])  # Use replicate map
    ) %>%
    dplyr::filter(!is.na(CTvalue)) %>%
    dplyr::mutate(chemical = tolower(gsub(" ", "", chemical)))
    
    return(formatted_data)
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
    
    message("\nFound main folder:", basename(main_folder))
    
    # Find all rep folders within the main folder
    rep_folders <- list.dirs(main_folder, full.names = TRUE, recursive = FALSE)
    rep_folders <- rep_folders[!grepl("__MACOSX", rep_folders)]
    
    if (length(rep_folders) == 0) {
      stop("No valid rep folders found in the uploaded zip file")
    }
    
    message("\nFound rep folders:")
    message(paste(basename(rep_folders), collapse = "\n"))
    
    # Process each rep folder
    all_results <- lapply(rep_folders, function(rep_folder) {
      message("\nProcessing folder: ", basename(rep_folder))
      
      # Get all plate directories within this rep
      plate_dirs <- list.dirs(rep_folder, full.names = TRUE, recursive = FALSE)
      
      if (length(plate_dirs) == 0) {
        message("No plate directories found in folder ", basename(rep_folder))
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
        repmap_file <- list.files(plate_dir, pattern = ".*repmap\\.csv$", 
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
        message("No valid results for folder ", basename(rep_folder))
        return(NULL)
      }
      
      dplyr::bind_rows(valid_results)
    })
    
    # Remove NULL results and combine all reps
    valid_results <- Filter(Negate(is.null), all_results)
    if (length(valid_results) == 0) {
      stop("No valid results found in any rep folder")
    }
    
    # Combine all results
    all_formatted_data <- dplyr::bind_rows(valid_results)
    
    if (!is.data.frame(all_formatted_data) || nrow(all_formatted_data) == 0) {
      stop("Failed to create valid formatted data")
    }
    
    # Clean up
    unlink(qpcr_folder_dir, recursive = TRUE)
    message("Processing complete. Found ", nrow(all_formatted_data), " rows of data")
    
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
          input$wellmap_file$datapath,
          input$genemap_file$datapath,
          input$ct_file$datapath,
          input$repmap_file$datapath  # Keep repmap_file input for individual file processing
        )
      }

      # Save to temp directory
      output_file <- file.path(temp_dir, "formatted_qpcr_data.csv")
      write.csv(formatted_data, file = output_file, row.names = FALSE)

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
