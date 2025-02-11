temp_dir <- tempdir()
output_dir <<- temp_dir
source("qPCR_western_Analysis/qPCR analysis Script_updated.R")
stats_only_server <- function(input, output, session) {
  
  # Ensure temp directory is deleted when session ends
  session$onSessionEnded(function() {
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
    }
  })
  
  # Function to process the RQ data file
  process_rq_data <- function(file_path, control_to_check) {
    tryCatch({
      # Read the data
      data <- read.csv(file_path)
      # Convert text columns to lowercase and remove whitespace
      text_cols <- sapply(data, is.character)
      data[text_cols] <- lapply(data[text_cols], function(x) {
        tolower(gsub("\\s+", "", x))
      })
      # Validate required columns
      required_cols <- c("chemical", "gene")
      rq_pattern <- "^RQ\\d*$"
      
      if (!all(required_cols %in% colnames(data))) {
        stop("Missing required columns: chemical and/or gene")
      }
      
      rq_cols <- grep(rq_pattern, colnames(data), value = TRUE)
      if (length(rq_cols) == 0) {
        stop("No RQ columns found in the data")
      }
      
      # Process the data through calculate_pvalues
      result <- process_CTdf(data, stats_on = "RQ", control_to_check = tolower(gsub("\\s+", "", input$control_condition)), analysistype = "qPCR", unit = "uM", gene_or_protein_label = "gene")
      
      return(result)
      
    }, error = function(e) {
      # Print debugging info when error occurs
      cat("\nError in data processing. Debugging information:\n")
      cat("Error message:", conditionMessage(e), "\n")
      return(NULL)
    })
  }
  
  # Handle the Run Statistics button click
  observeEvent(input$run_stats, {
    id <- shiny::showNotification("Running statistical analysis...", duration = NULL, type = "message")
    
    tryCatch({
      if (is.null(input$rq_data_file)) {
        stop("Please select a data file")
      }
      
      # Process the data
      results <- process_rq_data(
        input$rq_data_file$datapath,
        input$control_condition
      )
      
      if (is.null(results)) {
        stop("Error processing data")
      }
      
      # Save results to temp directory
      output_file <- file.path(temp_dir, "statistical_results.csv")
      write.csv(results, output_file, row.names = FALSE)
      
      # Display results table
      output$stats_results_table <- DT::renderDataTable({
        DT::datatable(
          results,
          options = list(
            pageLength = 25,
            scrollX = TRUE
          )
        )
      })
      
      # Enable download button
      output$stats_download_section <- shiny::renderUI({
        shiny::downloadButton("download_stats", "Download Results")
      })
      
      shiny::removeNotification(id)
      shiny::showNotification("Statistical analysis complete!", type = "message")
      
    }, error = function(e) {
      shiny::removeNotification(id)
      shiny::showNotification(
        paste("Error:", e$message),
        type = "error"
      )
    })
  })
  
  output$download_stats <- shiny::downloadHandler(
    filename = "results.zip",
    content = function(file) {
      zip_file <- tempfile(fileext = ".zip")
      zip(zip_file, c(
        file.path(temp_dir, "statistical_results.csv"),
        file.path(temp_dir, "analysis_log.txt")
      ))
      file.copy(zip_file, file)
      unlink(zip_file)  # Clean up the temporary zip file
    }
  )
}
