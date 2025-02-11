qpcr_western_analysis_server <- function(input, output, session) {
  # Create a temporary directory to store the extracted files
  temp_dir <- tempdir()
  results_dir <- paste0(temp_dir, "/output/")
  # Check if results directory exists, and create it if not
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  # Ensure temp directory is deleted when session ends
  session$onSessionEnded(function() {
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
      print("temp dir deleted")
      temp_dir <- tempdir()
    }
  })
  shiny::observeEvent(input$submit_qpcr, {
    shiny::req(input$qPCR_western_files)
    print("qPCR/Western analysis submit button pressed")
    id <- shiny::showNotification(
      "Processing, please wait...",
      duration = NULL,
      closeButton = FALSE,
      type = "message"
    )
    
    qpcr_western_files_dir <- paste0(temp_dir, "/qpcr_western_rep_dir/")
    if (!dir.exists(qpcr_western_files_dir)) {
      dir.create(qpcr_western_files_dir, recursive = TRUE)
    }
    print(paste("Temp directory path:", qpcr_western_files_dir))
    
    # Get original file name and create path
    original_filename <- input$qPCR_western_files$name
    files_path <- file.path(qpcr_western_files_dir, original_filename)
    
    # Copy file and verify
    file.copy(input$qPCR_western_files$datapath, files_path, overwrite = TRUE)
  
    print(input$Control_name_input)
    rmarkdown::render(
      "qPCR_western_Analysis/NoShinyqPCR results.Rmd",
      params = list(
        formatted_data_path = files_path,
        output_dir = temp_dir,
        control_gene = tolower(gsub("\\s+", "", input$Control_gene_input)),
        control_to_check = tolower(gsub("\\s+", "", input$Control_name_input)),
        analysis_type = input$qPCR_or_Western,
        #chemicalsymbols = strsplit(input$chemical_symbols, ",")[[1]],
        #chemicalnames = strsplit(input$chemical_name_full, ",")[[1]],
        stats_on = input$RQ_or_DDCT
      ),
      quiet = FALSE,
      envir = globalenv(),
      output_dir = temp_dir
    )
    print(list.files(temp_dir))
    ui_elements <- list()
    output$tables2 <- shiny::renderUI({
      ui_elements[[length(ui_elements) + 1]] <- shiny::downloadButton(
        "download_results_data",
        "Download Results"
      )
      shiny::tagList(ui_elements)
    })
    shiny::removeNotification(id)
  })
  output$download_results_data <- shiny::downloadHandler(
    filename = "results.zip",
    content = function(file) {
      zip_file <- tempfile(fileext = ".zip")
      zip(zip_file, c(
        file.path(temp_dir, "NoShinyqPCR-results.html"),
        file.path(temp_dir, "analysis_log.txt"),
        file.path(temp_dir, "input_data_formatted.csv"),
        file.path(temp_dir, "RQs.csv"),
        file.path(temp_dir, "Final_results_with_stats.csv"),
        file.path(temp_dir, "Data_averaged.csv"),
        file.path(temp_dir, "Calculated_RQ_by_rep.csv")
      ))
      file.copy(zip_file, file)
      unlink(zip_file)  # Clean up the temporary zip file
    }
  )
}