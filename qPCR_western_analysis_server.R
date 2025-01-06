qPCRWesternAnalysisServer <- function(input, output, session) {
  # Create a temporary directory to store the extracted files
  temp_dir <- tempdir()
  results_dir <<- paste0(temp_dir, "/output/")
  
  # Check if results directory exists, and create it if not
  if (!dir.exists(results_dir)) {
    dir.create(results_dir, recursive = TRUE)
  }
  
  # Ensure temp directory is deleted when session ends
  session$onSessionEnded(function() {
    if (dir.exists(temp_dir)) {
      unlink(temp_dir, recursive = TRUE)
      print("temp dir deleted")
      # Reset temp_dir after deletion to ensure a fresh one is created
      temp_dir <<- tempdir()
    }
  })
  
  observeEvent(input$submit_qpcr, {
    req(input$qPCR_western_files)
    id <- showNotification("Processing, please wait...", duration = NULL, closeButton = FALSE)
    qPCR_western_files_dir <- paste0(temp_dir, "/qPCR_western_rep_dir/")
    if (!dir.exists(qPCR_western_files_dir)) {
      dir.create(qPCR_western_files_dir, recursive = TRUE)
    }
    unzip(input$qPCR_western_files$datapath, exdir = qPCR_western_files_dir)
    
    print(list.dirs(file.path(qPCR_western_files_dir, file_path_sans_ext(basename(input$qPCR_western_files$name))),recursive = TRUE))
    print(file.path(qPCR_western_files_dir, file_path_sans_ext(basename(input$qPCR_western_files$name))))
    print(paste0(qPCR_western_files_dir, file_path_sans_ext(basename(input$qPCR_western_files$name))))
    rmarkdown::render("qPCR_western_Analysis/NoShinyqPCR results.Rmd",
                      params = list(rep_files_dir= file.path(qPCR_western_files_dir, file_path_sans_ext(basename(input$qPCR_western_files$name))),
                                      output_dir= temp_dir,
                                      control_gene= input$Control_gene_input,
                                      control_to_check= input$Control_name_input,
                                      analysis_type= input$qPCR_or_Western, #qPCR or western
                                      chemicalsymbols= strsplit(input$chemical_symbols, ",")[[1]],
                                      chemicalnames= strsplit(input$chemical_name_full, ",")[[1]],
                                      stats_on= input$RQ_or_DDCT),  #RQ or DeltaCT
                                      quiet = FALSE,envir = globalenv(), output_dir = temp_dir)
 
    ui_elements <- list()
    output$tables2 <- renderUI({
      ui_elements[[length(ui_elements) + 1]] <- downloadButton("download_results_data", "Download Results")
      tagList(ui_elements)
    })
    #browser() 
    removeNotification(id)
     })

output$download_results_data <- downloadHandler(
  filename = "results.html",
  
  content = function(file) {
    # Copy the file from the temp directory to the download location
    file.copy(file.path(temp_dir, "NoShinyqPCR-results.html"), file)
  }
)
}