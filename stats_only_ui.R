stats_only_ui <- function() {
  shiny::tabPanel(
    "Statistics Only",
    shiny::sidebarLayout(
      shiny::sidebarPanel(
        # File input for RQ data
        shiny::fileInput(
          "rq_data_file", 
          "Select RQ data file (CSV)",
          accept = c(
            "text/csv",
            "text/comma-separated-values,text/plain",
            ".csv"
          )
        ),
        
        # Text input for control condition
        shiny::textInput(
          "control_condition",
          "Control condition (e.g., 'mir+vehcon')",
          value = "mir+vehcon"
        ),
        
        shiny::actionButton("run_stats", "Run Statistical Analysis")
      ),
      shiny::mainPanel(
        shiny::uiOutput("stats_download_section"),
        DT::dataTableOutput("stats_results_table")
      )
    )
  )
}
