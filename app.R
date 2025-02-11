# Load R packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(readr)
library(DT)
library("tools")
library(rmarkdown)
library(dplyr)
library(markdown)

options(shiny.maxRequestSize = 800*1024^2)

source("bodipy_columbus_analysis_ui.R")
source("bodipy_columbus_analysis_server.R")
source("qPCR_western_analysis_ui.R")
source("qPCR_western_analysis_server.R")
source("qPCR_file_formatter_ui.R")
source("qPCR_file_formatter_server.R")
source("stats_only_ui.R")
source("stats_only_server.R")

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    "Atlas Lab",
    bodipyColumbusAnalysisUI(),
    qpcr_file_formatter_ui(),
    qPCRWesternAnalysisUI(),
    stats_only_ui(),
    tabPanel("How to use",
      fluidRow(
        column(12,
          includeMarkdown("README.md")
        )
      )
    )
  )
)

# Define server
server <- function(input, output, session) {
  bodipyColumbusAnalysisServer(input, output, session)
  qpcr_western_analysis_server(input, output, session)
  qpcr_file_formatter_server(input, output, session)
  stats_only_server(input, output, session)
  output$download_example1 <- downloadHandler(
    filename = function() {
      "Example1_qPCR_input_data.zip"
    },
    content = function(file) {
      file.copy("Examples/Example1_qPCR_input_data.zip", file)
    }
  )
  output$download_example2 <- downloadHandler(
    filename = function() {
      "Example2_qPCR_input_data.zip"
    },
    content = function(file) {
      file.copy("Examples/Example2_qPCR_input_data.zip", file)
    }
  )
  output$download_formatted_results1 <- downloadHandler(
    filename = function() {
      "Example_Formatted_Results1.csv"
    },
    content = function(file) {
      file.copy("Examples/Example_Formatted_Results1.csv", file)
    }
  )
  output$download_formatted_results2 <- downloadHandler(
    filename = function() {
      "Example_Formatted_Results2.csv"
    },
    content = function(file) {
      file.copy("Examples/Example_Formatted_Results2.csv", file)
    }
  )
  output$download_final_results1 <- downloadHandler(
    filename = function() {
      "Example_Results_Final1.zip"
    },
    content = function(file) {
      file.copy("Examples/Example_Results_Final1.zip", file)
    }
  )
  output$download_final_results2 <- downloadHandler(
    filename = function() {
      "Example_Results_Final2.zip"
    },
    content = function(file) {
      file.copy("Examples/Example_Results_Final2.zip", file)
    }
  )
  output$download_stats_format <- downloadHandler(
    filename = function() {
      "Stats_Only_Format.csv"
    },
    content = function(file) {
      file.copy("Examples/Stats_Only_Format.csv", file)
    }
  )
}

# Run the app
shinyApp(ui = ui, server = server)