# Load R packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(readr)
library(DT)
library("tools")
library(rmarkdown)
library(dplyr)

options(shiny.maxRequestSize = 800*1024^2)

source("bodipy_columbus_analysis_ui.R")
source("bodipy_columbus_analysis_server.R")
source("qPCR_western_analysis_ui.R")
source("qPCR_western_analysis_server.R")
source("qPCR_file_formatter_ui.R")
source("qPCR_file_formatter_server.R")

# Define UI
ui <- fluidPage(
  theme = shinytheme("cerulean"),
  navbarPage(
    "Atlas Lab",
    bodipyColumbusAnalysisUI(),
    qPCRWesternAnalysisUI(),
    qpcr_file_formatter_ui()
  )
)

# Define server
server <- function(input, output, session) {
  bodipyColumbusAnalysisServer(input, output, session)
  qpcr_western_analysis_server(input, output, session)
  qpcr_file_formatter_server(input, output, session)
}

# Run the app
shinyApp(ui = ui, server = server)