# Load R packages
library(shiny)
library(shinythemes)
library(shinyWidgets)
library(readr)
library(DT)
library("tools")
library(rmarkdown)

options(shiny.maxRequestSize = 800*1024^2)

source("bodipy_columbus_analysis_ui.R")
source("bodipy_columbus_analysis_server.R")

  # Define UI
ui <- fluidPage(theme = shinytheme("cerulean"),
                navbarPage(
                  "Atlas Lab",
                  bodipyColumbusAnalysisUI(),
                  tabPanel("Another Tab", "Content for another tab")
                ))

# Define server
server <- function(input, output, session) {
  bodipyColumbusAnalysisServer(input, output, session)
}

# Run the app
shinyApp(ui = ui, server = server)