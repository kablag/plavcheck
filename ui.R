# Sys.setlocale("LC_ALL","Russian_Russia.20866")

library(shiny)
library(shinythemes)

shinyUI(
  tags$div(
    fluidPage(
      title = "ПлавЧек",
      theme = shinytheme("cerulean"),
      
      h1("ПлавЧек"),
      tabsetPanel(
        tabPanel("Главна\u044F",
                 wellPanel(
                   fluidRow(
                     column(2,
                            fileInput("rdmlFile",
                                      h4("Загрузите RDML:"))),
                     column(2,
                            numericInput("distanceLimit",
                                         "Порог рассто\u044Fни\u044F, > ",
                                         0.10, min = 0, max = 0.30, step = 0.01)),
                     column(2, offset = 2,
                            downloadButton("createReport",
                                           "Создать отчёт..."))))
        ),
        tabPanel("Подробно",
                 wellPanel(
                   h4("Показать:"),
                   fluidRow(
                     column(3,
                            uiOutput("targetViewSelector")),
                     column(3,
                            uiOutput("sampleViewSelector"))),
                   plotOutput("hclustPlot"),
                   plotOutput("meltPlot")))
      ),
      wellPanel(
        dataTableOutput("resultsTable")
      )
    )
  )
)