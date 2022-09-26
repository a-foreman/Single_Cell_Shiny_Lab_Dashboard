
library(shiny)
library(shinythemes)
library(Seurat)
library(uwot)
library(ggplot2)
library(dplyr)

#set mandatory fields
fieldsMandatory <- c("Name","Experiment.Name", "Type", "Genotype", "Stage", "Time","Tissue", "Protocol", "Replicates", "Pre.process.analysis", "Analysis", "RDS.name", "Folder.Name", "Raw.Data.ID")
fieldsAll <- c("Name","Experiment.Name", "Type", "Genotype", "Stage", "Time","Tissue", "Protocol", "Replicates", "Pre.process.analysis", "Analysis", "RDS.name", "Folder.Name", "Raw.Data.ID", "Publications.Document","Additional.Comments")


#load genenames.csv (update to your file!)
genenames <- read.csv("Annotation_Files/Homo_sapiens.GRCh38.107.csv", header = F)

#Read in datanames
files <- list.files(file.path("Responses"), full.names = TRUE)
data <- lapply(files, read.csv, stringsAsFactors = FALSE)
data <- do.call(rbind, data)
datanames <- data$RDS.name
datanames <- na.omit(datanames)
datanames <- unique(unlist(strsplit(datanames, " ")))

#Set names list for select input in choices
named_file_list = as.list(datanames)
names(named_file_list) <- as.list(datanames)

#Loop to read in files

file_contents <- list()

for(i in 1:length(datanames)){
  file_contents[[i]] <-readRDS(
    file = datanames[[i]]
  )
}
file_contents <- setNames(file_contents, datanames)

#Functions to set .csv filenames
responsesDir <- file.path("Responses")
epochTime <- function() {
  as.integer(Sys.time())
}
labelMandatory <- function(label) {
  tagList(
    label,
    span("*", class = "mandatory_star")
  )
}

#Set * on mandatory fields.
appCSS <- ".mandatory_star { color: red; }"
humanTime <- function() format(Sys.time(), "%Y%m%d-%H%M%OS")

#User interface. 
ui <-fluidPage(
  shinyjs::useShinyjs(),
  shinyjs::inlineCSS(appCSS), navbarPage("Brand Lab Dashboard", theme = shinytheme("cosmo"), 
                                         
                                         tabPanel(title = div(icon("chart-line", lib = "font-awesome"), "Single Cell"),
                                                  navlistPanel(
                                                    id = "tabset",
                                                    "Protocols", 
                                                    tabPanel("Analysis",
                                             #Insert your protocols text here                
                                                             p("INSERT YOUR TEXT HERE"), tags$a(href="https://satijalab.org/seurat/", "here!"),
                                                             p("INSERT YOUR TEXT HERE"),
                                                    ),
                                                    "Datasets",
                                                    tabPanel("Single Cell Lab Datasets",
                                                             mainPanel(
                                                               DT::dataTableOutput("responsesTable1"))
                                                    ),
                                                    tabPanel("Add new dataset",
                                                             mainPanel(
                                                               div(
                                                                 id = "form",
                                                                 textInput("Name", labelMandatory("Name"), ""),
                                                                 textInput("Experiment.Name", labelMandatory("Experiment Name"), ""),
                                                                 textInput("Type", labelMandatory("Type"),""),
                                                                 textInput("Genotype", labelMandatory("Genotype"),""),
                                                                 textInput("Stage", labelMandatory("Stage"),""),
                                                                 textInput("Time", labelMandatory("Time"),""),
                                                                 textInput("Tissue", labelMandatory("Tissue"),""),
                                                                 textInput("Protocol", labelMandatory("Protocol"),""),
                                                                 textInput("Replicates", labelMandatory("Replicates"),""),
                                                                 textInput("Pre.process.analysis", labelMandatory("Pre-process analysis"),""),
                                                                 textInput("Analysis", labelMandatory("Analysis"),""),
                                                                 textInput("RDS.name", labelMandatory("RDS Name"),""),
                                                                 textInput("Folder.Name", labelMandatory("Folder Name"),""),
                                                                 textInput("Raw.Data.ID", labelMandatory("Raw Data ID"),""),
                                                                 textInput("Publications.Document","Publications Document", ""),
                                                                 textInput("Additional.Comments", "Additional Comments", ""),
                                                                 actionButton("submit", "Submit", class = "btn-primary")),
                                                               shinyjs::hidden(
                                                                 div(
                                                                   id = "thankyou_msg",
                                                                   h3("Thanks, your response was submitted successfully!"),
                                                                   actionLink("submit_another", "Submit another response")
                                                                 )
                                                               ),  
                                                               DT::dataTableOutput("responsesTable")
                                                             )),
                                                    "sc Explorer",
                                                    tabPanel("Umap Browser", 
                                                             titlePanel("Umap Browser"),
                                                             fluidRow(
                                                               column(8,
                                                                      selectInput("dataset", label = h3("Dataset"),
                                                                                  choices = named_file_list,
                                                                                  selected = "pbmc.rds"),
                                                                      helpText("Gene names must be exact."),
                                                                      selectizeInput("genesearch1", label = "Gene Name", choices = NULL, multiple = TRUE, options = NULL)
                                                               )),
                                                             fluidRow(
                                                               column(12, 
                                                                      plotOutput("dimPlot1", width = "200%"))),
                                                             fluidRow(column(12,         
                                                                             plotOutput("umapPlot", width = "200%"))),
                                                             fluidRow(column(3, 
                                                                             downloadButton("downloadumap", "Download.png")))),
                                                   
                                                    
#This section is to show Violin plots of selected gene names. UPDATE: the selected = ".rds" to your favorite .rds file. This will be automatically selected when the app loads.                   

                                                 tabPanel("Violin Plots", 
                                                             titlePanel("Violin Plots"),
                                                             fluidRow(
                                                               column(8,
                                                                      selectInput("dataset", label = h3("Dataset"),
                                                                                  choices = named_file_list,
                                                                                  selected = "pbmc"),
                                                                      helpText("Gene names must be exact."),
                                                                      selectizeInput("genesearch2", label = "Gene Name", choices = NULL, multiple = FALSE, options = NULL)
                                                               )),
                                                             fluidRow(
                                                               column(8, 
                                                                      plotOutput("vlnPlot", width = "200%"))),
                                                             fluidRow(column(4,
                                                                             downloadButton("downloadvln", "Download.png")))),
                                                    tabPanel("Mulit-gene Browser",
                                                             titlePanel("Mulit-gene Browser"),
                                                             fluidRow(
                                                               column(8,
                                                                      selectInput("dataset", label = h3("Dataset"),
                                                                                  choices = named_file_list,
                                                                                  selected = "pbmc.rds"),
                                                                      helpText("Enter Gene names must be exact."),
                                                                      selectizeInput("genesearch3", label = "Gene Name", choices = NULL, multiple = TRUE, options = NULL)
                                                               )),
                                                             fluidRow(
                                                               column(8, 
                                                                      plotOutput("DotPlot", width = "200%"))),
                                                             fluidRow(column(4,
      
#Server Section                                                                              
                                                                                                                                                    downloadButton("downloaddot", "Download.png"))))))))
server <- function(input, output, session) {
  thematic::thematic_shiny()
  #Load prev responses
  loadData <- function() {
    files <- list.files(file.path(responsesDir), full.names = TRUE)
    data <- lapply(files, read.csv, stringsAsFactors = FALSE)
    data <- do.call(rbind, data)
    data
  }
  #Plot functions
  datasetInput <- reactive({
    df <- (input$dataset)
  })
  output$vlnPlot <- renderPlot({
    VlnPlot(datasetInput(), features = (input$genesearch2))
  })
output$DotPlot <- renderPlot({
  DotPlot(datasetInput(), features = (input$genesearch3))
}) 
  output$dimPlot1 <- renderPlot({
    DimPlot(datasetInput(), reduction = "umap")
  })
  
  output$umapPlot <- renderPlot({
    FeaturePlot(datasetInput(), features = (input$genesearch1), reduction = "umap")
  })
  
  #Select dataset
  datasetInput <- reactive({
    df <- file_contents[[input$dataset]]
  })
  
  #Gene selection
  updateSelectizeInput(session, 'genesearch1', choices = genenames$V1, server = TRUE)
  updateSelectizeInput(session, 'genesearch2', choices = genenames$V1, server = TRUE)
  updateSelectizeInput(session, 'genesearch3', choices = genenames$V1, server = TRUE)
    #Set mandatory fields.
  observe({
    mandatoryFilled <-
      vapply(fieldsMandatory,
             function(x) {
               !is.null(input[[x]]) && input[[x]] != ""
             },
             logical(1))
    
    mandatoryFilled <- all(mandatoryFilled)
    
    shinyjs::toggleState(id = "submit", condition = mandatoryFilled)
  })  
    #Input form data
  formData <- reactive({
    data <- sapply(fieldsAll, function(x) input[[x]])
    data <- c(data, timestamp = epochTime())
    data <- t(data)
    data
  })
  
  #Save form data.
  saveData <- function(data) {
    fileName <- sprintf("%s_%s.csv",
                        humanTime(),
                        digest::digest(data))
    #Write form data into a .csv file.
    write.csv(x = data, file = file.path("Responses", fileName),
              row.names = FALSE, quote = TRUE)
  }
  # action to take when submit button is pressed
  observeEvent(input$submit, {
    shinyjs::disable("submit")
    shinyjs::show("submit_msg")
    shinyjs::hide("error")
    
    tryCatch({
      saveData(formData())
      shinyjs::reset("form")
      shinyjs::hide("form")
      shinyjs::show("thankyou_msg")
    },
    error = function(err) {
      shinyjs::html("error_msg", err$message)
      shinyjs::show(id = "error", anim = TRUE, animType = "fade")
    },
    finally = {
      shinyjs::enable("submit")
      shinyjs::hide("submit_msg")
    })
  })
  #Load form data
  loadData <- function() {
    files <- list.files(file.path("Responses"), full.names = TRUE)
    data <- lapply(files, read.csv, stringsAsFactors = FALSE)
    data <- do.call(rbind, data)
    data
  }
  #Output form responses.
  output$responsesTable <- DT::renderDataTable(
    loadData(),
    rownames = FALSE,
    options = list(searching = FALSE, lengthChange = FALSE)
  ) 
  #Output table
  output$responsesTable1 <- DT::renderDataTable(
    loadData(),
    rownames = FALSE,
    options = list(searching = FALSE, lengthChange = FALSE)
  ) 
}

shinyApp(ui, server)
