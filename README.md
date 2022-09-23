# Shiny_Lab_Dashboard
How to create a shiny lab dash board.

MetadataSingleCellBrowser
2022-09-22
#About The R Shiny Single Cell App is a Single Cell Browser is an interactive app that allows you to upload and explore previously submitted single cell datasets and their metadata. The app shows a umap, violin plots and gene-expression dot plots as well as a metadata table and a submission form.

#Set up The app is compatible with single cell data analysed in seurat. To work, you will need to run a Umap during your analysis. Save your seurat object as an .rds file using the code below

saveRDS(object = x, file = "yourobjectname.rds")
Create a folder for your app (e.g. “single cell app folder”) and save the shiny script here. Then create two subdirectories - “Responses” and “Database”. The “Responses” directory will contain all of the submitted metadata responses as .csv files. The Database directory will contain the .csv file of genenames for your species of interest (as extracted from your species gtf file). You will need to save your .rds seurat files in the main folder (“single cell app folder”).

There are two parts of the app that need to be adjusted to your preferences. These are the genenames.csv file and the select input in the ui. Instrucutions on how to do this are below.

#Required packages

library(shiny)
library(shinythemes)
library(Seurat)
library(uwot)
library(ggplot2)
library(dplyr)
The below sets the Mandatory fields in the metadata form.

#set mandatory fields
fieldsMandatory <- c("Name","Experiment.Name", "Type", "Genotype", "Stage", "Time","Tissue", "Protocol", "Replicates", "Pre.process.analysis", "Analysis", "RDS.name", "Folder.Name", "Raw.Data.ID")
fieldsAll <- c("Name","Experiment.Name", "Type", "Genotype", "Stage", "Time","Tissue", "Protocol", "Replicates", "Pre.process.analysis", "Analysis", "RDS.name", "Folder.Name", "Raw.Data.ID", "Publications.Document","Additional.Comments")
#Read .csv files There are three parts of the R script that need updating to your preferences.

The first is the genenames.csv file. You will need to create this file which will contain all of the genes annotated in your species of interests genome in a single column (see the drosophila example in the repository). Genenames need to match the format used in Seurat. Save your genenames file in the “Database_files” directory and update the below with your genenames .csv name.

genenames <- read.csv("Database_files/fullgenelist_Drosophila.csv")
#Functions When metadata is saved in the browser, it includes the .rds filename. The .rds files can be read into the app using the loop below. This takes the RDS.name column from responses, removes duplicates and N.A’s and saves as the datanames object. The object is then used in the following loop to read in all .rds files in the single cell app folder.

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
Functions
The functions below relate to submitting new entries into the single cell database. “Responses” is a directory in the SingleCell.RDS folder. epoch and human time are used as filenames for metadata saved in the form so they can be easily deleted/updated if required based on the submission date. Each form submission is an individual .csv file.

Label Mandatory sets the ’*’ icon to show which fields are mandatory fields = if all mandatory fields are not filled in, then the form response will not submit. The row headers for ‘Mandatory fields’ as well as all the form fields are set as objects; these correspond to the form headers set in ‘Submit New Datasets’ - if you decided to change the metadata collected, you need to update Mandatory, all fields and the textInput headers.

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

appCSS <- ".mandatory_star { color: red; }"
humanTime <- function() format(Sys.time(), "%Y%m%d-%H%M%OS")
#ui section

The user interface script is broken down into navigation panels. The first section sets the theme (graphical display used by the app = look at shinythemes options if you wish to change the theme). It also tells the app to use json and CSS which allows you to set mandatory fields and hidden messages in the metadata submission form section.

ui <-fluidPage(
  shinyjs::useShinyjs(),
  shinyjs::inlineCSS(appCSS), navbarPage("Brand Lab Dashboard", theme = shinytheme("cosmo"), 
#Single Cell Protocol This sets the first section of the app where you can provide any analysis, protocol or method notes that would helpful to users. Most of this information is coded using the html tags feature in R (for example, to add a hyperlink = tags$a(href=““).All text needs to be written between”“. Use p() either side of text to foramt paragraphs.

 tabPanel(title = div(icon("chart-line", lib = "font-awesome"), "Single Cell"),
                          navlistPanel(
                            id = "tabset",
                            "Protocols", 
                            tabPanel("Analysis",
                                     
                                     p("INSERT YOUR TEXT HERE"), tags$a(href="https://satijalab.org/seurat/", "here!"),
                                     p("INSERT YOUR TEXT HERE"),
),
#Single Cell This section allows users to submit new datasets using the form. Each response is saved in the “Responses” directory as a .csv file. The .csv responses are “bound” and rendered to a table to show all the datasets. Mandatory fields are stated at the beginig of the script (fieldsMandatory). The action button placed at the end of UI is pressed to submit responses in the form and save these as a .csv - see server section.

"Datasets",
            tabPanel("Single Cell Lab Datasets",
             mainPanel(
              dataTableOutput(outputId = "SingleCell_results"))
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
A “Thank you for submitting” message is displayed when the submit button is pressed and a return to prev page to submit another responses is shown as a link (this uses JSON - make sure shinyjs::useShinyjs() and shinyjs::inlineCSS(appCSS) are loaded at start of the script or this hidden section will not run.

    shinyjs::hidden(
       div(
        id = "thankyou_msg",
        h3("Thanks, your response was submitted successfully!"),
        actionLink("submit_another", "Submit another response")
      )
    ),  
    DT::dataTableOutput("responsesTable")
  ))
This sections shows the Umap for the selected single cell dataset. Users need to ensure that a UMAP has been run in their analysis and that this feature is saved in the .rds object at the end of their seurat analysis script. UPDATE: the selected = “.rds” to your favorite .rds file. This will be automatically selected when the app loads.

           "sc Explorer",
           tabPanel("Umap Browser", 
                    titlePanel("Umap Browser"),
                    fluidRow(
                      column(8,
                             selectInput("dataset", label = h3("Dataset"),
                                         choices = named_file_list,
                                         selected = "st17_3replicates.rds"),
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
This section is to show Violin plots of selected gene names. UPDATE: the selected = “.rds” to your favorite .rds file. This will be automatically selected when the app loads.

           tabPanel("Violin Plots", 
                    titlePanel("Violin Plots"),
                    fluidRow(
                      column(8,
                             selectInput("dataset", label = h3("Dataset"),
                                         choices = named_file_list,
                                         selected = "st17_3replicates.rds"),
                             helpText("Gene names must be exact."),
                             selectizeInput("genesearch2", label = "Gene Name", choices = NULL, multiple = FALSE, options = NULL)
                      )),
                    fluidRow(
                      column(8, 
                             plotOutput("vlnPlot", width = "200%"))),
                    fluidRow(column(4,
                                    downloadButton("downloadvln", "Download.png")))),
This section shows gene expression per cluster as a dot plot.UPDATE: the selected = “.rds” to your favorite .rds file. This will be automatically selected when the app loads.

           tabPanel("Mulit-gene Browser",
                    titlePanel("Mulit-gene Browser"),
                    fluidRow(
                      column(8,
                             selectInput("dataset", label = h3("Dataset"),
                                         choices = named_file_list,
                                         selected = "st17_3replicates.rds"),
                             helpText("Enter Gene names must be exact."),
                             selectizeInput("genesearch3", label = "Gene Name", choices = NULL, multiple = TRUE, options = NULL)
                      )),
                    fluidRow(
                      column(8, 
                             plotOutput("DotPlot", width = "200%"))),
                    fluidRow(column(4,
                                    downloadButton("downloaddot", "Download.png"))))))))
##Server The server section is the ‘backend’ of the ui. thematic::thematic_shiny()’ corresponds with the ui theme selection to set the theme. renderDataTable shows data tables.ll output functions relate to an input e.g output$SingleCell - “SingleCell” is called in the ui.

server <- function(input, output, session) {
    thematic::thematic_shiny()
Load in previous submission data

   loadData <- function() {
    files <- list.files(file.path(responsesDir), full.names = TRUE)
    data <- lapply(files, read.csv, stringsAsFactors = FALSE)
    data <- do.call(rbind, data)
    data
  }
Each output feature called in the ui e.g. plotOutput(“DotPlot”) has corresponding code in the server (e.g. the code to create a plot etc). NOTE: the single cell dataset table is called later in the app.

  datasetInput <- reactive({
    df <- (input$dataset)
  })
  
  output$dimPlot1 <- renderPlot({
    DimPlot(datasetInput(), reduction = "umap")
  })
  
  output$umapPlot <- renderPlot({
FeaturePlot(datasetInput(), features = (input$genesearch1), reduction = "umap")
  })
Select dataset from user input.

  datasetInput <- reactive({
    df <- file_contents[[input$dataset]]
  })
User specifies genenames in the interface. This pulls the genenames from the genenames.csv file.

  updateSelectizeInput(session, 'genesearch1', choices = genenames$V1, server = TRUE)
  updateSelectizeInput(session, 'genesearch2', choices = genenames$V1, server = TRUE)
  updateSelectizeInput(session, 'genesearch3', choices = genenames$V1, server = TRUE)
Plot feature plots from the seurat object using the selected gene names.

  output$vlnPlot <- renderPlot({
      VlnPlot(datasetInput(), features = (input$genesearch2))
  })
  output$DotPlot <- renderPlot({
    DotPlot(datasetInput(), features = (input$genesearch3))
  })
The below ensures that all Mandatory fields are complete before the submit button can be selected. Mandatory fields are indetfied with * and stated at the beginning of the script.

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
The functions below save the metadata submitted in the form as a .csv file named with a timestamp. Each indvidual entry has a .csv file which can be found in the ‘Responses’ directory.

  formData <- reactive({
    data <- sapply(fieldsAll, function(x) input[[x]])
    data <- c(data, timestamp = epochTime())
    data <- t(data)
    data
  })

  saveData <- function(data) {
    fileName <- sprintf("%s_%s.csv",
                        humanTime(),
                        digest::digest(data))
    
    write.csv(x = data, file = file.path(responsesDir, fileName),
              row.names = FALSE, quote = TRUE)
  }
The below code saves the form entry data when the submit button is pressed. It shows the user a thank you message indicating that the metadata has been succesfully saved.

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
Data that has been submitted in the form is saved as individual .csv files in the ‘Responses’ directory. The below binds these .csv files together into a single object- data.

  loadData <- function() {
    files <- list.files(file.path(responsesDir), full.names = TRUE)
    data <- lapply(files, read.csv, stringsAsFactors = FALSE)
    data <- do.call(rbind, data)
    data
  }
This code renders the metadata submitted in the from as an interactive data table. There are two responses tables but these are the same (one is just to browse the metadata) the other is displayed on the form entry page as a reference.

  output$responsesTable <- DT::renderDataTable(
    loadData(),
    rownames = FALSE,
    options = list(searching = FALSE, lengthChange = FALSE)
  ) 
  
  output$responsesTable1 <- DT::renderDataTable(
    loadData(),
    rownames = FALSE,
    options = list(searching = FALSE, lengthChange = FALSE)
  ) 
    }
You are ready to launch the app!

shinyApp(ui, server)
