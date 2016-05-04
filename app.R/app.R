require(Mfuzz)
require(shiny)
require(miniUI)
require(data.table)

mFuzzOpt <- function(data) {
    
    ui <- miniPage(
        
        # Application title
        gadgetTitleBar("Mfuzz optimicer"),
        
        miniTabstripPanel(
            # Sidebar with a slider input for number of clusters and fuzziness
            # 
            miniTabPanel("Parameters", icon = icon("sliders"),
                         miniContentPanel(
                             sliderInput("cClusters",
                                         "Number of Clusters",
                                         step = 1,
                                         min = 2,
                                         max = 20,
                                         value = 4),
                             sliderInput("fuzzyness",
                                         "Fuzziness",
                                         step = 0.1,
                                         min = 1,
                                         max = 5,
                                         value = 1.4)
                         )
            ),
            # TODO data pre-processing .. or not ...
            #    filter missing ... or not ..
            #    fill missing ... or not ...
            #    standarization ... or not ...
            
            # TODO plot options panel
            # TODO : allow choosing other distance measurements
            # Plot Panel
            miniTabPanel("Plot", icon = icon("clusterPlot"),
                         miniContentPanel(padding = 0,
                                          plotOutput("clusters", height = "100%")
                         )
            ),
            
            miniTabPanel("Clustering Table", icon = icon("table"),
                         miniContentPanel(
                             selectInput("tableAttribute",
                                         "Desired Table:",
                                         c("Assay Data" = "assayData",
                                           "phenoData" = "phenoData",
                                           "featureData" = "featureData",
                                           "experimentData" = "experimentData",
                                           "annotation" = "annotation",
                                           "protocolData" = "protocolData",
                                           "Expression" = "exprs")
                             ),
                             
                             DT::dataTableOutput("table") 
                         )
            )
        )
    )
    
    # Define server logic 
    server <- function(input, output) {
        
        if (!class(data) == "ExpressionSet") {
            message("Data structure not supported, please use an ExpressionSet")
        }
        
        clustering <- reactive(
            mfuzz(data, 
                  centers = input$cClusters,
                  dist = "euclidean",
                  m = input$fuzzyness)
        ) 
        
        output$clusters <- renderPlot({
            collumns <- ceiling(sqrt(input$cClusters))
            rows <- ceiling(input$cClusters/collumns)
            # draw the plot with the specified parameters
            mfuzz.plot2(data, 
                        cl=clustering(),
                        x11 = FALSE, 
                        mfrow=c(rows,collumns),
                        centre=TRUE) # lines for cluster centres will be included
        })
        
        # output table of the clustering
        
        output$table <- DT::renderDataTable({
            
            memberships <- (clustering())$cluster 
            attributeTable <- data.table(
                do.call(what = input$tableAttribute, args = list(object = data)),
                keep.rownames = TRUE
            )
            attributeTable
        })
        
        # Handle the Done button being pressed.
        observeEvent(input$done, {
            stopApp( returnValue = list("fuzzyness" = input$cClusters ,
                                        "clusters" = input$fuzzyness)
            )
        })
    }
    
    # Run the application 
    runGadget(shinyApp(ui = ui, server = server), 
              viewer = browserViewer(browser = getOption("browser"))
    )
    
}

# TODO 
# Return parameters and code
# possibility of selecting members
# ggplot the hell out of the plots
# table delimiting the clusters
# possible annotationDB search and graph generation
# document the "package"
# make this become an actual package ...
