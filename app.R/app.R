

mFuzzOpt <- function(data) {
    require(shiny)
    require(miniUI)
    require(Mfuzz)
    require(dplyr)
    require(tidyr)
    require(ggplot2)
    require(DT)   
    ui <- miniPage(
        
        # Application title
        gadgetTitleBar("Mfuzz optimicer"),
        
        miniTabstripPanel(
            
            # TODO data pre-processing .. or not ...
            #    filter missing ... or not ..
            #    fill missing ... or not ...
            #    standarization ... or not ...
            
            # TODO plot options panel
            # TODO : allow choosing other distance measurements
            # Plot Panel
            miniTabPanel("Plot",
                icon = icon("plot"),
                miniContentPanel(padding = 0,
                    fillCol(flex = c(2, 10),
                        fillRow(
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
                                        value = 1.4),
                            height = "10%"
                            ),
                        plotOutput("clusters", 
                            height = "90%")
                    )
                )
            ),
            miniTabPanel("Clustering Table",
                icon = icon("table"),
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
        
        
        # Default mfuzz plot 
        
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
            # memberships <- (clustering())$cluster 
            attributeTable <- do.call(what = input$tableAttribute, args = list(object = data))
            attributeTable <- data.frame(attributeTable)
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
# Return code
# possibility of selecting members
# table delimiting the clusters
# possible annotationDB search and graph generation
# document the "package"
# make this become an actual package ... dependencies n stuff  ...
# choose appropiate and working icons for the tabs
# 
