require(shiny)
require(miniUI)
require(Mfuzz)
require(dplyr)
require(tidyr)
require(ggplot2)
require(DT)

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
            miniTabPanel("Plot", icon = icon("plot"),
                         miniContentPanel(padding = 0,
                                          plotOutput("clusters", height = "100%")
                         )
            ),
            miniTabPanel("ggPlot", icon = icon("plot"),
                         miniContentPanel(padding = 0,
                            plotOutput("ggClusterPlot", height = "90%")
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
        
        # ggMfuzzplot definition
        mFuzz.plotgg <- function(data, clustering) {
            clusterindex <- clustering$cluster
            memship <- clustering$membership 
            colnames(memship) <- paste("membership", 
                                       1:(dim(memship)[2]), 
                                       sep = ("")) 
            
            exp <- exprs(data) %>% 
                data.frame(. , rownames(data), clusterindex, memship) %>% 
                tbl_df()
            # Flag ()
            print(exp)
            exp <- exp %>% 
                gather(sample, expression ,
                       - contains("rownames"),
                       - contains("membership"),
                       - clusterindex) %>% 
                separate(col = sample, 
                         sep = "_", 
                         into = c("Something","time")) %>% 
                select(-Something) %>% 
                mutate(time = as.numeric(time))
            
            # Flag ()
            print(exp)
            memships <- grep("membership", 
                             unique(colnames(exp)), 
                             value = TRUE)
            
            exp[["maxMembership"]] <- exp %>%  
                select(contains("membership")) %>%
                apply(., 1, max) 
            
            # Flag ()
            print(exp)
            g <- ggplot(exp, aes(x = time, y = expression)) +
                geom_line(aes(group = rownames.data., 
                              colour = maxMembership, 
                              order = rank(maxMembership)
                )
                ) + 
                scale_colour_gradientn(colours = rainbow(3)) +
                facet_wrap(~clusterindex)    
            return(g)
        } 
        
        # Default mfuzz plot 
        
        output$clusters <- renderPlot({
            print("asdasdasdasdasd")
            collumns <- ceiling(sqrt(input$cClusters))
            rows <- ceiling(input$cClusters/collumns)
            # draw the plot with the specified parameters
            mfuzz.plot2(data, 
                        cl=clustering(),
                        x11 = FALSE, 
                        mfrow=c(rows,collumns),
                        centre=TRUE) # lines for cluster centres will be included
        })

        # ggMfuzzplot Rendering
        output$ggClusterPlot <- renderPlot({
            collumns <- ceiling(sqrt(input$cClusters))
            rows <- ceiling(input$cClusters/collumns)
            pl <- mFuzz.plotgg(data = data, clustering = clustering())
            pl
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
# Return parameters and code
# possibility of selecting members
# ggplot the hell out of the plots
# table delimiting the clusters
# possible annotationDB search and graph generation
# document the "package"
# make this become an actual package ...
