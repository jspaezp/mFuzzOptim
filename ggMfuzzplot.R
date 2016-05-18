# ggMfuzzplot definition
# This function generates a ggplot object with the clusters
mFuzz.ggplot <- function(data, clustering, centre = TRUE) {
    require(ggplot2)
    
    clusterindex <- clustering$cluster
    
    # data frame with Membership values
    memship <- clustering$membership 
    colnames(memship) <- paste("membership", 
        1:(dim(memship)[2]), 
        sep = ("")) 
    
    exp <- exprs(data) %>% 
        data.frame(. , 
            Identifier = rownames(data),
            clusterindex, memship) %>% 
        tbl_df()
    
    # Transform data frame into a ggplot-compatible format
    exp <- exp %>% 
        gather(sample, 
            expression ,
            - Identifier,
            - clusterindex,
            - contains("membership")) %>% 
        # mutate(sample = gsub("\\w+(?=([0-9]+$))", "", sample, perl = TRUE))
        #  this regular expression deletes all characters and numbers prior to 
        #  the last number in the string z.b. AA00AA00__00 -> 00
        separate(col = sample,  #actually ... this should not be here ... site it applies only to the yeasts data set
            sep = "_", 
            into = c("Something","Time")) %>% 
        select(-Something) %>% 
        mutate(Time = as.numeric(Time))
    
    exp[["maxMembership"]] <- exp %>%  
        select(contains("membership")) %>%
        apply(., 1, max) 
    
    g <- ggplot(data = exp,
            aes(x = Time, y = expression)) +
        geom_line(
            aes(group = Identifier,  
                colour = maxMembership, 
                order = rank(maxMembership)
            )
        ) + 
        scale_colour_gradientn(colours = rainbow(4, alpha = 0.1))
    # Addition of center plotting when centre == TRUE
    if (centre) {
        centers <- clustering$centers %>% 
            data.frame(., 
                clusterindex = rownames(.)) %>% 
            gather(sample, 
                Centre,
                - clusterindex) %>% 
            separate(col = sample, 
                sep = "_",
                into = c("Something","Time")) %>% 
            select(-Something) %>% 
            mutate(Time = as.numeric(Time)) 
        g <- g + geom_line(data = centers, aes(x = Time, y = Centre))
    }
    
    g <- g + facet_wrap(~clusterindex)
    return(g)
} 
# TODO make argument that queries if there is time values in the column names
# TODO Q:I like to display the number of genes in the cluster plots. Can I do this? 
    #  A: There is no functionality within the Mfuzz package for this task (yet). But you can plot the numbers (or any other information) directly in the plots, using the text function of the R graphics package. 


