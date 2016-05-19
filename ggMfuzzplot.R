# ggMfuzzplot definition
# This function generates a ggplot object with the clusters

mFuzz.ggplot <- function(data, clustering,
                         centre = TRUE, sort.collumns = TRUE) {
    require(ggplot2)
    require(dplyr)
    require(tidyr)
    
    clusterindex <- clustering$cluster
    
    # data frame with Membership values
    memship <- clustering$membership 
    colnames(memship) <- paste("membership", 
        seq_along(memship[1,]), 
        sep = ("")) 
    
    exp <- exprs(data)
    
    # This chunk replaces col names by numbers if 
    # more than 1 is character only 
    # or when sort.collumns is FALSE
    
    all.char.cols <- !grepl("\\d", colnames(exp))
    if ((sum(all.char.cols) > 1) | !sort.collumns) {
        colnames(exp) <- seq_along(all.char.cols)    
    }
    
    exp <- data.frame(exp , 
            Identifier = rownames(data),
            clusterindex, memship) 
    
    # Transform data frame into a ggplot-compatible format
    exp <- exp %>% 
        gather(sample, 
            expression ,
            - Identifier,
            - clusterindex,
            - contains("membership")) %>% 
        mutate(Time = gsub("(\\w*\\D+(?=([0-9]+)))|((?<=\\d)\\D+$)", 
                "", 
                sample,
                perl = TRUE)) %>%
        #  this regular expression deletes all characters and numbers prior to 
        #  the last number in the string z.b. AA00AA00__00 -> 00 else keeps the string
        mutate(Time = gsub("^\\D*$", # this needs to be fixed, bug when seveal character cols ...
                "0", 
                Time,
                perl = TRUE)) %>%
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
    
    # Center plotting when centre == TRUE
    if (centre) {
        centers <- clustering$centers %>% 
            data.frame(., 
                clusterindex = rownames(.)) %>% 
            gather(sample, 
                Centre,
                    - clusterindex) %>% 
            mutate(Time = gsub("(\\w*\\D+(?=([0-9]+)))|((?<=\\d)\\D+$)", 
                    "", 
                    sample,
                    perl = TRUE)) %>%
            #  this regular expression deletes all characters and numbers prior to 
            #  the last number in the string z.b. AA00AA00__00 -> 00 else keeps the string
            mutate(Time = gsub("^\\D*$", # this needs to be fixed, bug when all character names
                    "0", 
                    Time,
                    perl = TRUE)) %>%
            mutate(Time = as.numeric(Time))
        g <- g + geom_line(data = centers, aes(x = Time, y = Centre))
    }
    
    g <- g + facet_wrap(~clusterindex)
    return(g)
} 
# TODO Q:I like to display the number of genes in the cluster plots. Can I do this? 
    #  A: There is no functionality within the Mfuzz package for this task (yet). But you can plot the numbers (or any other information) directly in the plots, using the text function of the R graphics package. 
