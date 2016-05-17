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
    # Transform data frame into a ggplot-compatible format
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
    
    # memships <- grep("membership", 
    #     unique(colnames(exp)), 
    #     value = TRUE)
    
    exp[["maxMembership"]] <- exp %>%  
        select(contains("membership")) %>%
        apply(., 1, max) 
    
    g <- ggplot(data = exp,
            aes(x = time, y = expression)) +
        geom_line(
            aes(group = rownames.data., 
                colour = maxMembership, 
                order = rank(maxMembership)
            )
        ) + 
        scale_colour_gradientn(colours = rainbow(3)) +
        facet_wrap(~clusterindex)    
    return(g)
} 
