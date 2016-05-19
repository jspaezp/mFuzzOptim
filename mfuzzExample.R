# This is a personal scratchsheet

source("https://bioconductor.org/biocLite.R")
biocLite("Mfuzz")

require(Mfuzz)
data(yeast)
## Data pre-processing test
exprSet <- filter.NA(yeast, thres = 0.25)
exprSet <- fill.NA(exprSet, mode = "knn")
exprSet <- standardise(exprSet)

source("app.R/app.R")
mFuzzOpt(data = exprSet)

clustering <- mfuzz(exprSet, 4, 1.2)
mfuzz.plot2(exprSet, clustering)
source("ggMfuzzplot.R")
mFuzz.ggplot(exprSet, clustering)

## MANN lab data
##
##
require(Mfuzz)
require(dplyr)
phosProtData <- read.delim("data/2004856_TableS2.csv",
    na.strings = c("\\#NUM!","NaN","Nan")
    )
# View(phosProtData)

phosProtDataraw <- phosProtData %>% 
    select(contains("tgf")) %>% 
    select(-contains("median"))

phosProtDataMedian <- phosProtData %>% 
    select(matches("tgf|Gene")) %>% 
    select(matches("median|Gene"))

# rownames(phosProtDataMedian) <- phosProtDataMedian$Gene.names
phosProtDataMedian <- phosProtDataMedian %>% 
    select(-contains("Gene"))


# DONE make  expression annotations
# gorgeourly done by other group 

ppexprSet <- ExpressionSet(data.matrix(phosProtDataMedian))
ppexprSet <- filter.NA(ppexprSet, thres = 0.25)
ppexprSet <- fill.NA(ppexprSet, mode = "knn")
ppexprSet <- standardise(ppexprSet)

mFuzzOpt(ppexprSet)
## DONE

#exprSet <- ExpressionSet(assayData = data.matrix(phosProtDataMedian), exprs = M)

ppclustering <- mfuzz(ppexprSet, centers = 6, dist = "euclidean", m = 1.4)
mfuzz.plot2(ppexprSet, cl=ppclustering,mfrow=c(2,2),centre=TRUE) # lines for cluster centres will be included

mFuzz.ggplot(ppexprSet, ppclustering)
data <- ppexprSet
clustering <- ppclustering

############# GGPLOT TEST
# attempt to ggplot mfuzz.plot
# 

#Loading a test data set, normalizing and filling NAs
data(yeast)

yeast <- filter.NA(yeast, thres = 0.25)
yeast <- fill.NA(yeast, mode = "knn")
yeast <- standardise(yeast)

exprSet <- yeast 
cl <- mfuzz(yeast, centers = 4, m = 1.2)

#declaring arbitrary minimum membership
min.mem = 0

clusterindex <- cl$cluster
memship <- cl$membership 
memship[memship < min.mem] <- -1
colnames(memship) <- paste("membership", 1:(dim(memship)[2]), sep = ("")) 

require(tidyr)

expYeast <- exprs(yeast) %>% 
    data.frame(. , rownames(yeast), clusterindex, memship) %>% 
    tbl_df()

expYeast <- expYeast %>% 
    gather(sample, expression ,
        - rownames.yeast.,
        - contains("membership"),
        - clusterindex) %>% 
    separate(col = sample,
        sep = "_", 
        into = c("Something","time")) %>% 
    select(-Something) %>% 
    mutate(time = as.numeric(time))


g <- ggplot(expYeast, aes(x = time, y = expression)) +
    geom_line(aes(group = rownames.yeast., colour = membership1)) +
    facet_wrap(~clusterindex)

#This one workd gorgeously
memships <- grep("membership", unique(colnames(expYeast)), value = TRUE)
expYeast[["maxMembership"]] <- expYeast %>%  
    select(contains("membership")) %>%
    apply(., 1, max) 

g <- ggplot(expYeast, aes(x = time, y = expression)) +
    geom_line(aes(group = rownames.yeast., 
                  colour = maxMembership, 
                  order = rank(maxMembership)
    )
    ) + 
    scale_colour_gradientn(colours = rainbow(7)) +
    facet_wrap(~clusterindex)
g

g <- ggplot(expYeast, aes(x = time, y = expression)) 
memships <- grep("membership", unique(colnames(expYeast)), value = TRUE)


subsetMembership <- function(data, min.mem, col) {
    data[data[col] > min.mem,] %>% 
        select_("rownames.yeast.", "time", col, "expression") %>%
        mutate_("rank" = paste("rank(", m, ")", sep = "")) #hacky fix to iterate over the variables
}

ggTemplate <- 
    gglist <- list()

for (m in memships) {
    gglist[[m]] <- 
        ggplot(data = subsetMembership(expYeast, 
                                       min.mem = min.mem,
                                       col = m),
               aes_string(x = "time", y = "expression")) +
        geom_line(aes_string(group = "rownames.yeast.",
                             colour = m,
                             order = "rank"))
}

multiplot(plotlist = gglist)


############-------------- Time measurements ----------------#################

require(shiny)
require(miniUI)
require(Mfuzz)
require(dplyr)
require(tidyr)
require(ggplot2)
require(DT)
data("yeast")

# Data pre-processing
exprSet <- filter.NA(yeast, thres = 0.25)
exprSet <- fill.NA(exprSet, mode = "knn")
exprSet <- standardise(exprSet)

clustering <- mfuzz(exprSet, 4, 1.15)

flag <- Sys.time()

mfuzz.plot2(exprSet, cl=clustering,mfrow=c(2,2),centre=TRUE) # lines for cluster centres will be included

time <- Sys.time() - flag
print(time)


flag <- Sys.time()

pl <- mFuzz.ggplot(data = exprSet, clustering = clustering)
pl

time <- Sys.time() - flag
print(time)


#################### --------- Benchmarking ------------ ####################
require(microbenchmark)    

    f1 <- function() {
        1:(dim(memship)[2])
    }
    
    f2 <- function() {
        seq_along(memship[1,])
    }
    
compare <- microbenchmark(f1(), f2(), times = 10000)

library(ggplot2)
autoplot(compare)

