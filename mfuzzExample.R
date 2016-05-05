# This is a personal scratchsheet

require(Mfuzz)
require(dplyr)

## Data pre-processing test
##
exprSet <- filter.NA(yeast, thres = 0.25)
exprSet <- fill.NA(exprSet, mode = "knn")
exprSet <- standardise(exprSet)

source("app.R/app.R")
mFuzzOpt(data = exprSet)

## MANN lab data
##
##
phosProtData <- read.delim("data/2004856_TableS2.csv",
                                na.strings = c("\\#NUM!","NaN","Nan")
                                )
View(phosProtData)

phosProtDataraw <- phosProtData %>% 
    select(contains("tgf")) %>% 
    select(-contains("median"))

phosProtDataMedian <- phosProtData %>% 
    select(matches("tgf|Gene")) %>% 
    select(matches("median|Gene"))
View(phosProtDataMedian)

rownames(phosProtDataMedian) <- phosProtDataMedian$Gene.names
phosProtDataMedian <- phosProtDataMedian %>% 
    select(-contains("Gene"))

View(phosProtDataMedian)

# TODO make  expression annotations
# 
#exprSet <- ExpressionSet(data.matrix(phosProtDataMedian))

exprSet <- filter.NA(exprSet, thres = 0.25)
exprSet <- fill.NA(exprSet, mode = "knn")
exprSet <- standardise(exprSet)


## DONE

#exprSet <- ExpressionSet(assayData = data.matrix(phosProtDataMedian), exprs = M)

clustering <- mfuzz(exprSet, centers = 6, dist = "euclidean", m = 1.4)
mfuzz.plot2(exprSet, cl=clustering,mfrow=c(2,2),centre=TRUE) # lines for cluster centres will be included


############# GGPLOT TEST
# attempt to ggplot mfuzz.plot
# 

(eset, cl, mfrow = c(1, 1), colo, min.mem = 0, time.labels, 
time.points, ylim = c(0, 0), xlab = "Time", ylab = "Expression changes", 
x11 = TRUE, ax.col = "black", bg = "white", col.axis = "black", 
col.lab = "black", col.main = "black", col.sub = "black", 
col = "black", centre = FALSE, centre.col = "black", centre.lwd = 2, 
Xwidth = 5, Xheight = 5, single = FALSE, ...) 

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

expYeast <- expYeast %>% gather(sample, expression ,
                                - rownames.yeast.,
                                - contains("membership"),
                                - clusterindex) %>% 
    separate(col = sample, sep = "_", into = c("Something","time")) %>% 
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

## TODO loop the generation of cluster files, use grid package

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
