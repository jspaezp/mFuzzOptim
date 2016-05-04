# This is a personal scratchsheet

require(Mfuzz)
require(dplyr)

phosProtData <- read.delim("~/Documents/myFisrtShiny/data/2004856_TableS2.csv",
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

##
##
exprSet <- filter.NA(yeast, thres = 0.25)
exprSet <- fill.NA(exprSet, mode = "knn")
exprSet <- standardise(exprSet)
