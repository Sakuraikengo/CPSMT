##############################################
# Title : 4.1.gridSearch
# Author : Kengo Sakurai
# Date : 2023-08-16
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(plotrix)
library(stringr)
library(ggplot2)
source(file = "scripts/1.0.function.R")

# 1.1. Parameters
nGeneration <- 20 # the number of generations
nRep <- 30 # the number of replication
methodFactor <- c("PCS")
lambdaList <- c(0.001, 0.01, 0.1, 1, 10, 100)
typeCorList <- c("NoRelation", rep(c("Pleiotropy", "SpuriousPleiotropy"), each = 3))
corFactor <- c("NoRelation", rep(c("Positive", "NonLinear1_1", "NonLinear1_2"), 2))
lambdaColor <- c("#0055FF", "#61C5FF", "#BFFAFF", "#FFC55F", "#FF5500", "#FF0000")
target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c("effectiveTop", "nFail", "gainSum", 
            "geneticVar1", "geneticVar2")
generationInd1 <- paste0("C", formatC(seq(0, nGeneration, 2), width = 2, flag = "0"))

# 1.2. Save dir
dirSave <- "midstream/4.1.gridSearch/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

# 1.3. Read data
dirLambdaPath <- "midstream/4.0.PCS/"
lambdaPathList <- paste0(dirLambdaPath, 
                         rep(lambdaList, length(typeCorList)), "/", 
                         rep(typeCorList, each = length(lambdaList)), "/", 
                         rep(corFactor, each = length(lambdaList)), "/")

pathList <- list(list(method = "PCS", path = lambdaPathList))


###### 2. Analysis ######
SimSumLambda(pathList = pathList, 
             target = target1, 
             baseFileName = baseFileName1, 
             index = index1, 
             generationInd = generationInd1, 
             typeFactor = typeCorList, 
             corFactor = corFactor, 
             lambdaList = lambdaList, 
             lambdaColor = lambdaColor)

