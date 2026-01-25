##############################################
# Title : 5.0.comparison
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
methodFactor <- c("CPSMT", "PCS")
typeCorList <- c("NoRelation", rep(c("Pleiotropy", "SpuriousPleiotropy"), each = 3))
corFactor <- c("NoRelation", rep(c("Positive", "NonLinear1_1", "NonLinear1_2"), 2))
wGPList <- c(0.1, 1, 0.1, 1, 1, 0.1, 1)
wLPList <- c(1, 100, 0.001, 1, 0.1, 0.1, 100)
corColor <- c("red", "blue")
target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c("effectiveTop", "effectiveTop10Mean", 
            "geneticVar1", "geneticVar2", "geneticCov", "accuracy1", "accuracy2")
generationInd1 <- paste0("C", formatC(seq(0, nGeneration, 2), width = 2, flag = "0"))

target2 <- "RGS"
baseFileName2 <- "geneticValMat.csv"
index2 <- c("effectiveTop", "effectiveTop10Mean", 
            "geneticVar1", "geneticVar2", "geneticCov", "accuracy1", "accuracy2")
generationInd2 <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))

index3 <- c("accuracyVar1", "accuracyVar2", "accuracyCov", 
            "rmseVar1", "rmseVar2", "rmseCov", 
            "corMarker1", "corMarker2", "nFixedAllele")

index4 <- c("geneticCov")

# 1.2. Save dir
dirSave <- "midstream/5.0.comparison/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

# 1.3. Read data
dirCpPath <- "midstream/3.0.CPSMT/"
cpPathList <- paste0(dirCpPath, 
                     wLPList, "/", 
                     typeCorList, "/", 
                     corFactor, "/")

cpPathList3 <- paste0(dirCpPath, 
                      wLPList[c(1, 2, 5)], "/", 
                      typeCorList[c(1, 2, 5)], "/", 
                      corFactor[c(1, 2, 5)], "/")

dirGpPath <- "midstream/4.0.PCS/"
gpPathList <- paste0(dirGpPath, 
                     wGPList, "/", 
                     typeCorList, "/", 
                     corFactor, "/")


pathList <- list(list(method = "CPSMT", path = cpPathList), 
                 list(method = "PCS", path = gpPathList))

pathList3 <- list(list(method = "CPSMT", path = cpPathList3))

###### 2. Analysis ######
SimSum(pathList = pathList, 
        target = target1, 
        baseFileName = baseFileName1, 
        index = index1, 
        generationInd = generationInd1, 
        methodFactor = methodFactor, 
        typeFactor = typeCorList, 
        corFactor = corFactor, 
        corColor = corColor)

SimSum(pathList = pathList, 
        target = target2, 
        baseFileName = baseFileName2, 
        index = index2, 
        generationInd = generationInd2, 
        methodFactor = methodFactor, 
        typeFactor = typeCorList, 
        corFactor = corFactor, 
        corColor = corColor)

SimSum(pathList = pathList3, 
        target = target2, 
        baseFileName = baseFileName2, 
        index = index3, 
        generationInd = generationInd2, 
        methodFactor = methodFactor, 
        typeFactor = typeCorList, 
        corFactor = corFactor, 
        corColor = corColor)

SimSum(pathList = pathList3, 
        target = target1, 
        baseFileName = baseFileName1, 
        index = index4, 
        generationInd = generationInd1, 
        methodFactor = methodFactor, 
        typeFactor = typeCorList, 
        corFactor = corFactor, 
        corColor = corColor)
