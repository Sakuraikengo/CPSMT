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
corFactor <- c(
  "NoRelation",
  rep(c("Positive", "NonLinear1_1", "NonLinear1_2"), 2)
)
wGPList <- c(
  "0.9_0.5",
  "0.5_0.1",
  "0.1_0.5",
  "0.1_0.5",
  "0.5_0.1",
  "0.5_0.9",
  "0.9_0.1"
)
wLPList <- c(
  "0.1_0.5",
  "0.1_0.1",
  "0.1_0.1",
  "0.5_0.1",
  "0.5_0.1",
  "0.5_0.9",
  "0.9_0.1"
)
corColor <- c("red", "blue")

# corLty <- c("solid", "dashed", "longdash", "dotted")
# corPch <- c(19, 2, 6, 15)

target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c(
  "effectiveTop",
  "effectiveTop10Mean",
  "geneticVar1",
  "geneticVar2",
  "geneticCov",
  "accuracy1",
  "accuracy2"
)
generationInd1 <- paste0(
  "C",
  formatC(seq(0, nGeneration, 2), width = 2, flag = "0")
)

target2 <- "RGS"
baseFileName2 <- "geneticValMat.csv"
index2 <- c(
  "effectiveTop",
  "effectiveTop10Mean",
  "geneticVar1",
  "geneticVar2",
  "geneticCov",
  "accuracy1",
  "accuracy2"
)
generationInd2 <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))

index3 <- c(
  "accuracyVar1",
  "accuracyVar2",
  "accuracyCov",
  "rmseVar1",
  "rmseVar2",
  "rmseCov",
  "corMarker1",
  "corMarker2",
  "nFixedAllele"
)

index4 <- c("geneticCov")

# 1.2. Save dir
dirSave <- "midstream/5.0.comparison/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

# 1.3. Read data
dirCpPath <- "midstream/3.0.CPSMT/"
cpPathList <- paste0(dirCpPath, wLPList, "/", typeCorList, "/", corFactor, "/")

cpPathList3 <- paste0(
  dirCpPath,
  wLPList[c(1, 2, 5)],
  "/",
  typeCorList[c(1, 2, 5)],
  "/",
  corFactor[c(1, 2, 5)],
  "/"
)

dirGpPath <- "midstream/4.0.PCS/"
gpPathList <- paste0(dirGpPath, wGPList, "/", typeCorList, "/", corFactor, "/")


pathList <- list(
  list(method = "CPSMT", path = cpPathList),
  list(method = "PCS", path = gpPathList)
)

pathList3 <- list(list(method = "CPSMT", path = cpPathList3))

###### 2. Analysis ######
SimSum(
  pathList = pathList,
  target = target1,
  baseFileName = baseFileName1,
  index = index1,
  generationInd = generationInd1,
  methodFactor = methodFactor,
  typeFactor = typeCorList,
  corFactor = corFactor,
  corColor = corColor
)

SimSum(
  pathList = pathList,
  target = target2,
  baseFileName = baseFileName2,
  index = index2,
  generationInd = generationInd2,
  methodFactor = methodFactor,
  typeFactor = typeCorList,
  corFactor = corFactor,
  corColor = corColor
)

SimSum(
  pathList = pathList3,
  target = target2,
  baseFileName = baseFileName2,
  index = index3,
  generationInd = generationInd2,
  methodFactor = methodFactor,
  typeFactor = typeCorList,
  corFactor = corFactor,
  corColor = corColor
)

SimSum(
  pathList = pathList3,
  target = target1,
  baseFileName = baseFileName1,
  index = index4,
  generationInd = generationInd1,
  methodFactor = methodFactor,
  typeFactor = typeCorList,
  corFactor = corFactor,
  corColor = corColor
)

# 2.1. effectiveProp comparison
propDfList <- lapply(1:length(typeCorList), function(k) {
  # k <- 1
  typeEach <- typeCorList[k]
  corEach <- corFactor[k]
  scenarioName <- paste0(typeEach, "_", corEach)

  # read LP (3.1) data
  lpPropPath <- paste0(
    "midstream/3.1.gridSearch/",
    scenarioName,
    "_effectiveProp.csv"
  )
  lpSePath <- paste0(
    "midstream/3.1.gridSearch/",
    scenarioName,
    "_effectivePropSE.csv"
  )
  lpPropMat <- as.matrix(read.csv(lpPropPath, row.names = 1))
  lpSeMat <- as.matrix(read.csv(lpSePath, row.names = 1))
  lpRow <- lpPropMat[paste0("w_", wLPList[k]), ]
  lpSeRow <- lpSeMat[paste0("w_", wLPList[k]), ]
  generation <- as.numeric(gsub("C", "", names(lpRow)))

  # read GP (4.1) data
  gpPropPath <- paste0(
    "midstream/4.1.gridSearch/",
    scenarioName,
    "_effectiveProp.csv"
  )
  gpSePath <- paste0(
    "midstream/4.1.gridSearch/",
    scenarioName,
    "_effectivePropSE.csv"
  )
  gpPropMat <- as.matrix(read.csv(gpPropPath, row.names = 1))
  gpSeMat <- as.matrix(read.csv(gpSePath, row.names = 1))
  gpRow <- gpPropMat[paste0("w_", wGPList[k]), ]
  gpSeRow <- gpSeMat[paste0("w_", wGPList[k]), ]

  dfEach <- data.frame(
    Generation = rep(generation, 2),
    Value = c(lpRow, gpRow),
    SE = c(lpSeRow, gpSeRow),
    Method = rep(methodFactor, each = length(generation)),
    Type = typeEach,
    Factor = corEach
  )
  return(dfEach)
})
propDf <- do.call(rbind, propDfList)

propDf$Method <- factor(propDf$Method, levels = methodFactor)
propDf$Type <- factor(propDf$Type, levels = unique(typeCorList))
propDf$Factor <- factor(propDf$Factor, levels = unique(corFactor))

ylim <- range(propDf$Value[is.finite(propDf$Value)], na.rm = TRUE)
ylim[1] <- ylim[1] - ylim[1] * 0.1
ylim[2] <- ylim[2] + ylim[2] * 0.1

for (l in 1:length(unique(typeCorList))) {
  # l <- 1
  typeEach <- unique(typeCorList)[l]
  propDfSel <- propDf[propDf$Type == typeEach, ]

  g <- ggplot(propDfSel, aes(x = Generation, y = Value, color = Method)) +
    geom_line(size = 1.5) +
    facet_grid(Type ~ Factor) +
    scale_color_manual(values = corColor) +
    geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE), width = .3) +
    ylim(ylim) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.frame = element_blank()
    )

  if (length(unique(propDfSel$Factor)) > 2) {
    png(
      paste0(dirSave, "effectiveProp_", typeEach, ".png"),
      height = 720,
      width = 2880,
      res = 214
    )
  } else {
    png(
      paste0(dirSave, "effectiveProp_", typeEach, ".png"),
      height = 720,
      width = 960,
      res = 214
    )
  }
  print(g)
  dev.off()
}

# 2.2. effectivePropRGS comparison
propRGSDfList <- lapply(1:length(typeCorList), function(k) {
  # k <- 1
  typeEach <- typeCorList[k]
  corEach <- corFactor[k]
  scenarioName <- paste0(typeEach, "_", corEach)

  # read LP (3.1) data
  lpPropPath <- paste0(
    "midstream/3.1.gridSearch/",
    scenarioName,
    "_effectivePropRGS.csv"
  )
  lpSePath <- paste0(
    "midstream/3.1.gridSearch/",
    scenarioName,
    "_effectivePropRGSSE.csv"
  )
  lpPropMat <- as.matrix(read.csv(lpPropPath, row.names = 1))
  lpSeMat <- as.matrix(read.csv(lpSePath, row.names = 1))
  lpRow <- lpPropMat[paste0("w_", wLPList[k]), ]
  lpSeRow <- lpSeMat[paste0("w_", wLPList[k]), ]
  generation <- as.numeric(gsub("C", "", names(lpRow)))

  # read GP (4.1) data
  gpPropPath <- paste0(
    "midstream/4.1.gridSearch/",
    scenarioName,
    "_effectivePropRGS.csv"
  )
  gpSePath <- paste0(
    "midstream/4.1.gridSearch/",
    scenarioName,
    "_effectivePropRGSSE.csv"
  )
  gpPropMat <- as.matrix(read.csv(gpPropPath, row.names = 1))
  gpSeMat <- as.matrix(read.csv(gpSePath, row.names = 1))
  gpRow <- gpPropMat[paste0("w_", wGPList[k]), ]
  gpSeRow <- gpSeMat[paste0("w_", wGPList[k]), ]

  dfEach <- data.frame(
    Generation = rep(generation, 2),
    Value = c(lpRow, gpRow),
    SE = c(lpSeRow, gpSeRow),
    Method = rep(methodFactor, each = length(generation)),
    Type = typeEach,
    Factor = corEach
  )
  return(dfEach)
})
propRGSDf <- do.call(rbind, propRGSDfList)

propRGSDf$Method <- factor(propRGSDf$Method, levels = methodFactor)
propRGSDf$Type <- factor(propRGSDf$Type, levels = unique(typeCorList))
propRGSDf$Factor <- factor(propRGSDf$Factor, levels = unique(corFactor))

ylimRGS <- range(propRGSDf$Value, na.rm = TRUE)
ylimRGS[1] <- ylimRGS[1] - ylimRGS[1] * 0.1
ylimRGS[2] <- ylimRGS[2] + ylimRGS[2] * 0.1

for (l in 1:length(unique(typeCorList))) {
  # l <- 3
  typeEach <- unique(typeCorList)[l]
  propRGSDfSel <- propRGSDf[propRGSDf$Type == typeEach, ]

  g <- ggplot(propRGSDfSel, aes(x = Generation, y = Value, color = Method)) +
    geom_line(size = 1.5) +
    facet_grid(Type ~ Factor) +
    scale_color_manual(values = corColor) +
    geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE), width = .3) +
    ylim(ylimRGS) +
    theme(
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.text.y = element_text(size = 20),
      axis.title.y = element_blank(),
      legend.title = element_blank(),
      legend.frame = element_blank()
    )

  if (length(unique(propRGSDfSel$Factor)) > 2) {
    png(
      paste0(dirSave, "effectivePropRGS_", typeEach, ".png"),
      height = 720,
      width = 2880,
      res = 214
    )
  } else {
    png(
      paste0(dirSave, "effectivePropRGS_", typeEach, ".png"),
      height = 720,
      width = 960,
      res = 214
    )
  }
  print(g)
  dev.off()
}

# compare the result of total gain
pcsAll <- read.csv("midstream/4.1.gridSearch/totalGain.csv", row.names = 1)
pcs <- sapply(1:length(wGPList), function(w) {
  pcsAll[paste0("w_", wGPList[w]), w]
})
cpsAll <- read.csv("midstream/3.1.gridSearch/totalGain.csv", row.names = 1)
cps <- sapply(1:length(wLPList), function(w) {
  cpsAll[paste0("w_", wLPList[w]), w]
})
cps / pcs
mean(cps / pcs)
mean(cps[c(1, 3, 6)] / pcs[c(1, 3, 6)])
mean(cps[c(2, 4, 5, 7)] / pcs[c(2, 4, 5, 7)])

cat(paste0(
  "SpuriousPleiotropy_Positive vs Pleiotropy_Positive:\n",
  "  PCS +",
  round((pcs[5] / pcs[2] - 1) * 100, 0),
  "%\n",
  "  CPS-MT +",
  round((cps[5] / cps[2] - 1) * 100, 0),
  "%\n"
))

cat(paste0(
  "SpuriousPleiotropy_NonLinear1 vs Pleiotropy_NonLinear1:\n",
  "  PCS +",
  round((pcs[6] / pcs[3] - 1) * 100, 0),
  "%\n",
  "  CPS-MT +",
  round((cps[6] / cps[3] - 1) * 100, 0),
  "%\n"
))

cat(paste0(
  "SpuriousPleiotropy_NonLinear2 vs Pleiotropy_NonLinear2:\n",
  "  PCS +",
  round((pcs[7] / pcs[4] - 1) * 100, 0),
  "%\n",
  "  CPS-MT +",
  round((cps[7] / cps[4] - 1) * 100, 0),
  "%\n"
))
