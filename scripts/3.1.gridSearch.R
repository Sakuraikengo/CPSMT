##############################################
# Title : 3.1.gridSearch
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
wCombMat <- expand.grid(
  wInit = c(0.9, 0.5, 0.1),
  wFinal = c(0.9, 0.5, 0.1)
)
wLabel <- paste0(wCombMat$wInit, "_", wCombMat$wFinal)

typeCorList <- c("NoRelation", rep(c("Pleiotropy", "SpuriousPleiotropy"), each = 3))
corFactor <- c(
  "NoRelation",
  rep(c("Positive", "NonLinear1_1", "NonLinear1_2"), 2)
)

wColor <- c(
  "#FF0000",
  "#FF5500",
  "#FFC55F",
  "#00CC00",
  "#61C5FF",
  "#0055FF",
  "#800080",
  "#FF69B4",
  "#999999"
)

target1 <- "F8"
baseFileName1 <- "resultMat.csv"
index1 <- c(
  "effectiveProp",
  "effectiveTop",
  "nFail",
  "gainSum",
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
  "effectiveProp"
)
generationInd2 <- paste0("C", formatC(0:nGeneration, width = 2, flag = "0"))

# 1.2. Save dir
dirSave <- "midstream/3.1.gridSearch/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave, recursive = TRUE)
}

# 1.3. Read data
dirBasePath <- "midstream/3.0.CPSMT/"
mcPathList <- paste0(
  dirBasePath,
  rep(wLabel, length(typeCorList)),
  "/",
  rep(typeCorList, each = length(wLabel)),
  "/",
  rep(corFactor, each = length(wLabel)),
  "/"
)

###### 2. Analysis ######
# Read results and build data frame
resList <- lapply(mcPathList, function(path) {
  # path <- mcPathList[[1]]
  pathParts <- str_split(path, pattern = "/")[[1]]
  wEach <- pathParts[3]
  typeFactorEach <- pathParts[4]
  corFactorEach <- pathParts[5]

  # sort the file based on file name
  fileName <- list.files(path, pattern = baseFileName1)
  if (length(fileName) == 0) {
    return(NULL)
  }
  fileNum <- as.numeric(str_sub(
    fileName,
    start = 1,
    end = -nchar(fileName[1])
  ))
  fileOrder <- order(fileNum)

  resFileList <- list.files(path, pattern = baseFileName1, full.names = TRUE)[
    fileOrder
  ]
  resArray <- array(
    data = NA,
    dim = c(length(generationInd1), length(index1), nRep)
  )
  dimnames(resArray) <- list(generationInd1, index1, 1:nRep)

  for (i in 1:nRep) {
    # i <- 1
    resMat0 <- as.matrix(read.csv(resFileList[[i]], row.names = 1))
    for (k in 1:length(index1)) {
      # k <- 1
      traitName <- index1[k]
      if (traitName == "effectiveTop") {
        resArray[, traitName, i] <- (resMat0[, traitName] -
          resMat0[1, traitName]) /
          sqrt(resMat0[1, "geneticVar1"])
      } else if (traitName == "gainSum") {
        resArray[, traitName, i] <- sum(resArray[, "effectiveTop", i])
      } else if (traitName != "nFail") {
        resArray[, traitName, i] <- resMat0[, traitName]
      }
    }
  }

  nFail <- sum(apply(resArray[, "effectiveTop", ], 2, sum) < 0)
  resArray[resArray < -1e6] <- NA

  resMeanMat <- apply(resArray, c(1, 2), mean, na.rm = TRUE)
  resMeanMat[, "nFail"] <- nFail
  resSeMat <- apply(resArray, c(1, 2), std.error, na.rm = TRUE)

  generation <- as.numeric(str_sub(rownames(resMeanMat), start = 2))
  generationRep <- rep(generation, ncol(resMeanMat))
  indexRep <- rep(colnames(resMeanMat), each = nrow(resMeanMat))

  dfEach <- data.frame(
    Generation = generationRep,
    Index = indexRep,
    Type = typeFactorEach,
    Factor = corFactorEach,
    W = wEach,
    Value = c(resMeanMat),
    SE = c(resSeMat)
  )
  return(dfEach)
})

resDf <- do.call(rbind, resList)

# 2.1. Plot
for (i in 1:length(index1)) {
  # i <- 1
  resIndEach0 <- index1[i]
  resDfEach <- resDf[resDf$Index == resIndEach0, ]

  if (any(resIndEach0 == "nFail", resIndEach0 == "gainSum")) {
    next
  }

  finiteIdx <- is.finite(resDfEach$Value) & is.finite(resDfEach$SE)
  ylim <- range(
    c(
      resDfEach$Value[finiteIdx] - resDfEach$SE[finiteIdx],
      resDfEach$Value[finiteIdx] + resDfEach$SE[finiteIdx]
    )
  )
  margin <- abs(ylim[2] - ylim[1]) * 0.05
  ylim[1] <- ylim[1] - margin
  ylim[2] <- ylim[2] + margin

  for (k in 1:length(typeCorList)) {
    # k <- 1
    typeIndEach0 <- typeCorList[k]
    corIndEach0 <- corFactor[k]

    resDfSel <- resDfEach[resDfEach$Type == typeIndEach0, ]
    resDfW <- resDfSel[resDfSel$Factor == corIndEach0, ]

    resDfW$W <- factor(resDfW$W, levels = wLabel)

    g <- ggplot(resDfW, aes(x = Generation, y = Value, color = W)) +
      geom_line() +
      scale_color_manual(values = wColor) +
      geom_errorbar(aes(ymin = Value - SE, ymax = Value + SE), width = .3) +
      ylim(ylim) +
      labs(
        title = paste0(typeIndEach0, " / ", corIndEach0),
        y = resIndEach0
      ) +
      theme(
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 15),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        legend.frame = element_blank()
      )

    png(
      paste0(
        dirSave,
        target1,
        "_",
        resIndEach0,
        "_",
        typeIndEach0,
        "_",
        corIndEach0,
        ".png"
      ),
      width = 1080,
      height = 720,
      res = 214
    )
    print(g)
    dev.off()
  }
}

# 2.2. Summary tables
pattern <- paste0(typeCorList, "_", corFactor)
nFailMat <- matrix(NA, nrow = length(wLabel), ncol = length(pattern))
gainMat <- matrix(NA, nrow = length(wLabel), ncol = length(pattern))
rownames(nFailMat) <- rownames(gainMat) <- paste0("w_", wLabel)
colnames(nFailMat) <- colnames(gainMat) <- pattern

resDfGain <- resDf[resDf$Index == "gainSum", ]
resDfFail <- resDf[(resDf$Index == "nFail") & (resDf$Generation == 0), ]
resDfProp <- resDf[resDf$Index == "effectiveProp", ]

for (wEach in wLabel) {
  # wEach <- wLabel[1]
  dfGainW <- resDfGain[resDfGain$W == wEach, ]
  patternInd <- paste0(dfGainW$Type, "_", dfGainW$Factor)
  gain <- tapply(dfGainW[, "Value"], INDEX = patternInd, mean, na.rm = TRUE)
  gainMat[paste0("w_", wEach), names(gain)] <- gain

  resDfFailW <- resDfFail[resDfFail$W == wEach, ]
  patternInd <- paste0(resDfFailW$Type, "_", resDfFailW$Factor)
  nFail <- tapply(resDfFailW[, "Value"], INDEX = patternInd, mean)
  nFailMat[paste0("w_", wEach), names(nFail)] <- nFail
}
write.csv(gainMat, paste0(dirSave, "totalGain.csv"))
write.csv(nFailMat, paste0(dirSave, "numFail.csv"))

# 2.3. effectiveProp tables per scenario
for (k in 1:length(typeCorList)) {
  typeEach <- typeCorList[k]
  corEach <- corFactor[k]
  resDfPropSel <- resDfProp[
    resDfProp$Type == typeEach & resDfProp$Factor == corEach,
  ]
  propMat <- matrix(
    NA,
    nrow = length(wLabel),
    ncol = length(generationInd1)
  )
  propSeMat <- matrix(
    NA,
    nrow = length(wLabel),
    ncol = length(generationInd1)
  )
  rownames(propMat) <- rownames(propSeMat) <- paste0("w_", wLabel)
  colnames(propMat) <- colnames(propSeMat) <- generationInd1
  for (wEach in wLabel) {
    dfW <- resDfPropSel[resDfPropSel$W == wEach, ]
    matchIdx <- match(as.numeric(gsub("C", "", generationInd1)), dfW$Generation)
    propMat[paste0("w_", wEach), ] <- dfW$Value[matchIdx]
    propSeMat[paste0("w_", wEach), ] <- dfW$SE[matchIdx]
  }
  write.csv(
    propMat,
    paste0(dirSave, typeEach, "_", corEach, "_effectiveProp.csv")
  )
  write.csv(
    propSeMat,
    paste0(dirSave, typeEach, "_", corEach, "_effectivePropSE.csv")
  )
}

# 2.4. effectiveProp from RGS (geneticValMat.csv)
resListRGS <- lapply(mcPathList, function(path) {
  # path <- mcPathList[[1]]
  pathParts <- str_split(path, pattern = "/")[[1]]
  wEach <- pathParts[3]
  typeFactorEach <- pathParts[4]
  corFactorEach <- pathParts[5]

  fileName <- list.files(path, pattern = baseFileName2)
  if (length(fileName) == 0) {
    return(NULL)
  }
  fileNum <- as.numeric(str_sub(
    fileName,
    start = 1,
    end = -nchar(fileName[1])
  ))
  fileOrder <- order(fileNum)

  resFileList <- list.files(path, pattern = baseFileName2, full.names = TRUE)[
    fileOrder
  ]
  resArray <- array(
    data = NA,
    dim = c(length(generationInd2), length(index2), nRep)
  )
  dimnames(resArray) <- list(generationInd2, index2, 1:nRep)

  for (i in 1:nRep) {
    # i <- 1
    resMat0 <- as.matrix(read.csv(resFileList[[i]], row.names = 1))
    for (k in 1:length(index2)) {
      traitName <- index2[k]
      resArray[, traitName, i] <- resMat0[, traitName]
    }
  }

  resMeanMat <- apply(resArray, c(1, 2), mean, na.rm = TRUE)
  resSeMat <- apply(resArray, c(1, 2), std.error, na.rm = TRUE)

  generation <- as.numeric(str_sub(rownames(resMeanMat), start = 2))
  generationRep <- rep(generation, ncol(resMeanMat))
  indexRep <- rep(colnames(resMeanMat), each = nrow(resMeanMat))

  dfEach <- data.frame(
    Generation = generationRep,
    Index = indexRep,
    Type = typeFactorEach,
    Factor = corFactorEach,
    W = wEach,
    Value = c(resMeanMat),
    SE = c(resSeMat)
  )
  return(dfEach)
})

resDfRGS <- do.call(rbind, resListRGS)
resDfPropRGS <- resDfRGS[resDfRGS$Index == "effectiveProp", ]

for (k in 1:length(typeCorList)) {
  typeEach <- typeCorList[k]
  corEach <- corFactor[k]
  resDfPropRGSSel <- resDfPropRGS[
    resDfPropRGS$Type == typeEach & resDfPropRGS$Factor == corEach,
  ]
  propMat <- matrix(
    NA,
    nrow = length(wLabel),
    ncol = length(generationInd2)
  )
  propSeMat <- matrix(
    NA,
    nrow = length(wLabel),
    ncol = length(generationInd2)
  )
  rownames(propMat) <- rownames(propSeMat) <- paste0("w_", wLabel)
  colnames(propMat) <- colnames(propSeMat) <- generationInd2
  for (wEach in wLabel) {
    dfW <- resDfPropRGSSel[resDfPropRGSSel$W == wEach, ]
    matchIdx <- match(
      as.numeric(gsub("C", "", generationInd2)),
      dfW$Generation
    )
    propMat[paste0("w_", wEach), ] <- dfW$Value[matchIdx]
    propSeMat[paste0("w_", wEach), ] <- dfW$SE[matchIdx]
  }
  write.csv(
    propMat,
    paste0(dirSave, typeEach, "_", corEach, "_effectivePropRGS.csv")
  )
  write.csv(
    propSeMat,
    paste0(dirSave, typeEach, "_", corEach, "_effectivePropRGSSE.csv")
  )
}

# check the best combination
paraSuccess <- nFailMat == 0
for (i in 1:ncol(gainMat)) {
  # i <- 1
  gainEach <- gainMat[paraSuccess[, i], i, drop = FALSE]
  bestIdx <- which.max(gainEach)
  print(paste0(
    colnames(gainMat)[i],
    ": ",
    rownames(gainEach)[bestIdx],
    " (",
    round(gainEach[bestIdx], 2),
    ")"
  ))
}
