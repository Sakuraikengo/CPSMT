##############################################
# Title : 7.0.CPSMT3
# Author : Kengo Sakurai
# Date : 2026-05-23
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(breedSimulatR)
library(stringr)
library(RAINBOWR)
library(MASS)
library(lpSolve)
library(abind)
source(file = "scripts/1.0.function.R")

# 1.1. Parameters (breeding simulation)
nTrait <- 3
nProg <- 15
nProgSelf <- 50
probSelf <- 0.9
probCross <- 0.99333333
capa <- 2
nCrosses <- 10
nSelf <- 5
nGeneration <- 20
nPhenoField <- 100
nMCCross <- 15000
nMCSelf <- 10000

# Single w combination (demo)
wInit <- 0.9
wFinal <- 0.1

# Bounds for essential traits (length = nTrait - 1)
# trait2: [-1.08, 1.08], trait3: [-0.5, 0.5]
lVec <- c(-1.08, -0.5)
hVec <- c(1.08, 0.5)

f <- 8 # generation for evaluation (F8)

seed <- 123

###### 2. Load pre-built setting data (from 6.0.simulationSettingMT.R) ######
i <- 1
gCorName <- "NoRelation"
coeff <- NA
dirSetting <- paste0("midstream/6.0.simulationSettingMT/", gCorName, "/")

initPop <- readRDS(paste0(dirSetting, i, "_initPop.rds"))
map <- read.csv(paste0(dirSetting, i, "_map.csv"), row.names = 1)
initRes <- readRDS(paste0(dirSetting, i, "_initResList.rds"))
trueU <- initRes$uMat
pheno <- initRes$pheno
genoMat <- initRes$genoMat
sigmaE <- initRes$sigmaE

markerEffectTrueMat <- as.matrix(
  read.csv(paste0(dirSetting, i, "_markerEffects.csv"), row.names = 1)
)

###### 3. Output directory ######
dirSave <- paste0(
  "midstream/7.0.CPSMT3/",
  wInit,
  "_",
  wFinal,
  "/",
  gCorName,
  "/"
)
if (!dir.exists(dirSave)) {
  dir.create(dirSave, recursive = TRUE)
}

# Indices for variance / covariance columns of sigmaMat (CalcVarMat ordering)
covPairsIdx <- combn(nTrait, 2)
covLabels <- apply(covPairsIdx, 2, paste, collapse = "")
nVarCovCol <- nTrait + ncol(covPairsIdx)

###### 4. Recurrent genomic selection with multi-trait CPS-MT (in-memory) ######
set.seed(seed)

# resHistory[[g+1]] holds the per-generation state (g = 0, 1, ..., nGeneration)
# f8History[[k]] holds the F8 evaluation for even generations
resHistory <- vector("list", nGeneration + 1)
f8History <- list()

generation <- 0

model <- GP(phenoMat = pheno, genoMat = genoMat - 1)
predU <- model$uPredMat
markerEffectEstimatedMat <- model$markerEffPredMat
beta <- model$betaEstimated

# Population statistics at generation 0
AlleleFreq <- apply(genoMat, 2, function(eachAllele) {
  p <- sum(eachAllele) / (2 * length(eachAllele))
  return(4 * p * (1 - p))
})
genicVar <- apply(markerEffectTrueMat, 2, function(x) AlleleFreq %*% x^2)
geneticVar <- apply(trueU, 2, function(x) var(x) * (length(x) - 1) / length(x))
geneticCov <- sapply(1:ncol(covPairsIdx), function(cIdx) {
  t1 <- covPairsIdx[1, cIdx]
  t2 <- covPairsIdx[2, cIdx]
  var(trueU[, t1], trueU[, t2]) * (nrow(trueU) - 1) / nrow(trueU)
})
bulmer <- geneticVar / genicVar

D2 <- CalcD2(
  ObjectMap = map,
  markerEffectMat = markerEffectEstimatedMat,
  nTrait = nTrait,
  k = f - 1
)
D2True <- CalcD2(
  ObjectMap = map,
  markerEffectMat = markerEffectTrueMat,
  nTrait = nTrait,
  k = f - 1
)
D2Self <- CalcD2(
  ObjectMap = map,
  markerEffectMat = markerEffectEstimatedMat,
  nTrait = nTrait,
  k = f - 2
)

predProgMean <- CalcProgMean(predU = predU)
indName <- rownames(predU)
comb <- combn(indName, 2)
combInd <- combn(1:length(indName), 2)
combName <- apply(comb, 2, function(x) paste0(x[1], "_", x[2]))

gametArray <- ExtractGamet(
  nowPop = initPop,
  indName = indName,
  qtn = rownames(map)
)
dimnames(gametArray)[[2]] <- colnames(genoMat)

sigmaMat <- CalcVarMat(
  D2 = D2,
  combInd = combInd,
  gametArray = gametArray,
  nTrait = nTrait
)
rownames(sigmaMat) <- rownames(predProgMean) <- combName
sigmaTrueMat <- CalcVarMat(
  D2 = D2True,
  combInd = combInd,
  gametArray = gametArray,
  nTrait = nTrait
)
rownames(sigmaTrueMat) <- combName

w <- wInit + (wFinal - wInit) * generation / nGeneration
meanPredU <- apply(predU, 2, mean)
sdPredU <- apply(predU, 2, sd)

mcResult <- EvalByMC(
  meanMat = predProgMean,
  varMat = sigmaMat,
  w = w,
  lVec = lVec,
  hVec = hVec,
  prob = probCross,
  nMC = nMCCross,
  meanPredU = meanPredU,
  sdPredU = sdPredU,
  nTrait = nTrait
)
predProgMat <- matrix(mcResult, ncol = 1)
rownames(predProgMat) <- rownames(predProgMean)
colnames(predProgMat) <- "MC"

crossTable <- LP(
  predProgMat = predProgMat,
  indName = indName,
  comb = comb,
  capa = capa,
  nCrosses = nCrosses,
  nProg = 1,
  generation = generation
)
selectedCross <- unique(c(crossTable$ind1, crossTable$ind2))

resHistory[[generation + 1]] <- list(
  pheno = pheno,
  genoMat = genoMat,
  predU = predU,
  trueU = trueU,
  markerEffectEstimatedMat = markerEffectEstimatedMat,
  beta = beta,
  sigmaMat = sigmaMat,
  sigmaTrueMat = sigmaTrueMat,
  alleleFreq = AlleleFreq,
  genicVar = genicVar,
  geneticVar = geneticVar,
  geneticCov = geneticCov,
  bulmer = bulmer,
  selectedCross = selectedCross,
  selectedSelf = NULL,
  w = w,
  generation = generation
)

# F8 entry for generation 0 (no selfing decision yet; use base population stats)
f8History[[as.character(generation)]] <- list(
  pheno = pheno,
  genoMat = genoMat,
  predU = predU,
  trueU = trueU,
  genicVar = genicVar,
  geneticVar = geneticVar,
  geneticCov = geneticCov,
  bulmer = bulmer,
  generation = generation
)

nextPop <- population$new(
  name = "C0 offspring",
  inds = makeCrosses(crosses = crossTable, pop = initPop)
)

for (generation in 1:nGeneration) {
  nowPop <- nextPop
  genoMatNow <- nowPop$genoMat

  # Update GP every 2 generations after gen 3, using accumulated past data
  if ((generation > 3) & (generation %% 2 == 0)) {
    resPast <- resHistory[[generation]] # generation - 1, index +1
    phenoNewAll <- resPast$pheno
    genoMatNewAll <- resPast$genoMat
    indAll <- rownames(genoMatNewAll)

    selectedInd <- unique(c(resPast$selectedCross, resPast$selectedSelf))
    noSelectedInd <- indAll[!(indAll %in% selectedInd)]

    set.seed(phenoNewAll[1, 1])
    trialSelectedInd <- sort(sample(
      x = noSelectedInd,
      size = nPhenoField,
      replace = FALSE
    ))
    phenoNew <- phenoNewAll[trialSelectedInd, ]
    genoMatNew <- genoMatNewAll[trialSelectedInd, ]

    genoMat <- rbind(genoMat, genoMatNew)
    pheno <- rbind(pheno, phenoNew)
    newGP <- GP(phenoMat = pheno, genoMat = genoMat - 1)
    markerEffectEstimatedMat <- newGP$markerEffPredMat
    beta <- newGP$betaEstimated

    D2 <- CalcD2(
      ObjectMap = map,
      markerEffectMat = markerEffectEstimatedMat,
      nTrait = nTrait,
      k = f - 1
    )
    D2Self <- CalcD2(
      ObjectMap = map,
      markerEffectMat = markerEffectEstimatedMat,
      nTrait = nTrait,
      k = f - 2
    )
  }

  trueU <- (genoMatNow - 1) %*% markerEffectTrueMat
  trueU <- matrix(trueU, nrow = nrow(genoMatNow), ncol = nTrait)
  rownames(trueU) <- rownames(genoMatNow)

  # Population statistics for this generation
  AlleleFreq <- apply(genoMatNow, 2, function(eachAllele) {
    p <- sum(eachAllele) / (2 * length(eachAllele))
    return(4 * p * (1 - p))
  })
  genicVar <- apply(markerEffectTrueMat, 2, function(x) AlleleFreq %*% x^2)
  geneticVar <- apply(trueU, 2, function(x) {
    var(x) * (length(x) - 1) / length(x)
  })
  geneticCov <- sapply(1:ncol(covPairsIdx), function(cIdx) {
    t1 <- covPairsIdx[1, cIdx]
    t2 <- covPairsIdx[2, cIdx]
    var(trueU[, t1], trueU[, t2]) * (nrow(trueU) - 1) / nrow(trueU)
  })
  bulmer <- geneticVar / genicVar

  betaMat <- matrix(beta, nrow = nrow(genoMatNow), ncol = nTrait, byrow = TRUE)
  predU <- (genoMatNow - 1) %*% markerEffectEstimatedMat + betaMat
  predU <- matrix(predU, nrow = nrow(genoMatNow), ncol = nTrait)
  rownames(predU) <- rownames(genoMatNow)

  phenoNow <- GetPheno(uMat = trueU, sigmaE = sigmaE, nPheno = 1)

  predProgMean <- CalcProgMean(predU = predU)
  indName <- rownames(predU)
  comb <- combn(indName, 2)
  combInd <- combn(1:length(indName), 2)
  combName <- apply(comb, 2, function(x) paste0(x[1], "_", x[2]))

  gametArray <- ExtractGamet(
    nowPop = nowPop,
    indName = indName,
    qtn = rownames(map)
  )
  dimnames(gametArray)[[2]] <- colnames(genoMat)

  sigmaMat <- CalcVarMat(
    D2 = D2,
    combInd = combInd,
    gametArray = gametArray,
    nTrait = nTrait
  )
  rownames(sigmaMat) <- rownames(predProgMean) <- combName

  sigmaTrueMat <- CalcVarMat(
    D2 = D2True,
    combInd = combInd,
    gametArray = gametArray,
    nTrait = nTrait
  )
  rownames(sigmaTrueMat) <- combName

  w <- wInit + (wFinal - wInit) * generation / nGeneration
  meanPredU <- apply(predU, 2, mean)
  sdPredU <- apply(predU, 2, sd)

  mcResult <- EvalByMC(
    meanMat = predProgMean,
    varMat = sigmaMat,
    w = w,
    lVec = lVec,
    hVec = hVec,
    prob = probCross,
    nMC = nMCCross,
    meanPredU = meanPredU,
    sdPredU = sdPredU,
    nTrait = nTrait
  )
  predProgMat <- matrix(mcResult, ncol = 1)
  rownames(predProgMat) <- rownames(predProgMean)
  colnames(predProgMat) <- "MC"

  crossTable <- LP(
    predProgMat = predProgMat,
    indName = indName,
    comb = comb,
    capa = capa,
    nCrosses = nCrosses,
    nProg = nProg,
    generation = generation
  )

  nextPop <- population$new(
    name = paste0("C", generation, " offspring"),
    inds = makeCrosses(crosses = crossTable, pop = nowPop)
  )

  # Selfing evaluation (Inbred8) every 2 generations
  if (generation %% 2 == 0) {
    selfInd <- matrix(rep(1:nrow(predU), 2), nrow = 2, byrow = TRUE)
    varSelfMat <- CalcVarMat(
      D2 = D2Self,
      combInd = selfInd,
      gametArray = gametArray,
      nTrait = nTrait
    )
    mcSelfResult <- EvalByMC(
      meanMat = predU,
      varMat = varSelfMat,
      w = w,
      lVec = lVec,
      hVec = hVec,
      prob = probSelf,
      nMC = nMCSelf,
      meanPredU = meanPredU,
      sdPredU = sdPredU,
      nTrait = nTrait,
      hardConstraint = TRUE
    )
    selfSelected <- rownames(predU)[order(mcSelfResult, decreasing = TRUE)[
      1:nSelf
    ]]

    f8Res <- CreateF8(
      selfSelected = selfSelected,
      nowPop = nowPop,
      nProgSelf = nProgSelf,
      markerEffectEstimated = markerEffectEstimatedMat,
      markerEffectTrue = markerEffectTrueMat,
      marker = colnames(genoMatNow),
      beta = beta,
      sigmaE = sigmaE,
      nPheno = 1,
      nTrait = nTrait,
      coeff = coeff,
      generation = generation
    )
    f8History[[as.character(generation)]] <- f8Res
    selectedSelf <- selfSelected
  } else {
    selectedSelf <- NULL
  }
  selectedCross <- unique(c(crossTable$ind1, crossTable$ind2))

  resHistory[[generation + 1]] <- list(
    pheno = phenoNow,
    genoMat = genoMatNow,
    predU = predU,
    trueU = trueU,
    markerEffectEstimatedMat = markerEffectEstimatedMat,
    beta = beta,
    sigmaMat = sigmaMat,
    sigmaTrueMat = sigmaTrueMat,
    alleleFreq = AlleleFreq,
    genicVar = genicVar,
    geneticVar = geneticVar,
    geneticCov = geneticCov,
    bulmer = bulmer,
    selectedCross = selectedCross,
    selectedSelf = selectedSelf,
    w = w,
    generation = generation
  )
}

###### 5. Per-generation summary (extended to J = 3) ######

# Helper: pairwise scatter of true u with bound lines on essential traits
FigTrueUMulti <- function(plotPath, trueU, lVec, hVec, generation) {
  pairsIdx <- combn(nTrait, 2)
  png(plotPath, width = 480 * ncol(pairsIdx), height = 480)
  par(mfrow = c(1, ncol(pairsIdx)))
  for (cIdx in 1:ncol(pairsIdx)) {
    t1 <- pairsIdx[1, cIdx]
    t2 <- pairsIdx[2, cIdx]
    plot(
      trueU[, t1],
      trueU[, t2],
      xlab = paste0("Trait", t1),
      ylab = paste0("Trait", t2),
      main = paste0("Generation ", generation, ": Trait", t1, " vs Trait", t2)
    )
    if (t1 == 1 && t2 > 1) {
      abline(h = c(lVec[t2 - 1], hVec[t2 - 1]), lty = 2, col = "red")
    }
    if (t2 == 1 && t1 > 1) {
      abline(v = c(lVec[t1 - 1], hVec[t1 - 1]), lty = 2, col = "red")
    }
    if (t1 > 1 && t2 > 1) {
      abline(v = c(lVec[t1 - 1], hVec[t1 - 1]), lty = 2, col = "red")
      abline(h = c(lVec[t2 - 1], hVec[t2 - 1]), lty = 2, col = "red")
    }
  }
  dev.off()
}

resList <- lapply(resHistory, function(r) {
  generation <- r$generation
  trueU <- r$trueU
  predU <- r$predU
  markerEffectEstimatedMat <- r$markerEffectEstimatedMat
  sigmaMat <- r$sigmaMat
  sigmaTrueMat <- r$sigmaTrueMat

  # Figures
  FigTrueUMulti(
    paste0(dirSave, i, "_u_", generation, ".png"),
    trueU,
    lVec,
    hVec,
    generation
  )

  top10 <- sort(trueU[, 1], decreasing = TRUE)[1:10]
  topVal <- max(trueU[, 1])
  top10Mean <- mean(top10)

  inRange <- rep(TRUE, nrow(trueU))
  for (eIdx in 1:(nTrait - 1)) {
    inRange <- inRange &
      (trueU[, eIdx + 1] >= lVec[eIdx] & trueU[, eIdx + 1] <= hVec[eIdx])
  }
  effectiveProp <- mean(inRange)
  if (any(inRange)) {
    geneticValEff <- trueU[inRange, 1]
    effectiveTop <- max(geneticValEff)
    effTop10 <- sort(geneticValEff, decreasing = TRUE)[
      1:min(10, length(geneticValEff))
    ]
    effectiveTop10Mean <- mean(effTop10)
  } else {
    effectiveTop <- NA
    effectiveTop10Mean <- NA
  }

  accuracy <- sapply(1:nTrait, function(tr) cor(trueU[, tr], predU[, tr]))

  genoMat <- r$genoMat
  nMarker <- ncol(genoMat)
  alleleState <- apply(genoMat, 2, sum)
  alleleFixed <- alleleState == 0 | alleleState == nrow(genoMat) * 2
  nFixedAllele <- sum(alleleFixed)

  corMarker <- sapply(1:nTrait, function(tr) {
    cor(
      markerEffectTrueMat[!alleleFixed, tr],
      markerEffectEstimatedMat[!alleleFixed, tr]
    )
  })

  accuracyVarCov <- sapply(1:nVarCovCol, function(cc) {
    cor(sigmaTrueMat[, cc], sigmaMat[, cc])
  })
  rmseVarCov <- sapply(1:nVarCovCol, function(cc) {
    sqrt(mean((sigmaTrueMat[, cc] - sigmaMat[, cc])^2))
  })

  res <- c(
    topVal,
    top10Mean,
    effectiveProp,
    effectiveTop,
    effectiveTop10Mean,
    r$geneticVar,
    r$geneticCov,
    r$genicVar,
    r$bulmer,
    accuracy,
    accuracyVarCov,
    rmseVarCov,
    nFixedAllele,
    corMarker
  )
  names(res) <- c(
    "top",
    "top10Mean",
    "effectiveProp",
    "effectiveTop",
    "effectiveTop10Mean",
    paste0("geneticVar", 1:nTrait),
    paste0("geneticCov", covLabels),
    paste0("genicVar", 1:nTrait),
    paste0("bulmer", 1:nTrait),
    paste0("accuracy", 1:nTrait),
    paste0("accuracyVar", 1:nTrait),
    paste0("accuracyCov", covLabels),
    paste0("rmseVar", 1:nTrait),
    paste0("rmseCov", covLabels),
    "nFixedAllele",
    paste0("corMarker", 1:nTrait)
  )
  return(res)
})
resMat <- do.call(rbind, resList)
rownames(resMat) <- paste0(
  "C",
  formatC(0:nGeneration, width = 2, flag = "0")
)
write.csv(resMat, file = paste0(dirSave, i, "_geneticValMat.csv"))

# F8 summary (Inbred8 evaluation at even generations)
f8ResList <- lapply(f8History, function(r) {
  generation <- r$generation
  trueU <- r$trueU
  predU <- r$predU

  FigTrueUMulti(
    paste0(dirSave, i, "_F", f, "_", generation, ".png"),
    trueU,
    lVec,
    hVec,
    generation
  )

  top10 <- sort(trueU[, 1], decreasing = TRUE)[1:10]
  topVal <- max(trueU[, 1])
  top10Mean <- mean(top10)

  inRange <- rep(TRUE, nrow(trueU))
  for (eIdx in 1:(nTrait - 1)) {
    inRange <- inRange &
      (trueU[, eIdx + 1] >= lVec[eIdx] & trueU[, eIdx + 1] <= hVec[eIdx])
  }
  effectiveProp <- mean(inRange)
  if (any(inRange)) {
    geneticValEff <- trueU[inRange, 1]
    effectiveTop <- max(geneticValEff)
    effTop10 <- sort(geneticValEff, decreasing = TRUE)[
      1:min(10, length(geneticValEff))
    ]
    effectiveTop10Mean <- mean(effTop10)
  } else {
    effectiveTop <- NA
    effectiveTop10Mean <- NA
  }

  accuracy <- sapply(1:nTrait, function(tr) cor(trueU[, tr], predU[, tr]))

  res <- c(
    topVal,
    top10Mean,
    effectiveProp,
    effectiveTop,
    effectiveTop10Mean,
    r$geneticVar,
    r$geneticCov,
    r$genicVar,
    r$bulmer,
    accuracy
  )
  names(res) <- c(
    "top",
    "top10Mean",
    "effectiveProp",
    "effectiveTop",
    "effectiveTop10Mean",
    paste0("geneticVar", 1:nTrait),
    paste0("geneticCov", covLabels),
    paste0("genicVar", 1:nTrait),
    paste0("bulmer", 1:nTrait),
    paste0("accuracy", 1:nTrait)
  )
  return(res)
})
f8ResMat <- do.call(rbind, f8ResList)
rownames(f8ResMat) <- paste0(
  "C",
  formatC(as.numeric(names(f8History)), width = 2, flag = "0")
)
write.csv(f8ResMat, file = paste0(dirSave, i, "_resultMat.csv"))

# Time-series plot of effectiveTop (F8 / Inbred8 population)
f8Gen <- as.numeric(names(f8History))
effectiveTopSeries <- f8ResMat[, "effectiveTop"]
effectiveTopPath <- paste0(dirSave, i, "_F", f, "_effectiveTop.png")
png(effectiveTopPath, width = 720, height = 480)
plot(
  f8Gen,
  effectiveTopSeries,
  type = "b",
  pch = 19,
  col = "steelblue",
  xlab = "Generation",
  ylab = "effectiveTop (max trait1 within essential bounds)",
  main = paste0("F", f, " effectiveTop over generations")
)
grid()
dev.off()
