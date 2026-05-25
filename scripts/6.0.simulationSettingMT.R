##############################################
# Title : 6.0.simulationSettingMT
# Author : Kengo Sakurai
# Date : 2026-05-23
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(breedSimulatR)
library(stringr)
library(gaston)
library(rrBLUP)
library(glmnet)
library(MASS)
source(file = "scripts/1.0.function.R")

# 1.1. Parameters
thresholdMAF <- 0.05 # minimum MAF
thresholdLD <- 0.95 # maximum LD
nSNP <- 200 # number of SNPs in each chr
nChr <- 20 # number of chromosomes
lchr <- 1e9 # physical length of each chr
lchrCm <- 100 # map length of each chr
nQtn <- 10 # number of QTN in each chr (per trait)
nPheno <- 5 # number of replication in the init pop
nTrait <- 3 # number of trait
h2 <- 0.6 # heritability in the init pop
nLine <- 150 # number of accessions
nRep <- 1 # demo: single replicate
alpha <- 0 # 0 corresponds to ridge regression
lambda <- 0.1 # penalty value for ridge regression

# Only the NoRelation case (3 x 3 identity correlation)
gCorListList <- list(
  list(
    type = "NoRelation",
    value = diag(nTrait),
    coeff = NA
  )
)
typeCorList <- c("NoRelation")

# 1.2. Save dir
dirSave <- "midstream/6.0.simulationSettingMT/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave, recursive = TRUE)
}

seedIndCsv <- "midstream/seedInd.csv"
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(seedIndCsv, row.names = 1, header = TRUE)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:50000, 1000, replace = FALSE)
  write.csv(x = seedInd, file = seedIndCsv)
}

# 1.3. Build (or load cached) genome data
genomeDataListPath <- "midstream/2.0.simulationSetting/genomeDataList.rds"
genomeDataList <- readRDS(genomeDataListPath)
genomeMat <- genomeDataList$genomeMat
snpInfoDf <- genomeDataList$snpInfoDf

# Genetic relationship matrix (used to draw correlated init u)
K <- A.mat(X = genomeMat - 1)

###### 2. QTN, marker effects, and initial population per rep ######
for (typeCor in typeCorList) {
  dirSaveEach <- paste0(dirSave, typeCor, "/")
  if (!dir.exists(dirSaveEach)) {
    dir.create(dirSaveEach, recursive = TRUE)
  }

  for (i in 1:nRep) {
    specie_statEx <- specie$new(
      specName = "Soybean",
      nChr = nChr,
      lchr = lchr,
      lchrCm = lchrCm,
      verbose = FALSE
    )

    markersAll <- colnames(genomeMat)

    # NoRelation: three non-overlapping QTN sets
    set.seed(seedInd[i])
    qtn1 <- unlist(tapply(snpInfoDf$SNPid, snpInfoDf$chr, sample, nQtn))

    markersNoQtn1 <- markersAll[!(markersAll %in% qtn1)]
    markersNoQtn1Df <- snpInfoDf[markersNoQtn1, ]
    qtn2 <- unlist(tapply(
      markersNoQtn1Df$SNPid,
      markersNoQtn1Df$chr,
      sample,
      nQtn
    ))

    markersNoQtn12 <- markersAll[!(markersAll %in% c(qtn1, qtn2))]
    markersNoQtn12Df <- snpInfoDf[markersNoQtn12, ]
    qtn3 <- unlist(tapply(
      markersNoQtn12Df$SNPid,
      markersNoQtn12Df$chr,
      sample,
      nQtn
    ))

    qtnList <- list(qtn1, qtn2, qtn3)
    qtnMatList <- list(
      genomeMat[, qtn1],
      genomeMat[, qtn2],
      genomeMat[, qtn3]
    )
    qtnSel <- sort(unique(c(qtn1, qtn2, qtn3)))

    # Fill remaining markers up to nSNP per chromosome
    markersOthers <- markersAll[!(markersAll %in% qtnSel)]
    markersOthersDf <- snpInfoDf[markersOthers, ]
    nMarkers <- nSNP - length(qtnSel) / nChr
    snpSel <- unlist(tapply(
      markersOthersDf$SNPid,
      markersOthersDf$chr,
      sample,
      nMarkers,
      replace = FALSE
    ))
    snpSel <- sort(snpSel)
    markers <- sort(c(qtnSel, snpSel))
    snpCoord <- snpInfoDf[snpInfoDf$SNPid %in% markers, ]
    rownames(snpCoord) <- snpCoord$SNPid
    genomeMatSel <- genomeMat[, markers]

    if (ncol(genomeMatSel) != nSNP * nChr) {
      warning(paste0("Rep ", i, ": number of markers is ", ncol(genomeMatSel)))
    }

    SNPs <- SNPinfo$new(SNPcoord = snpCoord, specie = specie_statEx)
    map <- SNPs$SNPcoord
    map <- map[order(map$SNPid), ]
    write.csv(map, file = paste0(dirSaveEach, i, "_map.csv"))

    initPop <- createPop(
      geno = genomeMatSel,
      SNPinfo = SNPs,
      popName = "Initial population",
      verbose = FALSE
    )
    saveRDS(initPop, file = paste0(dirSaveEach, i, "_initPop.rds"))

    for (gCorList in gCorListList) {
      gCorName <- gCorList$type
      gCor <- gCorList$value
      coeff <- gCorList$coeff

      # NoRelation only: no nested subdirectory needed for this demo
      dirSaveCorEach <- dirSaveEach

      # True u from MVN with identity correlation (NoRelation)
      set.seed(seedInd[i])
      uInitVec <- mvrnorm(
        n = 1,
        mu = rep(0, nTrait * nrow(K)),
        Sigma = gCor %x% K
      )
      uInitMat <- matrix(uInitVec, ncol = nTrait, byrow = FALSE)

      # Ridge-regression marker effects per trait
      markerEffSelList <- vector("list", nTrait)
      for (tr in 1:nTrait) {
        fit <- glmnet(
          intercept = 0,
          x = qtnMatList[[tr]] - 1,
          y = uInitMat[, tr],
          alpha = alpha,
          lambda = lambda
        )
        coefTr <- coef(fit)
        markerEffSelList[[tr]] <- coefTr[2:nrow(coefTr), ]

        png(filename = paste0(dirSaveCorEach, i, "_markerEffect", tr, ".png"))
        hist(
          markerEffSelList[[tr]],
          main = paste0("Marker Effects for Trait", tr),
          xlab = "Marker Effect"
        )
        dev.off()
      }

      markerEffMat <- matrix(0, nrow = ncol(genomeMatSel), ncol = nTrait)
      rownames(markerEffMat) <- colnames(genomeMatSel)
      for (tr in 1:nTrait) {
        markerEffMat[qtnList[[tr]], tr] <- markerEffSelList[[tr]]
      }
      write.csv(
        markerEffMat,
        file = paste0(dirSaveCorEach, i, "_markerEffects.csv")
      )

      # True genotypic values
      uMat <- sapply(1:nTrait, function(tr) {
        (qtnMatList[[tr]] - 1) %*% markerEffSelList[[tr]]
      })
      rownames(uMat) <- rownames(genomeMatSel)

      # Pairwise scatter of true u (3 panels for J = 3)
      pairsIdx <- combn(nTrait, 2)
      png(
        filename = paste0(dirSaveCorEach, i, "_trueU.png"),
        width = 480 * ncol(pairsIdx),
        height = 480
      )
      par(mfrow = c(1, ncol(pairsIdx)))
      for (cIdx in 1:ncol(pairsIdx)) {
        t1 <- pairsIdx[1, cIdx]
        t2 <- pairsIdx[2, cIdx]
        plot(
          uMat[, t1],
          uMat[, t2],
          xlab = paste0("Trait", t1),
          ylab = paste0("Trait", t2),
          main = paste0(
            "g.Cor(t",
            t1,
            ", t",
            t2,
            ") = ",
            round(cor(uMat[, t1], uMat[, t2]), 2)
          )
        )
      }
      dev.off()

      # Field trial / residual variance
      sigmaG <- diag(var(uMat)) * (nrow(uMat) - 1) / nrow(uMat)
      sigmaE <- sigmaG / h2 - sigmaG
      phenoMat <- GetPheno(uMat = uMat, sigmaE = sigmaE, nPheno = nPheno)

      initResList <- list(
        pheno = phenoMat,
        uMat = uMat,
        genoMat = genomeMatSel,
        sigmaE = sigmaE
      )
      saveRDS(initResList, paste0(dirSaveCorEach, i, "_initResList.rds"))
    }
  }
}
