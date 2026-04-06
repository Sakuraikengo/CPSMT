##############################################
# Title : 4.0.PCS
# Author : Kengo Sakurai
# Date : 2023-07-21
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(RAINBOWR)
library(gaston)
library(lme4)
library(breedSimulatR)
library(lpSolve)
library(stringr)
library(doParallel)
library(tictoc)
library(abind)
library(MASS)
packages <- c(
  "RAINBOWR",
  "gaston",
  "lme4",
  "breedSimulatR",
  "lpSolve",
  "stringr",
  "abind",
  "MASS"
)
source(file = "scripts/1.0.function.R")

# 1.1. Parameters
nProg <- 15 # the number of progeny from each cross
nProgSelf <- 50 # the number of progeny for selfing
probSelf <- 0.9 # probability (which will decide the selection intensity for selfing)
capa <- 2 # how many times can we use each genotype for cross
nCrosses <- 10 # the number of crosses
nSelf <- 5 # the number of selfing genotypes for releasing varieties
nGeneration <- 20 # the number of generations
nTrait <- 2
nRep <- 30
nPhenoField <- 100
nMCSelf <- 10000 # the number of Monte Carlo simulation for selfing evaluation

# seven scenarios
gCorListList <- list(
  list(
    typeCor = "NoRelation",
    type = "NoRelation",
    value = matrix(c(1, 0, 0, 1), nrow = 2),
    coeff = NA
  ),
  list(
    typeCor = "Pleiotropy",
    type = "Positive",
    value = matrix(c(1, 0.6, 0.6, 1), nrow = 2),
    coeff = NA
  ),
  list(
    typeCor = "Pleiotropy",
    type = "NonLinear1_1",
    value = matrix(c(1, 0, 0, 1), nrow = 2),
    coeff = c(-1, 0)
  ),
  list(
    typeCor = "Pleiotropy",
    type = "NonLinear1_2",
    value = matrix(c(1, 0, 0, 1), nrow = 2),
    coeff = c(-1, 2)
  ),
  list(
    typeCor = "SpuriousPleiotropy",
    type = "Positive",
    value = matrix(c(1, 0.6, 0.6, 1), nrow = 2),
    coeff = NA
  ),
  list(
    typeCor = "SpuriousPleiotropy",
    type = "NonLinear1_1",
    value = matrix(c(1, 0, 0, 1), nrow = 2),
    coeff = c(-1, 0)
  ),
  list(
    typeCor = "SpuriousPleiotropy",
    type = "NonLinear1_2",
    value = matrix(c(1, 0, 0, 1), nrow = 2),
    coeff = c(-1, 2)
  )
)

wCombMat <- expand.grid(
  wInit = c(0.9, 0.5, 0.1),
  wFinal = c(0.9, 0.5, 0.1)
)
h <- 1.08 # upper line
l <- -1.08 # lower line
f <- 8 # the generation for evaluation (F8)

for (wIdx in 1:nrow(wCombMat)) {
  wInit <- wCombMat[wIdx, "wInit"]
  wFinal <- wCombMat[wIdx, "wFinal"]
  # wInit <- wInitList[1]; wFinal <- wFinalList[1]
  print(paste0("wInit=", wInit, ", wFinal=", wFinal))

  for (gCorList in gCorListList) {
    # gCorList <- gCorListList[[1]]
    # 1.2. Save dir
    typeCor <- gCorList$typeCor
    gCorName <- gCorList$type
    coeff <- gCorList$coeff
    dirSave <- paste0(
      "midstream/4.0.PCS/",
      wInit,
      "_",
      wFinal,
      "/",
      typeCor,
      "/",
      gCorName,
      "/"
    )
    if (!dir.exists(dirSave)) {
      dir.create(dirSave, recursive = T)
    }

    cl <- makeCluster(30)
    registerDoParallel(cl)
    foreach(
      i = 1:nRep,
      .export = c(
        "dirSave",
        "typeCor",
        "gCorName",
        "coeff",
        "nTrait",
        "nProg",
        "nProgSelf",
        "nCrosses",
        "nSelf",
        "nGeneration",
        "nPhenoField",
        "nMCSelf",
        "probSelf",
        "capa",
        "f",
        "l",
        "h",
        "wInit",
        "wFinal"
      ),
      .packages = packages,
      .combine = ignore
    ) %dopar%
      {
        # i <- 1
        resPath <- paste0(dirSave, i, "_resultMat.csv")
        if (!file.exists(resPath)) {
          # 1.3. Read data
          generation <- 0
          simPath <- paste0(
            "midstream/2.0.simulationSetting/",
            typeCor,
            "/"
          )
          initPop <- readRDS(file = paste0(simPath, i, "_initPop.rds"))
          genoMat <- initPop$genoMat
          map <- read.csv(file = paste0(simPath, i, "_map.csv"), row.names = 1)
          initRes <- readRDS(
            file = paste0(simPath, gCorName, "/", i, "_initResList.rds")
          )
          trueU <- initRes$uMat
          pheno <- initRes$pheno
          # blup <- CalcBlup(phenoMat = pheno)
          sigmaE <- initRes$sigmaE

          markerEffectTrueMat <- read.csv(
            file = paste0(simPath, gCorName, "/", i, "_markerEffects.csv"),
            row.names = 1
          )
          markerEffectTrueMat <- as.matrix(markerEffectTrueMat)
          # building a GP model
          model <- GP(phenoMat = pheno, genoMat = genoMat - 1)
          predU <- model$uPredMat
          markerEffectEstimatedMat <- model$markerEffPredMat
          beta <- model$betaEstimated

          # Computing bulmer effect
          AlleleFreq <- apply(genoMat, 2, function(eachAllele) {
            p <- sum(eachAllele) / (2 * length(eachAllele))
            return(4 * p * (1 - p))
          })
          genicVar <- apply(markerEffectTrueMat, 2, function(x) {
            return(AlleleFreq %*% x^2)
          })
          geneticVar <- apply(trueU, 2, function(x) {
            var(x) * (length(x) - 1) / length(x)
          })
          geneticCov <- var(trueU[, 1], trueU[, 2]) *
            (nrow(trueU) - 1) /
            nrow(trueU)
          bulmer <- geneticVar / genicVar

          if (i <= 5) {
            # draw a prediction accuracy figure of genotypic values
            gpPath <- paste0(dirSave, i, "_GP_", generation, ".png")
            FigPredU(gpPath, trueU, predU)
          }

          # Calculating D2 matrix
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

          ###### 2. Analysis ######
          #### Do the simulation ####
          # Calculate w for this generation
          w <- wInit + (wFinal - wInit) * generation / nGeneration

          # Standardize predU
          meanPredU <- apply(predU, 2, mean)
          sdPredU <- apply(predU, 2, sd)

          # producing all combinations
          indName <- rownames(predU)
          comb <- combn(indName, 2)
          combName <- apply(comb, 2, function(x) {
            paste0(x[1], "_", x[2])
          })

          # calculating the mean of genotypic value for each cross
          predProgMean <- CalcProgMean(predU = predU)
          rownames(predProgMean) <- combName

          # add the penalty term based on the distance from the desired range for each cross
          penalty <- ifelse(
            predProgMean[, 2] >= l & predProgMean[, 2] <= h,
            0,
            pmax(
              (predProgMean[, 2] - h) / sdPredU[2],
              (l - predProgMean[, 2]) / sdPredU[2]
            )
          )

          # calculating the evaluated values
          eval <- w *
            (predProgMean[, 1] - meanPredU[1]) /
            sdPredU[1] -
            (1 - w) * penalty
          predProgMat <- matrix(eval, ncol = 1)
          rownames(predProgMat) <- combName

          # deciding crossing pairs
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

          resEachList <- list(
            pheno = pheno,
            genoMat = genoMat,
            predU = predU,
            trueU = trueU,
            markerEffectEstimatedMat = markerEffectEstimatedMat,
            beta = beta,
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

          saveRDS(
            resEachList,
            paste0(dirSave, i, "_resList_", generation, ".rds")
          )
          rm(resEachList)
          gc()

          resF8List <- list(
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

          saveRDS(
            resF8List,
            paste0(dirSave, i, "_F", f, "_", generation, ".rds")
          )
          rm(resF8List)
          gc()

          # creating the next generation
          nextPop <- population$new(
            name = "C0 offspring",
            inds = makeCrosses(crosses = crossTable, pop = initPop)
          )

          ###### creating the next pop ########
          # do the GP
          for (generation in 1:nGeneration) {
            # generation <- 1
            nowPop <- nextPop
            genoMatNow <- nowPop$genoMat

            # update the GP model using remained individuals
            if ((generation > 3) & (generation %% 2 == 0)) {
              # add the field trial data
              pathPast <- paste0(
                dirSave,
                i,
                "_resList_",
                (generation - 1),
                ".rds"
              )
              resPast <- readRDS(pathPast)
              phenoNewAll <- resPast$pheno
              genoMatNewAll <- resPast$genoMat
              indAll <- rownames(genoMatNewAll)

              # removing the individuals used in cross or self
              selectedInd <- unique(c(
                resPast$selectedCross,
                resPast$selectedSelf
              ))
              noSelectedInd <- indAll[!(indAll %in% selectedInd)]

              # selecting the individuals for field trial
              set.seed(phenoNewAll[1, 1])
              trialSelectedInd <- sample(
                x = noSelectedInd,
                size = nPhenoField,
                replace = F
              )
              trialSelectedInd <- sort(trialSelectedInd)

              phenoNew <- phenoNewAll[trialSelectedInd, ]
              genoMatNew <- genoMatNewAll[trialSelectedInd, ]

              genoMat <- rbind(genoMat, genoMatNew)
              pheno <- rbind(pheno, phenoNew)
              newGP <- GP(phenoMat = pheno, genoMat = genoMat - 1)
              markerEffectEstimatedMat <- newGP$markerEffPredMat
              beta <- newGP$betaEstimated

              # Calculating D2 matrix
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

            trueU <- c((genoMatNow - 1) %*% markerEffectTrueMat)
            trueU <- matrix(trueU, nrow = nrow(genoMatNow), ncol = nTrait)
            rownames(trueU) <- rownames(genoMatNow)

            if (all(!is.na(coeff))) {
              # converting the genetic correlation to non-linear
              u1NonLinear <- coeff[1] * ((trueU[, 2] - coeff[2])^2) / 4
              trueU[, 1] <- trueU[, 1] + u1NonLinear
            }

            # Computing bulmer effect
            AlleleFreq <- apply(genoMatNow, 2, function(eachAllele) {
              p <- sum(eachAllele) / (2 * length(eachAllele))
              return(4 * p * (1 - p))
            })
            genicVar <- apply(markerEffectTrueMat, 2, function(x) {
              return(AlleleFreq %*% x^2)
            })
            geneticVar <- apply(trueU, 2, function(x) {
              var(x) * (length(x) - 1) / length(x)
            })
            geneticCov <- var(trueU[, 1], trueU[, 2]) *
              (nrow(trueU) - 1) /
              nrow(trueU)
            bulmer <- geneticVar / genicVar

            betaMat <- matrix(
              beta,
              nrow = nrow(genoMatNow),
              ncol = nTrait,
              byrow = T
            )
            predU <- (genoMatNow - 1) %*% markerEffectEstimatedMat + betaMat
            predU <- matrix(predU, nrow = nrow(genoMatNow), ncol = nTrait)
            rownames(predU) <- rownames(genoMatNow)

            if (i <= 5) {
              # draw a prediction accuracy figure of genotypic values
              gpPath <- paste0(dirSave, i, "_GP_", generation, ".png")
              FigPredU(gpPath, trueU, predU)
            }

            # collect phenotypic data
            phenoNow <- GetPheno(uMat = trueU, sigmaE = sigmaE, nPheno = 1)

            # Calculate w for this generation
            w <- wInit + (wFinal - wInit) * generation / nGeneration

            # Standardize predU
            meanPredU <- apply(predU, 2, mean)
            sdPredU <- apply(predU, 2, sd)

            # producing all combinations
            indName <- rownames(predU)
            comb <- combn(indName, 2)
            combInd <- combn(1:length(indName), 2)
            combName <- apply(comb, 2, function(x) {
              paste0(x[1], "_", x[2])
            })

            # calculating the mean of genotypic value for each cross
            predProgMean <- CalcProgMean(predU = predU)
            rownames(predProgMean) <- combName

            # add the penalty term based on the distance from the desired range for each cross
            penalty <- ifelse(
              predProgMean[, 2] >= l & predProgMean[, 2] <= h,
              0,
              pmax(
                (predProgMean[, 2] - h) / sdPredU[2],
                (l - predProgMean[, 2]) / sdPredU[2]
              )
            )

            # calculating the evaluated values
            eval <- w *
              (predProgMean[, 1] - meanPredU[1]) /
              sdPredU[1] -
              (1 - w) * penalty
            predProgMat <- matrix(eval, ncol = 1)
            rownames(predProgMat) <- combName

            # deciding crossing pairs
            crossTable <- LP(
              predProgMat = predProgMat,
              indName = indName,
              comb = comb,
              capa = capa,
              nCrosses = nCrosses,
              nProg = nProg,
              generation = generation
            )

            # extracting the gamet information
            gametArray <- ExtractGamet(
              nowPop = nowPop,
              indName = indName,
              qtn = rownames(map)
            )
            dimnames(gametArray)[[2]] <- colnames(genoMatNow)

            # creating the next generation
            nextPop <- population$new(
              name = paste0("C", generation, " offspring"),
              inds = makeCrosses(crosses = crossTable, pop = nowPop)
            )

            ###################################################
            # selecting the genotypes for releasing varieties
            if (generation %% 2 == 0) {
              selfInd <- matrix(rep(1:nrow(predU), 2), nrow = 2, byrow = T)
              varSelfMat <- CalcVarMat(
                D2 = D2Self,
                combInd = selfInd,
                gametArray = gametArray,
                nTrait = nTrait
              )
              # Monte Carlo simulation for selfing
              mcSelfResult <- EvalByMC(
                meanMat = predU,
                varMat = varSelfMat,
                w = w,
                l = l,
                h = h,
                prob = probSelf,
                nMC = nMCSelf,
                meanPredU = meanPredU,
                sdPredU = sdPredU,
                hardConstraint = TRUE
              )
              selectedInd <- order(mcSelfResult, decreasing = T)[1:nSelf]
              selfSelected <- rownames(predU)[selectedInd]

              # list contains "phenotypic data", "genome data",
              # "estimated genetic value", and "true genetic value"
              resF8List <- CreateF8(
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
              selectedSelf <- selfSelected
              # resEachList <- c(resEachList, list(resEachListNew))
              saveRDS(
                resF8List,
                file = paste0(dirSave, i, "_F", f, "_", generation, ".rds")
              )
            } else {
              selectedSelf <- NULL
            }
            selectedCross <- unique(c(crossTable$ind1, crossTable$ind2))
            resEachList <- list(
              pheno = phenoNow,
              genoMat = genoMatNow,
              predU = predU,
              trueU = trueU,
              markerEffectEstimatedMat = markerEffectEstimatedMat,
              beta = beta,
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
            saveRDS(
              resEachList,
              paste0(dirSave, i, "_resList_", generation, ".rds")
            )
            rm(resEachList)
            gc()
          }
          # summary of the recurrent genomic selection
          resFilePathList <- paste0(
            dirSave,
            i,
            "_resList_",
            0:nGeneration,
            ".rds"
          )

          resList <- lapply(resFilePathList, function(resFilePathEach) {
            # resFilePathEach <- resFilePathList[[1]]
            resGenereationEach <- readRDS(resFilePathEach)
            generation <- resGenereationEach$generation
            trueU <- resGenereationEach$trueU
            predU <- resGenereationEach$predU
            markerEffectEstimatedMat <- resGenereationEach$markerEffectEstimatedMat
            if (i == 1) {
              # draw a figure of true genotypic values
              plotPath0 <- gsub(
                resFilePathEach,
                pattern = "resList",
                replacement = "u"
              )
              plotPath <- gsub(
                plotPath0,
                pattern = ".rds",
                replacement = ".png"
              )
              FigTrueU(plotPath, trueU, h, l, generation)
            }

            top10 <- sort(trueU[, 1], decreasing = T)[1:10]
            top <- sort(trueU[, 1], decreasing = T)[1]
            top10Mean <- mean(top10)

            effectiveInd <- (l <= trueU[, 2]) & (trueU[, 2] <= h)
            effectiveProp <- sum(effectiveInd) / length(effectiveInd)
            geneticValEffective <- trueU[effectiveInd, 1]
            effectiveTop10 <- sort(geneticValEffective, decreasing = T)[1:10]
            effectiveTop <- max(geneticValEffective)
            effectiveTop10Mean <- mean(effectiveTop10)

            gVariance <- resGenereationEach$geneticVar
            gCov <- resGenereationEach$geneticCov

            genicVar <- resGenereationEach$genicVar
            bulmer <- resGenereationEach$bulmer

            accuracy1 <- cor(trueU[, 1], predU[, 1])
            accuracy2 <- cor(trueU[, 2], predU[, 2])

            # count the number of fixed alleles
            genoMat <- resGenereationEach$genoMat
            nMarker <- ncol(genoMat)
            alleleState <- apply(genoMat, 2, sum)
            alleleFixed <- alleleState == 0 | alleleState == nrow(genoMat) * 2
            nFixedAllele <- sum(alleleFixed)

            # compute the prediction accuracy of marker effects
            corMarker1 <- cor(
              markerEffectTrueMat[!alleleFixed, 1],
              markerEffectEstimatedMat[!alleleFixed, 1]
            )
            corMarker2 <- cor(
              markerEffectTrueMat[!alleleFixed, 2],
              markerEffectEstimatedMat[!alleleFixed, 2]
            )

            if (i == 1) {
              markerPredictPath <- paste0(
                dirSave,
                i,
                "_MarkerPrediction_",
                generation,
                ".png"
              )

              # draw a figure of marker effects prediction accuracy
              FigPredMarker(
                markerPredictPath,
                markerEffectEstimatedMat,
                markerEffectTrueMat,
                alleleFixed,
                corMarker1,
                corMarker2,
                nMarker,
                nFixedAllele
              )
            }
            res <- c(
              top,
              top10Mean,
              effectiveProp,
              effectiveTop,
              effectiveTop10Mean,
              gVariance,
              gCov,
              genicVar,
              bulmer,
              accuracy1,
              accuracy2,
              nFixedAllele,
              corMarker1,
              corMarker2
            )
            names(res) <- c(
              "top",
              "top10Mean",
              "effectiveProp",
              "effectiveTop",
              "effectiveTop10Mean",
              "geneticVar1",
              "geneticVar2",
              "geneticCov",
              "genicVar1",
              "genicVar2",
              "bulmer1",
              "bulmer2",
              "accuracy1",
              "accuracy2",
              "nFixedAllele",
              "corMarker1",
              "corMarker2"
            )
            # plot(geneticValueEstimated, trueU)
            return(res)
          })
          resMat <- do.call(rbind, resList)
          rownames(resMat) <- paste0(
            "C",
            formatC(0:nGeneration, width = 2, flag = "0")
          )
          write.csv(resMat, file = paste0(dirSave, i, "_geneticValMat.csv"))

          # summary of the Inbred8 population
          f8FilePathList <- paste0(
            dirSave,
            i,
            "_F8_",
            seq(0, nGeneration, 2),
            ".rds"
          )

          f8ResList <- lapply(f8FilePathList, function(resFilePathEach) {
            # resFilePathEach <- f8FilePathList[[10]]
            resGenereationEach <- readRDS(resFilePathEach)
            generation <- resGenereationEach$generation
            trueU <- resGenereationEach$trueU
            predU <- resGenereationEach$predU

            if (i == 1) {
              # draw a figure of true genotypic values
              plotPath <- gsub(
                resFilePathEach,
                pattern = ".rds",
                replacement = ".png"
              )
              FigTrueU(plotPath, trueU, h, l, generation)
            }

            top10 <- sort(trueU[, 1], decreasing = T)[1:10]
            top <- sort(trueU[, 1], decreasing = T)[1]
            top10Mean <- mean(top10)

            effectiveInd <- (l <= trueU[, 2]) & (trueU[, 2] <= h)
            effectiveProp <- sum(effectiveInd) / length(effectiveInd)
            geneticValEffective <- trueU[effectiveInd, 1]
            effectiveTop10 <- sort(geneticValEffective, decreasing = T)[1:10]
            effectiveTop <- max(geneticValEffective)
            effectiveTop10Mean <- mean(effectiveTop10)

            gVariance <- resGenereationEach$geneticVar
            gCov <- resGenereationEach$geneticCov

            genicVar <- resGenereationEach$genicVar
            bulmer <- resGenereationEach$bulmer

            accuracy1 <- cor(trueU[, 1], predU[, 1])
            accuracy2 <- cor(trueU[, 2], predU[, 2])

            res <- c(
              top,
              top10Mean,
              effectiveProp,
              effectiveTop,
              effectiveTop10Mean,
              gVariance,
              gCov,
              genicVar,
              bulmer,
              accuracy1,
              accuracy2
            )
            names(res) <- c(
              "top",
              "top10Mean",
              "effectiveProp",
              "effectiveTop",
              "effectiveTop10Mean",
              "geneticVar1",
              "geneticVar2",
              "geneticCov",
              "genicVar1",
              "genicVar2",
              "bulmer1",
              "bulmer2",
              "accuracy1",
              "accuracy2"
            )
            # plot(geneticValueEstimated, trueU)
            return(res)
          })
          f8ResMat <- do.call(rbind, f8ResList)
          rownames(f8ResMat) <- paste0(
            "C",
            formatC(seq(0, nGeneration, 2), width = 2, flag = "0")
          )
          write.csv(f8ResMat, file = paste0(dirSave, i, "_resultMat.csv"))
        }
      }
    stopCluster(cl)
  }
}
