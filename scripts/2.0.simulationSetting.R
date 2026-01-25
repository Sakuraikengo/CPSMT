##############################################
# Title : 2.0.simulationSetting
# Author : Kengo Sakurai
# Date : 2025-01-18
##############################################

###### 1. Settings ######
# 1.0. Libraries
library(breedSimulatR)
library(stringr)
library(gaston)
library(lme4)
library(RAINBOWR)
library(rrBLUP)
library(glmnet)
library(MASS)
library(qtl)
# read functions
source(file = "scripts/1.0.function.R") 

# 1.1. Parameters
thresholdMAF <- 0.05 # for minimum MAF
thresholdLD <- 0.95 # for maximum LD
nSNP <- 200 # the number of SNPs in each chr
nChr <- 20 # the number of chromosomes
lchr <- 1e9 # physical length of each chr (no meaning, just for setting)
lchrCm <- 100 # map length of each chr
nQtn <- 10 # the number of QTN in each chr
nPheno <- 5 # the number of replication in the init pop
nTrait <- 2 # the number of trait
h2 <- 0.6 # heritability in the init pop
nLine <- 150 # the number of accessions
nRep <- 30 # the number of replications for breeding simulations

# four types of genetic correlations
gCorListList <- list(list(type = "NoRelation", 
                          value = matrix(c(1, 0, 0, 1), nrow = 2), 
                          coeff = NA), 
                     list(type = "Positive", 
                          value = matrix(c(1, 0.6, 0.6, 1), nrow = 2), 
                          coeff = NA), 
                     list(type = "NonLinear1_1", 
                          value = matrix(c(1, 0, 0, 1), nrow = 2), 
                          coeff = c(-1, 0)), 
                     list(type = "NonLinear1_2", 
                          value = matrix(c(1, 0, 0, 1), nrow = 2), 
                          coeff = c(-1, 2)))

# three types of genetic causality
typeCorList <- c("NoRelation", "Pleiotropy", "SpuriousPleiotropy")
linkTight <- 5 # map distance between 2 QTNs in the case of "SpuriousPleiotropy"
linkTightLD <- 0.8 # LD between 2 QTNs in the case of "SpuriousPleiotropy"
alpha <- 0  # 0 corresponds to ridge regression
lambda <- 0.1 # penalty value for ridge regression

# 1.2. Save dir
dirSave <- "midstream/2.0.simulationSetting/"
if (!dir.exists(dirSave)) {
  dir.create(dirSave)
}

# set the seed
seedIndCsv <- paste0("midstream/seedInd.csv")
if (file.exists(seedIndCsv)) {
  seedInd <- read.csv(paste0(seedIndCsv), row.names = 1, header = T)
  seedInd <- c(as.matrix(seedInd))
} else {
  seedInd <- sample(1:50000, 1000, replace = F)
  write.csv(x = seedInd, file = paste0(seedIndCsv))
}

# read genome data
genomeDataListPath <- paste0(dirSave, "genomeDataList.rds")
if (!file.exists(genomeDataListPath)) {
  genoInfoPath <- "raw_data/genotype/Gm198_HCDB_190207.fil.snp.remHet.MS0.95_bi_MQ20_DP3-1000.MAF0.025.imputed.v2.chrnum.vcf.gz"
  genoInfoRaw0 <- read.vcf(genoInfoPath)
  
  genoInfoRaw0@snps$id <- paste0("Chr", 
                                 formatC(genoInfoRaw0@snps$chr, width = 2, flag = "0"), 
                                 "_", 
                                 formatC(genoInfoRaw0@snps$pos, width = 8, flag = "0"))
  
  genomeMatRaw <- as.matrix(genoInfoRaw0)
  
  # Selecting the founder lines
  set.seed(seedInd[1])
  founderLine <- sample(x = rownames(genomeMatRaw), size = nLine, replace = F)
  write.csv(x = sort(founderLine), file = paste0(dirSave, "founderLine", nLine, ".csv"))
  genomeMatFounder <- genomeMatRaw[founderLine, ]
  # dim(genomeMatFounder)
  
  fam <- data.frame(famid = rownames(genomeMatFounder), 
                    id = rownames(genomeMatFounder), 
                    father = NA, 
                    mother = NA, 
                    sex = NA, 
                    pheno = NA)
  bim <- data.frame(chr = genoInfoRaw0@snps$chr, 
                    id = genoInfoRaw0@snps$id, 
                    dist = genoInfoRaw0@snps$dist, 
                    pos = genoInfoRaw0@snps$pos, 
                    A1 = genoInfoRaw0@snps$A1, 
                    A2 = genoInfoRaw0@snps$A2)
  
  bedMat <- as.bed.matrix(genomeMatFounder, fam, bim)
  genomeFil <- select.snps(bedMat, condition = maf >= thresholdMAF)
  genomeFil <- LD.thin(genomeFil, threshold = thresholdLD)
  genomeMat <- as.matrix(genomeFil)
  
  snpInfoList <- str_split(colnames(genomeMat), pattern = "_")
  snpInfoMat <- do.call(rbind, snpInfoList)
  snpInfoDf0 <- data.frame(chr = as.character(snpInfoMat[, 1]), 
                           SNPid = colnames(genomeMat), 
                           physPos = as.numeric(snpInfoMat[, 2]))
  snpInfoDf0 <- snpInfoDf0[order(snpInfoDf0$SNPid), ]
  
  # calculating map position
  mapPosList <- tapply(X = snpInfoDf0$physPos, INDEX = snpInfoDf0$chr, function(eachChr) {
    # eachChr <- snpInfoDf0$PysPos[snpInfoDf0$Chr == snpInfoDf0$Chr[1]]
    pysPosMax <- max(eachChr)
    pysPosMin <- min(eachChr)
    mapPos <- (eachChr -  pysPosMin) * 100 / (pysPosMax - pysPosMin)
    return(mapPos)
  }) 
  mapPos <- unlist(mapPosList)
  snpInfoDf <- data.frame(snpInfoDf0, 
                          linkMapPos = mapPos)
  rownames(snpInfoDf) <- snpInfoDf$SNPid
  genomeDataList <- list(genomeMat = genomeMat, 
                         snpInfoDf = snpInfoDf)
  saveRDS(genomeDataList, file = genomeDataListPath)
} else {
  genomeDataList <- readRDS(genomeDataListPath)
  genomeMat <- genomeDataList$genomeMat
  snpInfoDf <- genomeDataList$snpInfoDf
}

# computing the genetic relationship matrix
K <- A.mat(X = genomeMat - 1)

# simulate marker effects
for (typeCor in typeCorList) {
  # typeCor <- typeCorList[1]
  # 1.2. Save dir
  dirSaveEach <- paste0(dirSave, 
                        typeCor, "/")
  if (!dir.exists(dirSaveEach)) {
    dir.create(dirSaveEach, recursive = T)
  }
  
  for (i in 1:nRep) {
    # i <- 1
    # #### Creating the initial population ####
    # # create specie object
    specie_statEx <- specie$new(specName = "Soybean",
                                nChr = nChr,
                                lchr = lchr,
                                lchrCm = lchrCm, 
                                verbose = F)
    #### set the marker effect for trait1 ####
    # selecting the QTNs (totally nChr * nQtn)
    markersAll <- colnames(genomeMat)
    if (typeCor == "NoRelation") {
      qtn1 <- unlist(tapply(snpInfoDf$SNPid, snpInfoDf$chr, sample, nQtn))
      
      # Remove QTN markers
      markersNoQtn1 <- markersAll[!(markersAll %in% qtn1)]
      markersNoQtn1Df <- snpInfoDf[markersNoQtn1, ]
      qtn2 <- unlist(tapply(markersNoQtn1Df$SNPid, markersNoQtn1Df$chr, sample, nQtn))
      
      qtn <- cbind(qtn1, qtn2)
    } else {
      
      l <- 0
      # qtn2 <- NA
      qtnList <- lapply(unique(snpInfoDf$chr), function(chrEach) {
        # chrEach <- unique(snpInfoDf$chr)[1]
        qtn2Each <- rep(NA, nQtn)
        while (any(is.na(qtn2Each))) {
          l <- l + 1
          if (l == 10000) {
            stop()
          }
          
          set.seed(seedInd[i] + l)
          markerChrEach <- markersAll[snpInfoDf$chr == chrEach]
          qtn1Each <- sample(markerChrEach, nQtn, replace = F)
          qtn1Each <- sort(qtn1Each)
          
          #### set the marker effect for trait2 ####
          ##### Case1: Pleiotropy #####
          if (typeCor == "Pleiotropy") {
            qtn2Each <- qtn1Each
          } else if (typeCor == "SpuriousPleiotropy") {
            #### Case2: SpuriousPleiotropy ####
            for (qtnEachInd in 1:length(qtn1Each)) {
              # qtnEachInd <- 1
              iPlus <- 0
              qtnEach <- qtn1Each[qtnEachInd]
              qtlSelChr0 <- snpInfoDf[snpInfoDf$chr == str_split(qtnEach, "_")[[1]][1], ]
              qtlSelChr <- qtlSelChr0[!(qtlSelChr0$SNPid %in% qtn1Each), ]
              # calculating the distance between QTLs
              qtlDist <- abs(qtlSelChr$linkMapPos - qtlSelChr0[qtlSelChr0$SNPid == qtnEach, "linkMapPos"])
              
              qtlCand0 <- qtlSelChr[qtlDist <= linkTight, ]
              
              LD <- CalcLD(genoMat = genomeMat, 
                           qtnEach = qtnEach, 
                           qtlCand = qtlCand0$SNPid)
              qtlCand <- qtlCand0[LD > linkTightLD, ]
              if (nrow(qtlCand) < 1) {
                break
              }
              set.seed(seedInd[i])
              qtlSel <- sample(x = qtlCand$SNPid, size = 1)
              qtn2Each[qtnEachInd] <- qtlSel
              while (max(table(qtn2Each)) >= 2) {
                set.seed(seedInd[i] + iPlus)
                qtlSel <- sample(x = qtlCand$SNPid, size = 1)
                qtn2Each[qtnEachInd] <- qtlSel
                iPlus <- iPlus + 1
                if (iPlus > 100) {
                  qtn2Each[qtnEachInd] <- NA
                  print("Search again")
                  break
                }
              }
              # print(qtnEachInd)
            }
          }
        }
        return(cbind(qtn1Each, qtn2Each))
      })
      qtn <- do.call(rbind, qtnList)
    }
    
    if (all(!is.na(qtn))) {
      
      qtn1 <- qtn[, 1]
      qtn2 <- qtn[, 2]
      
      # extracting the qtn information
      qtnSel0 <- c(qtn1, qtn2)
      qtnSel <- sort(unique(qtnSel0))
      qtnMat1 <- genomeMat[, qtn1]
      qtnMat2 <- genomeMat[, qtn2]
      qtnMatSel <- genomeMat[, qtnSel]
      
      # Extracting the SNP markers (which don't have effect)
      markersOthers <- markersAll[!(markersAll %in% qtnSel)]
      markersOthersDf <- snpInfoDf[markersOthers, ]
      
      nMarkers <- nSNP - length(qtnSel) / nChr
      snpSel <- unlist(tapply(markersOthersDf$SNPid, markersOthersDf$chr, sample, nMarkers, replace = F))
      snpSel <- sort(snpSel)
      markers <- sort(c(qtnSel, snpSel))
      snpCoord <- snpInfoDf[snpInfoDf$SNPid %in% markers, ]
      rownames(snpCoord) <- snpCoord$SNPid
      genomeMatSel <- genomeMat[, markers]
      
      if (ncol(genomeMatSel) != 4000) {
        print("the number of markers is wrong")
      }
      
      # # create SNPinfo object
      SNPs <- SNPinfo$new(SNPcoord = snpCoord,
                          specie = specie_statEx)
      map <- SNPs$SNPcoord
      map <- map[order(map$SNPid), ]
      write.csv(map, file = paste0(dirSaveEach, i, "_map.csv"))
      
      # create population object
      initPop <- createPop(geno = genomeMatSel,
                           SNPinfo = SNPs,
                           popName = "Initial population", verbose = F)
      saveRDS(initPop, file = paste0(dirSaveEach, i, "_initPop.rds"))
      for (gCorList in gCorListList) {
        # gCorList <- gCorListList[[1]]
        
        gCorName <- gCorList$type
        gCor <- gCorList$value
        coeff <- gCorList$coeff
        
        dirSaveCorEach <- paste0(dirSaveEach, 
                                 gCorName, "/")
        if (!dir.exists(dirSaveCorEach)) {
          dir.create(dirSaveCorEach, recursive = T)
        }
        
        # make genetic correlation between 2 traits in the present population
        set.seed(seedInd[i])
        uInitVec <- mvrnorm(n = 1, mu = rep(0, nTrait*nrow(K)), Sigma = gCor %x% K)
        uInitMat <- matrix(uInitVec, ncol = nTrait, byrow = F)
        
        # estimate the marker effect using ridge regression
        fit1 <- glmnet(intercept = 0, x = qtnMat1 - 1, y = uInitMat[, 1], 
                       alpha = alpha, lambda = lambda)
        markerEffSel1 <- coef(fit1)[2:nrow(coef(fit1)), ]
        fit2 <- glmnet(intercept = 0, x = qtnMat2 - 1, y = uInitMat[, 2], 
                       alpha = alpha, lambda = lambda)
        markerEffSel2 <- coef(fit2)[2:nrow(coef(fit2)), ]
        
        # plot(uInitMat[, 1], (qtnMat1 - 1) %*% markerEffSel1, 
        #      main = round(cor(uInitMat[, 1], (qtnMat1 - 1) %*% markerEffSel1), 2))
        png(filename = paste0(dirSaveCorEach, i, "_markerEffect1.png"))
        hist(markerEffSel1, main = "Marker Effects for Trait1", 
             xlab = "Marker Effect")
        dev.off()
        
        png(filename = paste0(dirSaveCorEach, i, "_markerEffect2.png"))
        hist(markerEffSel2, main = "Marker Effects for Trait2", 
             xlab = "Marker Effect")
        dev.off()
        
        # save the marker effect
        markerEffMat <- matrix(0, nrow = ncol(genomeMatSel), ncol = nTrait)
        rownames(markerEffMat) <- colnames(genomeMatSel)
        markerEffMat[qtn1, 1] <- markerEffSel1
        markerEffMat[qtn2, 2] <- markerEffSel2
        
        write.csv(markerEffMat, file = paste0(dirSaveCorEach, i, "_markerEffects.csv"))
        if (all(!is.na(coeff))) {
          # converting the genetic correlation to non-linear
          u2 <- (qtnMat2 - 1) %*% markerEffSel2
          u1NonLinear <- coeff[1] * ((u2 - coeff[2]) ^ 2) / 4
          
          u1 <- ((qtnMat1 - 1) %*% markerEffSel1) + u1NonLinear
        } else {
          u1 <- (qtnMat1 - 1) %*% markerEffSel1
          u2 <- (qtnMat2 - 1) %*% markerEffSel2
        }
        
        png(filename = paste0(dirSaveCorEach, i, "_trueU.png"))
        
        plot(u1, u2, main = paste0("Genetic Correlation = ", 
                                   round(cor(u1, u2), 2)))
        dev.off()
        
        ##### Do a field trial and collect phenotypic data #####
        # Calculating the residual variance
        uMat <- cbind(u1, u2)
        sigmaG <- diag(var(uMat))* (nrow(uMat) - 1) / nrow(uMat)
        sigmaE <- sigmaG / h2 - sigmaG
        
        phenoMat <- GetPheno(uMat = uMat, 
                             sigmaE = sigmaE, 
                             nPheno = nPheno)
        
        initResList <- list(pheno = phenoMat,
                            uMat = uMat, 
                            genoMat = genomeMatSel, 
                            sigmaE = sigmaE)
        saveRDS(initResList, paste0(dirSaveCorEach, i, "_initResList.rds"))
        
        # save a genome data and phenotype data for visualizing the qtn data
        genoMat <- initPop$genoMat
        genoMat[genoMat == 0] <- "H"
        genoMat[genoMat == 2] <- "A"
        
        chr <- as.numeric(str_sub(map$chr, 4, 5))
        geno <- data.frame(chr, map$linkMapPos, t(genoMat))
        genoDf <- data.frame(id = c("", "", rownames(genoMat)), t(geno))
        write.csv(genoDf,
                  file = paste0(dirSaveCorEach, i, "_genoAll.csv"), row.names = F)
        
        pheno0 <- matrix(rnorm(2 * nrow(genoMat)), ncol = 2)
        pheno <- data.frame(id = rownames(genoMat), pheno0)
        colnames(pheno) <- c("id", "trait1", "trait2")
        rownames(pheno) <- rownames(genoMat)
        write.csv(pheno, file = paste0(dirSaveCorEach, i, "_phenoQTL.csv"))
        
        cross1 <- cross2 <- suppressWarnings(suppressMessages(
          read.cross(format = "csvs",
                     genfile = paste0(dirSaveCorEach, i, "_genoAll.csv"),
                     phefile = paste0(dirSaveCorEach, i, "_phenoQTL.csv")) |> 
            (\(x) {capture.output(x); x})()
        ))
        
        for (eachChrInd in 1:length(cross1$geno)) {
          # eachChrInd <- 1
          names(cross1$geno[[eachChrInd]]$map)[!(names(cross1$geno[[eachChrInd]]$map) %in% qtn1)] <- NA
          names(cross2$geno[[eachChrInd]]$map)[!(names(cross2$geno[[eachChrInd]]$map) %in% qtn2)] <- NA
        }
        
        # visualize the location of QTN1
        png(filename = paste0(dirSaveCorEach, i, "_QTN1s.png"))
        plotMap(cross1, col = "red", show.marker.names = T)
        dev.off()
        
        # visualize the location of QTN2
        png(filename = paste0(dirSaveCorEach, i, "_QTN2s.png"))
        plotMap(cross2, col = "green", show.marker.names = T)
        dev.off()
        
      }
    }
  }
}

# compute the mean variance of trait2
simName <- c("NoRelation/NoRelation", 
             "Pleiotropy/Positive", "Pleiotropy/NonLinear1_1", 
             "SpuriousPleiotropy/Positive", "SpuriousPleiotropy/NonLinear1_1")
resList <- paste0(dirSave, rep(simName, each = nRep), "/", 
                  rep(1:nRep, length(simName)), "_initResList.rds")

varList <- lapply(resList, function(resEachPath) {
  # resEachPath <- resList[[1]]
  res <- readRDS(resEachPath)
  uMat <- res$uMat
  return(var(uMat[, 2]) * (nrow(uMat)) / (nrow(uMat) - 1))
})
var2 <- mean(unlist(varList)) # 1.60 as the mean variance of trait2
qnorm(p = 0.75, mean = 0, sd = var2, lower.tail = T) # h = 1.08
qnorm(p = 0.25, mean = 0, sd = var2, lower.tail = T) # l = -1.08
