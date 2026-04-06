Sys.setenv(OPENBLAS_NUM_THREADS = "1")
Sys.setenv(OMP_NUM_THREADS      = "1")
Sys.setenv(MKL_NUM_THREADS      = "1")

ignore <- function(...) NULL

# compute the linkage disequilibrium (LD)
CalcLD <- function(genoMat, qtnEach, qtlCand) {
  snp <- c(qtnEach, qtlCand)
  genoSel <- genoMat[, snp]
  freq <- apply(genoSel, 2, function(x) {
    sum(x / 2) / length(x)
  })
  pAB <- (t(genoSel[, qtnEach]/2)) %*% (genoSel[, snp] / 2) / nrow(genoSel)
  D <- pAB - (freq[qtnEach] * t(freq))
  LD <- D^2 / ((freq[qtnEach] * t(freq)) * ((1 - freq[qtnEach]) * (1 - t(freq))))
  return(LD[2:length(LD)])
}

# Simulate the phenotype
GetPheno <- function(uMat, sigmaE, nPheno) {
  phenoList <- lapply(1:length(sigmaE), function(x) {
    residual <- mvrnorm(n = nPheno, mu = rep(0, nrow(uMat)), Sigma = diag(sigmaE[x], nrow = nrow(uMat)))
    phenoVec <- rep(uMat[, x], each = nPheno) + c(residual)
    pheno <- matrix(phenoVec, ncol = 1)
    rownames(pheno) <- rep(rownames(uMat), each = nPheno)
    return(pheno)
  })
  phenoMat <- do.call(cbind, phenoList)
  return(phenoMat)
}

# visualize the location of QTNs (modified plotMap in "qtl" package)
plotMap <- function (x, col = "green", map2, chr, horizontal = FALSE, shift = TRUE, show.marker.names = FALSE,
                     alternate.chrid = FALSE, ...)
{
  dots <- list(...)
  if ("main" %in% names(dots)) {
    themain <- dots$main
    usemaindefault <- FALSE
  }
  else usemaindefault <- TRUE
  if ("xlim" %in% names(dots)) {
    xlim <- dots$xlim
    usexlimdefault <- FALSE
  }
  else usexlimdefault <- TRUE
  if ("ylim" %in% names(dots)) {
    ylim <- dots$ylim
    useylimdefault <- FALSE
  }
  else useylimdefault <- TRUE
  if ("xlab" %in% names(dots))
    xlab <- dots$xlab
  else {
    if (horizontal)
      xlab <- "Location (cM)"
    else xlab <- "Chromosome"
  }
  if ("ylab" %in% names(dots))
    ylab <- dots$ylab
  else {
    if (horizontal)
      ylab <- "Chromosome"
    else ylab <- "Location (cM)"
  }
  map <- x
  if (inherits(map, "cross"))
    map <- pull.map(map)
  if (!missing(map2) && inherits(map2, "cross"))
    map2 <- pull.map(map2)
  if (!inherits(map, "map") || (!missing(map2) && !inherits(map2,
                                                            "map")))
    warning("Input should have class \"cross\" or \"map\".")
  if (!missing(map2) && is.matrix(map[[1]]) != is.matrix(map2[[1]]))
    stop("Maps must be both sex-specific or neither sex-specific.")
  if (!missing(chr)) {
    map <- map[matchchr(chr, names(map))]
    if (!missing(map2))
      map2 <- map2[matchchr(chr, names(map2))]
  }
  sex.sp <- FALSE
  if (is.matrix(map[[1]])) {
    one.map <- FALSE
    sex.sp <- TRUE
    if (!missing(map2)) {
      if (is.logical(map2)) {
        horizontal <- map2
        map2 <- lapply(map, function(a) a[2, ])
        map <- lapply(map, function(a) a[1, ])
      }
      else {
        Map1 <- lapply(map, function(a) a[1, , drop = TRUE])
        Map2 <- lapply(map, function(a) a[2, , drop = TRUE])
        Map3 <- lapply(map2, function(a) a[1, , drop = TRUE])
        Map4 <- lapply(map2, function(a) a[2, , drop = TRUE])
        old.mfrow <- par("mfrow")
        on.exit(par(mfrow = old.mfrow))
        par(mfrow = c(2, 1))
        class(Map1) <- class(Map2) <- class(Map3) <- class(Map4) <- "map"
        plotMap(Map1, Map3, horizontal = horizontal,
                shift = shift, show.marker.names = show.marker.names,
                alternate.chrid = alternate.chrid)
        plotMap(Map2, Map4, horizontal = horizontal,
                shift = shift, show.marker.names = show.marker.names,
                alternate.chrid = alternate.chrid)
        return(invisible(NULL))
      }
    }
    else {
      map2 <- lapply(map, function(a) a[2, ])
      map <- lapply(map, function(a) a[1, ])
    }
  }
  else {
    if (!missing(map2))
      one.map <- FALSE
    else one.map <- TRUE
  }
  if (one.map) {
    n.chr <- length(map)
    if (!show.marker.names) {
      chrpos <- 1:n.chr
      thelim <- range(chrpos) + c(-0.5, 0.5)
    }
    else {
      chrpos <- seq(1, n.chr * 2, by = 2)
      thelim <- range(chrpos) + c(-0.35, 2.35)
    }
    if (shift)
      map <- lapply(map, function(a) a - a[1])
    maxlen <- max(unlist(lapply(map, max)))
    if (horizontal) {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- c(0, maxlen)
      if (useylimdefault)
        ylim <- rev(thelim)
      plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
           yaxs = "i", xlab = xlab, ylab = ylab, yaxt = "n")
      a <- par("usr")
      for (i in 1:n.chr) {
        segments(min(map[[i]]), chrpos[i], max(map[[i]]),
                 chrpos[i])
        segments(map[[i]], chrpos[i] - 0.25, map[[i]],
                 chrpos[i] + 0.25)
        if (show.marker.names)
          text(map[[i]], chrpos[i] + 0.35, names(map[[i]]),
               srt = 90, adj = c(1, 0.5))
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 2,
                                            at = chrpos[i], labels = names(map)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    else {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- thelim
      if (useylimdefault)
        ylim <- c(maxlen, 0)
      plot(0, 0, type = "n", ylim = ylim, xlim = xlim,
           xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
      a <- par("usr")
      for (i in 1:n.chr) {
        segments(chrpos[i], min(map[[i]]), chrpos[i],
                 max(map[[i]]))
        segments(chrpos[i] - 0.25, map[[i]], chrpos[i] +
                   0.25, map[[i]])
        if (show.marker.names)
          if (any(!is.na(names(map[[i]]))))
            
            segments(chrpos[i] - 0.25, map[[i]][!(is.na(names(map[[i]])))], chrpos[i] +
                       0.25, map[[i]][!(is.na(names(map[[i]])))], col = col)
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 1,
                                            at = chrpos[i], labels = names(map)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    if (usemaindefault)
      title(main = "Genetic map")
    else if (themain != "")
      title(main = themain)
  }
  else {
    map1 <- map
    if (is.matrix(map2[[1]]))
      stop("Second map appears to be a sex-specific map.")
    if (length(map1) != length(map2))
      stop("Maps have different numbers of chromosomes.")
    if (any(names(map1) != names(map2))) {
      cat("Map1: ", names(map1), "\n")
      cat("Map2: ", names(map2), "\n")
      stop("Maps have different chromosome names.")
    }
    if (shift) {
      map1 <- lapply(map1, function(a) a - a[1])
      map2 <- lapply(map2, function(a) a - a[1])
    }
    n.mar1 <- sapply(map1, length)
    n.mar2 <- sapply(map2, length)
    markernames1 <- lapply(map1, names)
    markernames2 <- lapply(map2, names)
    if (any(n.mar1 != n.mar2)) {
      if (show.marker.names) {
        warning("Can't show marker names because of different numbers of markers.")
        show.marker.names <- FALSE
      }
    }
    else if (any(unlist(markernames1) != unlist(markernames2))) {
      if (show.marker.names) {
        warning("Can't show marker names because markers in different orders.")
        show.marker.names <- FALSE
      }
    }
    n.chr <- length(map1)
    maxloc <- max(c(unlist(lapply(map1, max)), unlist(lapply(map2,
                                                             max))))
    if (!show.marker.names) {
      chrpos <- 1:n.chr
      thelim <- range(chrpos) + c(-0.5, 0.5)
    }
    else {
      chrpos <- seq(1, n.chr * 2, by = 2)
      thelim <- range(chrpos) + c(-0.4, 2.4)
    }
    if (!horizontal) {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- thelim
      if (useylimdefault)
        ylim <- c(maxloc, 0)
      plot(0, 0, type = "n", ylim = ylim, xlim = xlim,
           xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
      a <- par("usr")
      for (i in 1:n.chr) {
        if (max(map2[[i]]) < max(map1[[i]]))
          map2[[i]] <- map2[[i]] + (max(map1[[i]]) -
                                      max(map2[[i]]))/2
        else map1[[i]] <- map1[[i]] + (max(map2[[i]]) -
                                         max(map1[[i]]))/2
        segments(chrpos[i] - 0.3, min(map1[[i]]), chrpos[i] -
                   0.3, max(map1[[i]]))
        segments(chrpos[i] + 0.3, min(map2[[i]]), chrpos[i] +
                   0.3, max(map2[[i]]))
        wh <- match(markernames1[[i]], markernames2[[i]])
        for (j in which(!is.na(wh))) segments(chrpos[i] -
                                                0.3, map1[[i]][j], chrpos[i] + 0.3, map2[[i]][wh[j]])
        if (any(is.na(wh)))
          segments(chrpos[i] - 0.4, map1[[i]][is.na(wh)],
                   chrpos[i] - 0.2, map1[[i]][is.na(wh)])
        wh <- match(markernames2[[i]], markernames1[[i]])
        if (any(is.na(wh)))
          segments(chrpos[i] + 0.4, map2[[i]][is.na(wh)],
                   chrpos[i] + 0.2, map2[[i]][is.na(wh)])
        if (show.marker.names)
          text(chrpos[i] + 0.35, map2[[i]], names(map2[[i]]),
               adj = c(0, 0.5))
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 1,
                                            at = chrpos[i], labels = names(map1)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map1)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 1, at = chrpos[i], labels = "")
          axis(side = 1, at = chrpos[i], labels = names(map1)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    else {
      old.xpd <- par("xpd")
      old.las <- par("las")
      par(xpd = TRUE, las = 1)
      on.exit(par(xpd = old.xpd, las = old.las))
      if (usexlimdefault)
        xlim <- c(0, maxloc)
      if (useylimdefault)
        ylim <- rev(thelim)
      plot(0, 0, type = "n", xlim = xlim, ylim = ylim,
           xlab = xlab, ylab = ylab, yaxt = "n", yaxs = "i")
      a <- par("usr")
      for (i in 1:n.chr) {
        if (max(map2[[i]]) < max(map1[[i]]))
          map2[[i]] <- map2[[i]] + (max(map1[[i]]) -
                                      max(map2[[i]]))/2
        else map1[[i]] <- map1[[i]] + (max(map2[[i]]) -
                                         max(map1[[i]]))/2
        segments(min(map1[[i]]), chrpos[i] - 0.3, max(map1[[i]]),
                 chrpos[[i]] - 0.3)
        segments(min(map2[[i]]), chrpos[i] + 0.3, max(map2[[i]]),
                 chrpos[[i]] + 0.3)
        wh <- match(markernames1[[i]], markernames2[[i]])
        for (j in which(!is.na(wh))) segments(map1[[i]][j],
                                              chrpos[i] - 0.3, map2[[i]][wh[j]], chrpos[i] +
                                                0.3)
        if (any(is.na(wh)))
          segments(map1[[i]][is.na(wh)], chrpos[i] -
                     0.4, map1[[i]][is.na(wh)], chrpos[i] - 0.2)
        wh <- match(markernames2[[i]], markernames1[[i]])
        if (any(is.na(wh)))
          segments(map2[[i]][is.na(wh)], chrpos[i] +
                     0.4, map2[[i]][is.na(wh)], chrpos[i] + 0.2)
        if (show.marker.names)
          text(map2[[i]], chrpos[i] + 0.35, names(map2[[i]]),
               srt = 90, adj = c(1, 0.5))
      }
      if (!alternate.chrid || length(chrpos) < 2) {
        for (i in seq(along = chrpos)) axis(side = 2,
                                            at = chrpos[i], labels = names(map1)[i])
      }
      else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map1)[i],
               line = -0.4, tick = FALSE)
        }
        for (i in even) {
          axis(side = 2, at = chrpos[i], labels = "")
          axis(side = 2, at = chrpos[i], labels = names(map1)[i],
               line = +0.4, tick = FALSE)
        }
      }
    }
    if (usemaindefault) {
      if (!sex.sp)
        title(main = "Comparison of genetic maps")
      else title(main = "Genetic map")
    }
    else if (themain != "")
      title(main = themain)
  }
  invisible()
}

GP <- function(phenoMat, genoMat) {
  amat <- calcGRM(genoMat = genoMat, methodGRM = "A.mat")
  Z <- design.Z(pheno.labels = rownames(phenoMat), geno.names = rownames(amat))
  ZETA <- list(A = list(Z = Z, K = amat))
  
  markerEffPredMat <- matrix(0, nrow = ncol(genoMat), ncol = ncol(phenoMat))
  uPredMat <- matrix(0, nrow = nrow(genoMat), ncol = ncol(phenoMat))
  rownames(markerEffPredMat) <- colnames(genoMat)
  rownames(uPredMat) <- rownames(genoMat)
  colnames(markerEffPredMat) <- colnames(uPredMat) <- colnames(phenoMat)
  beta <- matrix(NA, nrow = 1, ncol = ncol(phenoMat))
  
  for (traitInd in 1:ncol(phenoMat)) {
    # traitInd <- 1
    model <- EMM.cpp(y = phenoMat[, traitInd], ZETA = ZETA, n.core = 1)
    # h2 <- model$Vu / (model$Vu + model$Ve)
    
    # estimate marker effects from estimated genotypic values
    if (min(eigen(genoMat %*% t(genoMat))$values) < 1e-08) {
      mrkEffects <- t(genoMat) %*% solve(genoMat %*% t(genoMat) + diag(1e-06, nrow = nrow(genoMat))) %*% model$u
    } else {
      mrkEffects <- t(genoMat) %*% solve(genoMat %*% t(genoMat)) %*% model$u
    }
    mrkEffects <- c(mrkEffects)
    
    # input the result into matrix
    markerEffPredMat[, traitInd] <- mrkEffects
    uPredMat[, traitInd] <- model$u + model$beta
    beta[, traitInd] <- model$beta
  }
  return(list(uPredMat = uPredMat,
              markerEffPredMat = markerEffPredMat,
              betaEstimated = beta))
}

# calculating the D(2) matrix based on recombination frequency
CalcD2 <- function(ObjectMap, markerEffectMat, nTrait, k) {
  # creating the combinations among traits
  varInd <- matrix(rep(1:nTrait, each = 2), nrow = 2)
  covInd <- combn(x = 1:nTrait, m = 2)
  varCovInd <- cbind(varInd, covInd)
  
  # getting the number of chr
  nChr <- length(unique(ObjectMap$chr))
  
  # calculating the D2 for each combination
  D2List <- lapply(1:ncol(varCovInd), function(combInd) {
    # combInd <- 3
    comb <- varCovInd[, combInd]
    
    # separating to each chr
    D2EachChrList <- lapply(unique(ObjectMap$chr), function(eachChr) {
      # eachChr <- unique(ObjectMap$chr)[1]
      mapEach <- ObjectMap[ObjectMap$chr == eachChr, ]
      beta1 <- markerEffectMat[ObjectMap$chr == eachChr, comb[1]]
      beta2 <- markerEffectMat[ObjectMap$chr == eachChr, comb[2]]
      markerName <- rownames(mapEach)
      
      # calculating the distance of each marker set
      myDist <- sapply(1:nrow(mapEach),
                       function(x) abs(mapEach$linkMapPos[x] - mapEach$linkMapPos))
      
      # convert the distance to the recombination frequency (Haldane)
      c1 <- 0.5 * (1 - exp(-2 * (myDist / 100)))
      
      # r for the future F"k+1" generations
      ck <- 2*c1*(1 - ((0.5)^k)*(1 - 2*c1)^k) / (1 + 2*c1)
      
      D1_1 <- (1 - ck) * (1 - 2*c1) / 4
      D1_2 <- (1 - 2*ck - (0.5*(1 - 2*c1))^(k)) / 4
      
      D2_1 <- diag(beta1) %*% D1_1 %*% diag(beta2)
      D2_2 <- diag(beta1) %*% D1_2 %*% diag(beta2)
      
      rownames(D2_1) <- rownames(D2_2) <- markerName
      colnames(D2_1) <- colnames(D2_2) <- markerName
      return(D2 = list(D2_1, D2_2))
    })
    return(D2EachChrList)
    
  })
  return(D2List)
}

# Calculating the progeny variance for each cross
CalcVarMat <- function(D2, combInd, gametArray, nTrait) {
  sigmaMat <- sapply(D2, function(D2Each) {
    # D2Each <- D2[[1]]
    sigmaMatEachChr <- sapply(D2Each, function(D2EachChr) {
      # D2EachChr <- D2Each[[1]]
      D2_1_each <- D2EachChr[[1]]
      D2_2_each <- D2EachChr[[2]]
      markerName <- rownames(D2_1_each)
      
      markerMat_1 <- (gametArray[, markerName, 1] * 2) - 1
      markerMat_2 <- (gametArray[, markerName, 2] * 2) - 1
      
      mu1_1 <- diag(markerMat_1 %*% D2_1_each %*% t(markerMat_1)) / 4
      mu1_2 <- diag(markerMat_2 %*% D2_1_each %*% t(markerMat_2)) / 4
      gamma1 <- diag(markerMat_1 %*% D2_1_each %*% t(markerMat_2)) / 4
      gamma2 <- diag(markerMat_2 %*% D2_1_each %*% t(markerMat_1)) / 4
      
      gammaMat2_1 <- markerMat_1 %*% D2_2_each %*% t(markerMat_1) / 4
      gammaMat2_2 <- markerMat_2 %*% D2_2_each %*% t(markerMat_2) / 4
      mu2_1 <- diag(gammaMat2_1)
      mu2_2 <- diag(gammaMat2_2)
      gammaMat2_3 <- markerMat_1 %*% D2_2_each %*% t(markerMat_2) / 4
      gammaMat2_4 <- markerMat_2 %*% D2_2_each %*% t(markerMat_1) / 4
      
      p1 <- combInd[1, ]
      p2 <- combInd[2, ]
      
      eq1 <- mu1_1[p1] + mu1_2[p1] + mu1_1[p2] + mu1_2[p2]
      eq2 <- gamma1[p1] + gamma1[p2] + gamma2[p1] + gamma2[p2]
      eq3 <- mu2_1[p1] + mu2_2[p1] + mu2_1[p2] + mu2_2[p2]
      eq4 <- sapply(1:length(p1), function(x) {
        p1Each <- p1[x]
        p2Each <- p2[x]
        
        x1_1 <- gammaMat2_1[p1Each, p2Each]
        x1_2 <- gammaMat2_2[p1Each, p2Each]
        x1_3 <- gammaMat2_3[p1Each, p2Each]
        x1_4 <- gammaMat2_4[p1Each, p2Each]
        x2_1 <- gammaMat2_1[p2Each, p1Each]
        x2_2 <- gammaMat2_2[p2Each, p1Each]
        x2_3 <- gammaMat2_3[p2Each, p1Each]
        x2_4 <- gammaMat2_4[p2Each, p1Each]
        return(x1_1 + x1_2 + x1_3 + x1_4 + x2_1 + x2_2 + x2_3 + x2_4)
      })
      sigmaEach <- eq1 - eq2 + 2*eq3 - eq4
      return(sigmaEach)
    })
    sigmaVec <- apply(sigmaMatEachChr, 1, sum)
    return(sigmaVec)
  })
  return(sigmaMat)
}

CalcProgMean <- function(predU) {
  predProgMeanMat <- apply(predU, 2, function(predUeach) {
    predProgMeanAll <- (rep(predUeach, each = length(predUeach)) + rep(predUeach, length(predUeach))) / 2
    predProgMeanMat <- matrix(predProgMeanAll, nrow = length(predUeach), ncol = length(predUeach))
    rownames(predProgMeanMat) <- colnames(predProgMeanMat) <- names(predUeach)
    predProgMean <- predProgMeanMat[lower.tri(predProgMeanMat)]
    return(predProgMean)
  })
  return(predProgMeanMat)
  
}

ExtractGamet <- function(nowPop, indName, qtn) {
  gametArrayList <- lapply(nowPop$inds, function(x) {
    # x <- initPop$inds[[1]]
    gamet <- do.call(cbind, x$haplo$values)
    gamet <- gamet[, qtn]
    array(data = c(t(gamet)), dim = c(1, ncol(gamet), 2))
  })
  gametArrayAll <- do.call(abind, c(gametArrayList, along = 1))
  gametArray <- gametArrayAll[indName, , ]
  return(gametArray)
}


# Extracting the best crosses using Linear Programing(LP)
LP <- function(predProgMat, # the potential of each cross
               indName, # the name of individual
               comb, # the matrix of each cross
               capa, # capacity for each edge (individual)
               nCrosses, # the number of crosses
               nProg, # the number of progenies for each cross
               generation # count of the generation
) {
  # matrix to vector
  dataVec <- c(-predProgMat)
  crossNameSel <- rownames(predProgMat)
  
  # edges
  ed1 <- paste0("L_", indName)
  ed2 <- paste0("R_", indName)
  ed <- c(ed1, ed2)
  
  # nodes
  s_es <- paste0("S-", ed1)
  t_es <- paste0(ed2, "-T")
  c_es <- paste0("S-L_", comb[1, ],
                 "-R_",
                 comb[2, ], "-T")
  es <- c(s_es, t_es, c_es)
  
  # cost
  s_cs <- rep(0, length(s_es))
  t_cs <- rep(0, length(t_es))
  c_cs <- dataVec
  cs <- c(s_cs, t_cs, c_cs)
  
  # limitations
  # Normally, we make a design matrix, but it will be so large.
  # we made a M times 3 matrix (M depends on the limitation)
  # first column means the location of row
  # second column means the location of column
  # third column means the value of limitation
  stc_lim <- matrix(c(rep(1:length(cs), 2),
                      rep(1, length(cs))), ncol = 3)
  s_lim <- matrix(c(rep(length(cs) + 1, length(s_cs)),
                    1:length(s_cs),
                    rep(1, length(s_cs))), ncol = 3)
  t_lim <- matrix(c(rep(length(cs) + 2, length(s_cs)),
                    (length(s_cs)+1):(length(s_cs) + length(t_cs)),
                    rep(1, length(s_cs))), ncol = 3)
  t_lim_end <- t_lim[1, 1]
  
  # the amount of flow of input and output from each edge is must be the same
  node_num <- length(s_es) + length(t_es) + 2*length(c_es)
  ed_lim <- matrix(NA,
                   nrow = node_num,
                   ncol = 3)
  for (i in 1:length(ed)) {
    # i <- 1
    # extracting each edge
    ed_each <- ed[i]
    ed_lim_ind <- t_lim_end + i
    
    # index for searching the nodes (edge to edge) which come from "ed_each"
    c_ind <- paste0("-", ed_each, "-")
    
    # start to edge (edge to end) (-1)
    start <- min(which(is.na(ed_lim[, 1])))
    ed_lim[start, ] <- c(ed_lim_ind, i, -1)
    
    # edge to edge (1)
    ed_lim_each <- which(grepl(pattern = c_ind, x = es))
    if (length(ed_lim_each) < 1) {
      next
    }
    
    # add the limitation (the amount of flow is the same)
    ed_lim_mat <- matrix(c(rep(ed_lim_ind, length(ed_lim_each)),
                           ed_lim_each,
                           rep(1, length(ed_lim_each))), ncol = 3)
    start <- min(which(is.na(ed_lim[, 1])))
    ed_lim[start:(start + length(ed_lim_each) - 1), ] <- ed_lim_mat
  }
  
  ed_lim_end <- ed_lim[nrow(ed_lim), 1] + 1
  
  # the total amount of flow of specific edge
  # (start to the edge & the edge to end)
  # must be limited to "capa"
  sym_lim_row_ind <- ed_lim_end:(ed_lim_end + length(s_cs) - 1)
  sym_lim_start <- matrix(c(sym_lim_row_ind,
                            1:length(s_cs),
                            rep(1, length(s_cs))), ncol = 3)
  sym_lim_end_ind <- (1+length(s_cs)):(length(s_cs)+length(t_cs))
  sym_lim_end <- matrix(c(sym_lim_row_ind,
                          sym_lim_end_ind,
                          rep(1, length(t_cs))), ncol = 3)
  sym_lim <- rbind(sym_lim_start, sym_lim_end)
  lim <- rbind(stc_lim, s_lim, t_lim, ed_lim, sym_lim)
  
  # lim <- rbind(stc_lim, s_lim, t_lim, ed_lim, sym_lim)
  eq <- c(rep("<=", length(c(s_cs, t_cs))),
          rep("<=", length(c_cs)),
          rep("==", 2 + length(ed)),
          rep("<=", length(s_cs)))
  
  # capacity
  s_ws <- rep(capa, length(s_es))
  t_ws <- rep(capa, length(t_es))
  c_ws <- rep(1, length(c_es))
  flow <- rep(nCrosses, 2)
  ed_ws <- rep(0, length(ed))
  sym_ws <- rep(capa, length(s_cs))
  ws <- c(s_ws, t_ws, c_ws, flow, ed_ws, sym_ws)
  res <- lp(direction = "min",
            objective.in = cs,
            dense.const = lim,
            const.dir = eq,
            const.rhs = ws,
            all.int = T)
  ind <- as.character(round(res$solution))
  # sum(dataVec[ind[(2 * n + 1):length(cs)] == 1])
  selCrosses0 <- es[ind == 1]
  selCrosses0 <- selCrosses0[grepl(pattern = "^S.{1,}-T$", x = selCrosses0)]
  selCrosses <- str_sub(selCrosses0, 5, -3)
  selCrossesList <- str_split(selCrosses, pattern = "-R_")
  bestCrossMat <- do.call(rbind, selCrossesList)
  if (max(table(c(bestCrossMat))) > 2) {
    print("Each genotype must not be used more than nCrosses times")
  } else {
    crossTable <- data.frame(ind1 = bestCrossMat[, 1],
                             ind2 = bestCrossMat[, 2],
                             n = nProg,
                             names = paste0("C", generation, "N",
                                            formatC(1:nrow(bestCrossMat), width = 3, flag = "0")))
    return(crossTable)
  }
}

CreateF8 <- function(selfSelected, nowPop, nProgSelf, markerEffectEstimated, markerEffectTrue, marker, beta, sigmaE, nPheno, nTrait, coeff, generation) {
  
  f8GenerationInd <- str_split(selfSelected[1], pattern = "N")[[1]][1]
  f8Generation <- as.numeric(str_sub(f8GenerationInd, 2)) + 8
  f8Generation <- formatC(f8Generation, width = 2, flag = "0")
  selfTable <- data.frame(ind1 = selfSelected,
                          ind2 = selfSelected,
                          n = nProgSelf,
                          names = paste0("C", f8Generation, "F2N",
                                         formatC(1:length(selfSelected), width = 3, flag = "0")))
  f2Pop <- population$new(name = "F2 offspring",
                          inds = makeCrosses(crosses = selfTable, pop = nowPop))
  selfNowPop <- f2Pop
  for (l in 2:7) {
    selfTable <- data.frame(ind1 = rownames(selfNowPop$genoMat),
                            ind2 = rownames(selfNowPop$genoMat),
                            n = 1,
                            names = paste0("C", f8Generation, "F", l+1, "N",
                                           formatC(1:nrow(selfNowPop$genoMat), width = 3, flag = "0")))
    selfNextPop <- population$new(name = paste0("F", l+1, " offspring"),
                                  inds = makeCrosses(crosses = selfTable, pop = selfNowPop))
    selfNowPop <- selfNextPop
  }
  genoMatNew <- selfNowPop$genoMat
  betaMat <- matrix(beta, nrow = nrow(genoMatNew), ncol = nTrait, byrow = T)
  geneticValEstimatedNew <- (genoMatNew[, marker] - 1) %*% markerEffectEstimated + betaMat
  geneticValTrueNew <- (genoMatNew - 1) %*% markerEffectTrue
  
  # Computing bulmer effect
  AlleleFreq <- apply(genoMatNew, 2, function(eachAllele) {
    p <- sum(eachAllele) / (2 * length(eachAllele))
    return(4 * p * (1 - p))
  })
  genicVar <- apply(markerEffectTrue, 2, function(x) {
    return(AlleleFreq %*% x^2)
  })
  geneticVar <- apply(geneticValTrueNew, 2, function(x) {
    var(x) * (length(x) - 1) / length(x)
  })
  geneticCov <- var(geneticValTrueNew[, 1], geneticValTrueNew[, 2]) * (nrow(geneticValTrueNew) - 1) / nrow(geneticValTrueNew)
  bulmer <- geneticVar / genicVar
  
  if (all(!is.na(coeff))) {
    # converting the genetic correlation to non-linear
    u1NonLinear <- coeff[1] * ((geneticValTrueNew[, 2] - coeff[2]) ^ 2) / 4
    geneticValTrueNew[, 1] <- geneticValTrueNew[, 1] + u1NonLinear
  } 
  
  set.seed(geneticValTrueNew[1, 1])
  phenoNew <- GetPheno(uMat = geneticValTrueNew, 
                       sigmaE = sigmaE, 
                       nPheno = nPheno)
  return(list(pheno = phenoNew,
              genoMat = genoMatNew,
              predU = geneticValEstimatedNew,
              trueU = geneticValTrueNew,
              genicVar = genicVar, 
              geneticVar = geneticVar, 
              geneticCov = geneticCov, 
              bulmer = bulmer, 
              generation = generation))
}

# visualize the prediction accuracy of genotypic values
FigPredU <- function(gpPath, trueU, predU) {
  png(gpPath, width = 1440, height = 720)
  par(mfrow = c(1, 2))
  xlim <- ylim <- range(trueU[, 1], predU[, 1])
  plot(trueU[, 1], predU[, 1], 
       xlab = "True Genotypic Value of Trait1", ylab = "Predicted Genotypic Value of Trait1", 
       xlim = xlim, ylim = ylim, 
       main = paste0("r = ", round(cor(trueU[, 1], predU[, 1]), 2)))
  abline(a = 0, b = 1, lty = 2, col = "red")
  
  xlim <- ylim <- range(trueU[, 2], predU[, 2])
  plot(trueU[, 2], predU[, 2], 
       xlab = "True Genotypic Value of Trait2", ylab = "Predicted Genotypic Value of Trait2", 
       xlim = xlim, ylim = ylim, 
       main = paste0("r = ", round(cor(trueU[, 2], predU[, 2]), 2)))
  abline(a = 0, b = 1, lty = 2, col = "red")
  dev.off()
}

# draw a figure of true genotypic values
FigTrueU <- function(plotPath, trueU, h, l, generation) {
  png(plotPath)
  plot(
    trueU,
    xlab = "Trait1",
    ylab = "Trait2",
    xlim = c(-5, 40),
    ylim = c(-15, 15),
    main = paste0("Genetic Value of Generation ", generation)
  )
  abline(h = c(h, l), lty = 2, col = "red")
  dev.off()
}


# draw a figure of marker effects prediction accuracy
FigPredMarker <- function(markerPredictPath, 
                          markerEffectEstimatedMat, 
                          markerEffectTrueMat, 
                          alleleFixed, 
                          corMarker1, 
                          corMarker2, 
                          nMarker, 
                          nFixedAllele) {
  png(markerPredictPath, width = 1440, height = 720)
  par(mfrow = c(1, 2))
  xlim <- ylim <- range(markerEffectEstimatedMat[, 1],
                        markerEffectTrueMat[, 1])
  plot(markerEffectTrueMat[!alleleFixed, 1],
       markerEffectEstimatedMat[!alleleFixed, 1],
       xlab = "True Marker Effects of Trait1", ylab = "Estimated Marker Effects of Trait1",
       xlim = xlim, ylim = ylim,
       main = paste0("r = ", round(corMarker1, 2),
                     ", nMarker = ", nMarker - nFixedAllele))
  abline(a = 0, b = 1, lty = 2, col = "red")
  
  xlim <- ylim <- range(markerEffectEstimatedMat[, 2],
                        markerEffectTrueMat[, 2])
  plot(markerEffectTrueMat[!alleleFixed, 2],
       markerEffectEstimatedMat[!alleleFixed, 2],
       xlab = "True Marker Effects of Trait2", ylab = "Estimated Marker Effects of Trait2",
       xlim = xlim, ylim = ylim,
       main = paste0("r = ", round(corMarker2, 2),
                     ", nMarker = ", nMarker - nFixedAllele))
  abline(a = 0, b = 1, lty = 2, col = "red")
  dev.off()
}

# draw a figure of estimated genetic variance
FigPredGenVar <- function(plotPathVar1, 
                          plotPathVar2, 
                          plotPathCov, 
                          sigmaMat, 
                          sigmaTrueMat, 
                          rmseVar1, 
                          rmseVar2, 
                          rmseCov) {
  
  xlim <- ylim <- range(c(sigmaMat[, 1], sigmaTrueMat[, 1]))
  png(plotPathVar1)
  plot(sigmaTrueMat[, 1], sigmaMat[, 1],
       xlab = "Computed", ylab = "Estimated",
       xlim = xlim, ylim = ylim,
       main = paste0("gVar1 generation ", generation,
                     ", RMSE = ", round(rmseVar1, 3)))
  dev.off()
  
  xlim <- ylim <- range(c(sigmaMat[, 2], sigmaTrueMat[, 2]))
  png(plotPathVar2)
  plot(sigmaTrueMat[, 2], sigmaMat[, 2],
       xlab = "Computed", ylab = "Estimated",
       xlim = xlim, ylim = ylim,
       main = paste0("gVar2 generation ", generation,
                     ", RMSE = ", round(rmseVar2, 3)))
  dev.off()
  
  xlim <- ylim <- range(c(sigmaMat[, 3], sigmaTrueMat[, 3]))
  png(plotPathCov)
  plot(sigmaTrueMat[, 3], sigmaMat[, 3],
       xlab = "Computed", ylab = "Estimated",
       xlim = xlim, ylim = ylim,
       main = paste0("gCov generation ", generation,
                     ", RMSE = ", round(rmseCov, 3)))
  dev.off()
}

SimSum <- function(
  pathList,
  target,
  baseFileName,
  index,
  generationInd,
  methodFactor,
  typeFactor,
  corFactor,
  corColor
) {
  resDfList <- lapply(pathList, function(pathListEach) {
    # pathListEach <- pathList[[1]]
    method <- pathListEach$method
    resList <- lapply(pathListEach$path, function(path) {
      # path <- pathListEach$path[[1]]
      nameInd <- length(str_split(path, pattern = "/")[[1]])
      typeFactorEach <- str_split(path, pattern = "/")[[1]][nameInd - 2]
      corFactorEach <- str_split(path, pattern = "/")[[1]][nameInd - 1]
      # sort the file based on file name
      fileName <- list.files(path, pattern = baseFileName)
      fileNum <- as.numeric(str_sub(
        fileName,
        start = 1,
        end = -nchar(fileName[1])
      ))
      fileOrder <- order(fileNum)

      # gathering the results files of "F4" and "Recurrent genomic selection"
      resList <- list.files(path, pattern = baseFileName, full.names = T)[
        fileOrder
      ]
      resArray <- array(
        data = NA,
        dim = c(length(generationInd), length(index), nRep)
      )
      dimnames(resArray) <- list(generationInd, index, 1:nRep)

      for (i in 1:nRep) {
        # i <- 1
        resMat0 <- as.matrix(read.csv(resList[[i]], row.names = 1))
        resMat0 <- resMat0[1:dim(resArray)[1], index, drop = F]
        # resMat0[, c("top", "top10Mean", "effectiveTop", "effectiveTop10Mean")] <- resMat0[, c("top", "top10Mean", "effectiveTop", "effectiveTop10Mean")] + 10

        for (l in 1:length(index)) {
          traitName <- index[l]
          if (traitName == "effectiveTop") {
            resArray[, traitName, i] <- (resMat0[, traitName] -
              resMat0[1, traitName]) /
              sqrt(resMat0[1, "geneticVar1"])
          } else if (traitName == "nFixedAllele") {
            resArray[, traitName, i] <- 4000 - resMat0[, traitName]
          } else {
            resArray[, traitName, i] <- resMat0[, traitName]
          }
        }
      }

      resMeanMat <- apply(resArray, c(1, 2), mean, na.rm = T)
      resSeMat <- apply(resArray, c(1, 2), std.error, na.rm = T)

      generation <- as.numeric(str_sub(rownames(resMeanMat), start = 2))
      generationRep <- rep(generation, ncol(resMeanMat))
      methodRep <- rep(method, length(generationRep))
      indexRep <- rep(colnames(resMeanMat), each = nrow(resMeanMat))
      typeFactorRep <- rep(typeFactorEach, ncol(resMeanMat))
      corFactorRep <- rep(corFactorEach, ncol(resMeanMat))
      dfEach <- data.frame(
        Generation = generationRep,
        Method = methodRep,
        Index = indexRep,
        Type = typeFactorRep,
        Factor = corFactorRep,
        Value = c(resMeanMat),
        SE = c(resSeMat)
      )
      return(dfEach)
    })
    resDfListEach <- do.call(rbind, resList)
    return(resDfListEach)
  })
  resDf <- do.call(rbind, resDfList)
  for (i in 1:length(index)) {
    # i <- 1
    resIndEach0 <- index[i]
    resDfEach <- resDf[resDf$Index == resIndEach0, ]
    resDfEach$Method <- factor(resDfEach$Method, levels = methodFactor)
    resDfEach$Type <- factor(resDfEach$Type, levels = unique(typeFactor))
    resDfEach$Method <- factor(resDfEach$Method, levels = methodFactor)
    resDfEach$Factor <- factor(resDfEach$Factor, levels = unique(corFactor))
    col <- corColor

    if (grepl("^accuracy.{3}", resIndEach0)) {
      res0 <- resDfEach[resDfEach$Generation == 0, "Value"]
      res20 <- resDfEach[resDfEach$Generation == 20, "Value"]
      improve <- res20 / res0
      cat(paste0(
        "CPSMT Improvement Rate of ",
        resIndEach0,
        " of 3 Scenarios = ",
        mean(improve),
        " \n"
      ))
    }

    if (grepl("^rmse.*", resIndEach0)) {
      res1 <- resDfEach[resDfEach$Generation == 1, "Value"]
      res20 <- resDfEach[resDfEach$Generation == 20, "Value"]
      pr <- (res1 - res20) / res1 * 100
      cat(paste0(
        "CPSMT Reduction Rate of ",
        resIndEach0,
        " of 3 Scenarios = ",
        mean(pr),
        " %\n"
      ))
    }

    if (grepl("^corMarker.*", resIndEach0)) {
      res1 <- resDfEach[resDfEach$Generation == 1, "Value"]
      res20 <- resDfEach[resDfEach$Generation == 20, "Value"]
      improve <- res20 / res1
      cat(paste0(
        "CPSMT Improvement Rate of ",
        resIndEach0,
        " of 3 Scenarios = ",
        mean(improve),
        " \n"
      ))
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

    for (l in 1:length(unique(typeFactor))) {
      # l <- 3
      typeEach <- unique(typeFactor)[l]
      resDfSel <- resDfEach[resDfEach$Type == typeEach, ]
      g <- ggplot(resDfSel, aes(x = Generation, y = Value, color = Method)) +
        geom_line(size = 1.5) +
        facet_grid(Type ~ Factor) +
        scale_color_manual(values = col) +
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

      if (length(unique(resDfEach$Factor)) > 2) {
        png(
          paste0(dirSave, target, "_", resIndEach0, "_", typeEach, ".png"),
          height = 720,
          width = 2880,
          res = 214
        )
      } else {
        png(
          paste0(dirSave, target, "_", resIndEach0, "_", typeEach, ".png"),
          height = 720,
          width = 960,
          res = 214
        )
      }
      print(g)
      dev.off()
    }

    if (length(pathList) > 1) {
      # calculating the proportion that LP exceed GS
      pattern <- paste0(typeFactor, "_", corFactor)
      patternDf <- paste0(resDfEach$Type, "_", resDfEach$Factor)
      propList <- lapply(pattern, function(patternEach) {
        # patternEach <- pattern[1]
        newMethodVal <- resDfEach[
          (resDfEach$Method == "CPSMT" & patternDf == patternEach),
          "Value"
        ]
        oldMethodVal <- resDfEach[
          (resDfEach$Method == "PCS" & patternDf == patternEach),
          "Value"
        ]
        proportion <- newMethodVal / oldMethodVal
        return(proportion)
      })
      propMat <- do.call(rbind, propList)
      rownames(propMat) <- pattern
      write.csv(propMat, paste0(dirSave, target, "_", resIndEach0, "_Prop.csv"))

      # compute remained gentic variance
      if (grepl("geneticVar", resIndEach0)) {
        varInit <- resDfEach[
          (resDfEach$Method == "PCS" & resDfEach$Generation == 0),
        ]
        varFinalCpsmt <- resDfEach[
          (resDfEach$Method == "CPSMT" & resDfEach$Generation == 20),
        ]
        varFinalPcs <- resDfEach[
          (resDfEach$Method == "PCS" & resDfEach$Generation == 20),
        ]

        patternInit <- paste0(varInit$Type, "_", varInit$Factor)
        patternCpsmt <- paste0(varFinalCpsmt$Type, "_", varFinalCpsmt$Factor)
        patternPcs <- paste0(varFinalPcs$Type, "_", varFinalPcs$Factor)

        if (all(c(patternInit == patternCpsmt, patternInit == patternPcs))) {
          varRemainedCpsmt <- varFinalCpsmt[, "Value"] / varInit[, "Value"]
          varRemainedPcs <- varFinalPcs[, "Value"] / varInit[, "Value"]
          cat(paste0(
            "CPSMT remained ",
            resIndEach0,
            " = ",
            mean(varRemainedCpsmt) * 100,
            " %\n"
          ))
          cat(paste0(
            "PCS remained ",
            resIndEach0,
            " = ",
            mean(varRemainedPcs) * 100,
            " %\n"
          ))
        }
      }
      # compute prediction accuracy of the final generation
      if (grepl("accuracy", resIndEach0)) {
        accFinalCpsmt <- resDfEach[
          (resDfEach$Method == "CPSMT" & resDfEach$Generation == 20),
        ]
        accFinalPcs <- resDfEach[
          (resDfEach$Method == "PCS" & resDfEach$Generation == 20),
        ]

        excludeFactors <- c("NonLinear1_2", "Positive")
        accFinalSelCpsmt <- accFinalCpsmt[
          !(accFinalCpsmt$Factor %in% excludeFactors),
          "Value"
        ]
        accFinalSelPcs <- accFinalPcs[
          !(accFinalPcs$Factor %in% excludeFactors),
          "Value"
        ]

        cat(paste0(
          "CPSMT Pred ",
          resIndEach0,
          " of 3 Scenarios in Final = ",
          round(mean(accFinalSelCpsmt), 3),
          "\n"
        ))
        cat(paste0(
          "PCS Pred ",
          resIndEach0,
          " of 3 Scenarios in Final = ",
          round(mean(accFinalSelPcs), 3),
          "\n"
        ))
      }
    }
  }
}

# Monte Carlo evaluation of penalized cross/self potential
EvalByMC <- function(
  meanMat,
  varMat,
  w,
  l,
  h,
  prob,
  nMC,
  meanPredU,
  sdPredU,
  hardConstraint = FALSE
) {
  z <- matrix(rnorm(nMC * 2), nrow = nMC, ncol = 2)
  mcResult <- sapply(1:nrow(meanMat), function(j) {
    muEach <- meanMat[j, ]
    sigmaMat_j <- matrix(
      c(varMat[j, 1], varMat[j, 3], varMat[j, 3], varMat[j, 2]),
      nrow = 2
    )

    # if variance is nearly zero, return penalized mean directly
    if (sigmaMat_j[1, 1] < 1e-6 & sigmaMat_j[2, 2] < 1e-6) {
      inRange <- muEach[2] >= l & muEach[2] <= h
      if (hardConstraint) {
        if (!inRange) {
          return(-10000)
        }
        return((muEach[1] - meanPredU[1]) / sdPredU[1])
      }
      penalty <- ifelse(
        inRange,
        0,
        max((muEach[2] - h) / sdPredU[2], (l - muEach[2]) / sdPredU[2])
      )
      return(w * (muEach[1] - meanPredU[1]) / sdPredU[1] - (1 - w) * penalty)
    }

    # ensure positive definiteness
    sigmaMat_j[1, 1] <- max(sigmaMat_j[1, 1], 1e-6)
    sigmaMat_j[2, 2] <- max(sigmaMat_j[2, 2], 1e-6)
    maxCov <- sqrt(sigmaMat_j[1, 1] * sigmaMat_j[2, 2]) * (1 - 1e-6)
    sigmaMat_j[1, 2] <- max(min(sigmaMat_j[1, 2], maxCov), -maxCov)
    sigmaMat_j[2, 1] <- sigmaMat_j[1, 2]

    L <- chol(sigmaMat_j)
    x <- z %*% L
    x[, 1] <- x[, 1] + muEach[1]
    x[, 2] <- x[, 2] + muEach[2]

    if (hardConstraint) {
      inRange <- x[, 2] >= l & x[, 2] <= h
      x1hard <- ifelse(inRange, x[, 1], -10000)
      threshold <- quantile(x1hard, probs = prob)
      return(min(x1hard[x1hard >= threshold]))
    }

    penalty <- ifelse(
      x[, 2] >= l & x[, 2] <= h,
      0,
      pmax((x[, 2] - h) / sdPredU[2], (l - x[, 2]) / sdPredU[2])
    )
    penalizedVal <- w * (x[, 1] - meanPredU[1]) / sdPredU[1] - (1 - w) * penalty
    threshold <- quantile(penalizedVal, probs = prob)
    min(penalizedVal[penalizedVal >= threshold])
  })
  names(mcResult) <- rownames(meanMat)
  return(mcResult)
}
