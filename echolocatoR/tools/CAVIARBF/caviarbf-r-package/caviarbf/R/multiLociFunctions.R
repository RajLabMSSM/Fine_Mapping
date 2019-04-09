calculateHessianOfNegLik = function(annotation, gamma, J, lambda,
                                    fittingMethod = NULL, alpha = NULL) {
  # Hessian is 0 for lasso penalty right now
  # This needs to be further justified or modified
  
  nGamma = dim(annotation)[2]
  if (is.null(fittingMethod) || 
      fittingMethod == "PCEMin" ||
      fittingMethod == "glmnetRidgeMin") {
      S0Inverse = diag(c(0, rep(lambda, nGamma - 1)), nrow = nGamma)
  } else if (fittingMethod == "glmnetLASSOMin") {
    S0Inverse = 0
  } else if (fittingMethod == "glmnetENETMin") {
    if (is.null(alpha)) {
      stop("alpha should not be null")
    }
    lambda2 = lambda * (1 - alpha)
    S0Inverse = diag(c(0, rep(lambda2, nGamma - 1)), nrow = nGamma)
  }
#   print("S0Inverse")
#   print(S0Inverse)
  y_hat = 1 / (1 + exp(-(annotation %*% gamma)))
  w = c(y_hat * (1 - y_hat))
#   print("XTWX")
#   print(t(annotation) %*% (annotation * w))
  secondDerivOfNegLogLik = t(annotation) %*% (annotation * w) + S0Inverse
  if (!is.null(J)) {
    # make sure it is symmetric
    XTJ = t(annotation) %*% J
    XTJ = 0.5 * (XTJ + t(XTJ))
#     print("XTJ")
#     print(XTJ)
    secondDerivOfNegLogLik = secondDerivOfNegLogLik - XTJ
  }
  
  # make sure all eigenvalues is positive
  eigenResult = eigen(secondDerivOfNegLogLik, symmetric = T)
  eigenvalues = eigenResult$values
  if (any(eigenvalues < 0)) {
    print("eigenvalues of secondDerivOfNegLogLik (before fix):")
    print(eigenvalues)
    
    # adhoc fix. a better way is to derive the analytical version 
    # or better approximation
    indNegative = which(eigenvalues <= 0)
    indPositive = which(eigenvalues > 0)
    if (is.null(indPositive)) {
      stop("bad Hessian matrix")
    }
    minPositive = min(eigenvalues[indPositive])
    eigenvalues[indNegative] = minPositive
    
    print("eigenvalues of secondDerivOfNegLogLik (after fix):")
    print(eigenvalues)
  }
  v = eigenResult$vectors
  secondDerivOfNegLogLik = v %*% (t(v) * eigenvalues)
  
#   print("Hessian of the negative likelihood")
#   print(secondDerivOfNegLogLik)
  return(secondDerivOfNegLogLik)
}

predictivePrior = function(annotation, gamma, A, annotationNew = NULL) {
  # average over the distribution, even though this is not the correct way
  # in the Bayesian framework
  
  if (is.null(annotationNew)) {
    annotationNew = annotation
  }
  mu_a = c(annotationNew %*% gamma)
  nNew = dim(annotationNew)[1]
  sigma_a2 = numeric(nNew)
  for (i in 1:nNew) {
    sigma_a2[i] = t(annotationNew[i, ]) %*% ginv(A) %*% annotationNew[i, ]
  }
  k_sigma_a2 = (1 + 3.1415926 / 8 * sigma_a2)^(-1/2)
  priorAveraged = 1 / (1 + exp(-(k_sigma_a2 * mu_a)))
  return(priorAveraged)
}

aggregateAnnotation = function(lociListFile, 
                               annotationSuffix = ".annotations") {
  filePrefixes = as.vector(as.matrix(read.table(lociListFile, as.is = T)))
  dataDir = dirname(lociListFile)
  annotation = NULL
  lociID = NULL
  for (i in 1:length(filePrefixes)) {
    annotationFile = paste0(dataDir, "/", filePrefixes[i], annotationSuffix)
    annotationPerFile = as.matrix(read.table(annotationFile))
    lociIDPerFile = rep(i, dim(annotationPerFile)[1])
    annotation = rbind(annotation, annotationPerFile)
    lociID = c(lociID, lociIDPerFile)
  }
  return(list(annotation = annotation, lociID = lociID))
}

connectPIPFileToZFile = function(lociListFile, PIPFile, PIPZFile, 
                                 hasHeader = F, zFileSuffix = "") {
  filePrefixes = as.vector(as.matrix(read.table(lociListFile, as.is = T)))
  dataDir = dirname(lociListFile)
  lociID = NULL
  zInfo = NULL
  for (i in 1:length(filePrefixes)) {
    zFile = paste0(dataDir, "/", filePrefixes[i], zFileSuffix)
    zPerFile = as.matrix(read.table(zFile, header = hasHeader))
    zInfo = rbind(zInfo, zPerFile)
    lociIDPerFile = rep(i, dim(zPerFile)[1])
    lociID = c(lociID, lociIDPerFile)
  }
  PIPs = read.table(PIPFile)
  indices = PIPs[, 1]
  PIPZ = cbind(PIPs, lociID[indices], zInfo[indices, ])
  write.table(PIPZ, file = PIPZFile, row.names = F, col.names = F, quote = F)
}

# CAVIARBFWrapper = function(zFile, rFile, priorType01, priorValue, nSample, 
#                            maxCausal, BFFile, exact = T, eps = 0,
#                            useIdentityMatrix = F) {
#   # A wrapper of the C++ program caviarbf
#   
#   exactString = ""
#   if (!exact) {
#     exactString = "--appr"
#   }
#   identityString = ""
#   if (useIdentityMatrix) {
#     identityString = "-i"
#   } 
#   CAVIARBFExe = normalizePath("../../../tools/caviarbf/caviarbf")
#   consoleOutputFile = paste0(BFFile, ".console")
#   cmd = paste(CAVIARBFExe, "-z", zFile, "-r", rFile, "-t", priorType01, 
#               "-a", paste(priorValue, collapse = ","), 
#               "-n", nSample, "-c", maxCausal, 
#               exactString, "-e", eps, identityString,
#               "-o", BFFile, ">>", consoleOutputFile)
#   cat(cmd, "\n", file = consoleOutputFile)
#   stopifnot(system(cmd) == 0)
# }
# 
# modelSearchWrapper = function(BFFile, outputFilePrefix, nSNPs, priorFile) {
#   modelSearchExe = normalizePath("../../../tools/caviarbf/model_search")
#   consoleOutputFile = paste0(BFFile, ".model_search.console")
#   cmd = paste(modelSearchExe, "-i", BFFile, "-o", outputFilePrefix,
#               "-m", nSNPs, "-f", priorFile, ">>", consoleOutputFile)
#   cat("\n\n== command executed:\n", cmd, "\n", file = consoleOutputFile,
#       append = T)
#   stopifnot(system(cmd) == 0)  
# }

# calculatePIPsFromBFs = function(PIPsFromAnnotations, lociID, BFFileList, 
#                                 iter) {
#   # browser()
#   PIPsFromBFs = numeric(0)
#   PIPsFromBFs[ ] = NA
#   PIPFiles = NULL
#   nLoci = length(BFFileList)
#   for (i in 1:nLoci) {
#     BFFile = BFFileList[i]
#     outputFilePrefix = BFFile
#     locusIndice = which(lociID == i)
#     nSNPs = length(locusIndice)
#     priorOfEachSNP = PIPsFromAnnotations[locusIndice]
#     priorFile = paste0(BFFile, "_iter", iter, ".SNPPrior")
#     write.table(priorOfEachSNP, priorFile, 
#                 row.names = F, col.names = F, quote = F)
#     modelSearchWrapper(BFFile, outputFilePrefix, nSNPs, priorFile)
#     PIPFiles = c(PIPFiles, paste0(outputFilePrefix, ".marginal"))
#   }
#   # collect PIPs and likelihood (proportional to a constant)
#   likelihood = numeric(nLoci)
#   for (i in 1:nLoci) {
#     PIPFile = PIPFiles[i]
#     locusIndice = which(lociID == i)
#     SNPIDAndPIP = read.table(PIPFile)
#     SNPID = SNPIDAndPIP[, 1]
#     PIPsFromBFs[locusIndice[SNPID]] = SNPIDAndPIP[, 2]
#     likelihoodFile = paste0(BFFileList[i], ".statistics")
#     likelihood[i] = as.matrix(read.table(likelihoodFile))[1]
#   }
#   return(list(PIPsFromBFs = PIPsFromBFs, likelihood = likelihood))
# }

calculatePIPsFromBFs = function(PIPsFromAnnotations, lociID, BFFileList, 
                                iter) {
  # browser()
  PIPsFromBFs = numeric(0)
  PIPsFromBFs[ ] = NA
  PIPFiles = NULL
  nLoci = length(BFFileList)
  uniqueLociID = unique(lociID)
  for (i in 1:nLoci) {
    BFFile = BFFileList[i]
    outputFilePrefix = BFFile
    locusIndice = which(lociID == uniqueLociID[i])
    nSNPs = length(locusIndice)
    priorOfEachSNP = PIPsFromAnnotations[locusIndice]
    priorFile = paste0(BFFile, "_iter", iter, ".SNPPrior")
    write.table(priorOfEachSNP, priorFile, 
                row.names = F, col.names = F, quote = F)
    modelSearchWrapper(BFFile, outputFilePrefix, nSNPs, priorFile)
    PIPFiles = c(PIPFiles, paste0(outputFilePrefix, ".marginal"))
  }
  # collect PIPs and likelihood (proportional to a constant)
  likelihood = numeric(nLoci)
  for (i in 1:nLoci) {
    PIPFile = PIPFiles[i]
    locusIndice = which(lociID == uniqueLociID[i])
    SNPIDAndPIP = read.table(PIPFile)
    SNPID = SNPIDAndPIP[, 1]
    PIPsFromBFs[locusIndice[SNPID]] = SNPIDAndPIP[, 2]
    likelihoodFile = paste0(BFFileList[i], ".statistics")
    likelihood[i] = as.matrix(read.table(likelihoodFile))[1]
  }
  return(list(PIPsFromBFs = PIPsFromBFs, likelihood = likelihood))
}

calculatePIPsFromBFModels = function(PIPsFromAnnotations, lociID, BFModelList, 
                                iter = NULL) {
  # input: 
  # iter: is not used
  
  if (!useParallel) {
    PIPsFromBFs = numeric(length(lociID))
    PIPsFromBFs[ ] = NA
    
    uniqueLociID = unique(lociID)
    nLoci = length(uniqueLociID)
    
    likelihood = numeric(nLoci)
    likelihood[ ] = NA
    for (i in 1:nLoci) {
      locusIndice = which(lociID == uniqueLociID[i])
      priorOfEachSNP = PIPsFromAnnotations[locusIndice]
      BFModelList[[i]]$calculateLikelihoodAndPIPs(priorOfEachSNP)
      SNPIDs = BFModelList[[i]]$SNPIDs
      PIPsFromBFs[locusIndice[SNPIDs]] = BFModelList[[i]]$PIPs
      likelihood[i] = BFModelList[[i]]$logLikelihoodUpToAConstant
    }
    
    return(list(PIPsFromBFs = PIPsFromBFs, likelihood = likelihood))
  } else {
    return(calculatePIPsFromBFModelsParallel(PIPsFromAnnotations, lociID, 
                                                 BFModelList))
  }
}

calculatePIPsFromBFModelsParallel = function(PIPsFromAnnotations, lociID, 
                                             BFModelList, iter = NULL) {
  # input: 
  # iter: is not used
  
  library(parallel)
    
  PIPsFromBFs = numeric(length(lociID))
  PIPsFromBFs[ ] = NA
  
  uniqueLociID = unique(lociID)
  nLoci = length(uniqueLociID)
  
  likelihood = numeric(nLoci)
  likelihood[ ] = NA
  
  priorOfEachSNPList = list()
  locusIndiceList = list()
  for (i in 1:nLoci) {
    locusIndice = which(lociID == uniqueLociID[i])
    locusIndiceList[[i]] = locusIndice
    priorOfEachSNPList[[i]] = PIPsFromAnnotations[locusIndice]
  }
  
  # begin parallel calculation for PIPs
  result = mclapply(X = 1:nLoci, FUN = function(i, priorOfEachSNPList, 
                                              BFModelList) {
        BFModelList[[i]]$calculateLikelihoodAndPIPs(priorOfEachSNPList[[i]])
        SNPIDs = BFModelList[[i]]$SNPIDs
        PIPs = BFModelList[[i]]$PIPs
        logLikelihood = BFModelList[[i]]$logLikelihoodUpToAConstant
        return(list(SNPIDs = SNPIDs, PIPs = PIPs, 
                    logLikelihood = logLikelihood))
    }, mc.cores = ncores, priorOfEachSNPList, BFModelList)
  # end parallel
  
  for (i in 1:nLoci) {
    SNPIDs = result[[i]]$SNPIDs
    locusIndice = locusIndiceList[[i]]
    # check before use
    if (sum(is.na(result[[i]]$PIPs)) > 0 || 
        length(result[[i]]$logLikelihood) == 0) {
      print("error locus:", i)
      print("SNPIDs")
      print(SNPIDs)
      print("current prior:")
      print(priorOfEachSNPList[[i]])
      priorFile = paste0("locus", i, ".error.prior")
      write.table(priorOfEachSNPList[[i]], priorFile, row.names = F, 
                  col.names = F, quote = F)
      print("PIPs:")
      print(result[[i]]$PIPs)
      print("logLikelihood") 
      print(result[[i]]$logLikelihood)
      quit(save = "yes")
    }
    PIPsFromBFs[locusIndice[SNPIDs]] = result[[i]]$PIPs
    likelihood[i] = result[[i]]$logLikelihood
  }
  
  return(list(PIPsFromBFs = PIPsFromBFs, likelihood = likelihood))
}

crossEntropyError = function(y, yHat) {
  # the cross entropy betwewn probabilities y and its fitted value yHat
  
  # input
  # y: the observed probabilities
  # yHat: the fitted probabilities
  
  # output
  # the cross entropy
  
  return(-sum(y * log(yHat) + (1 - y) * log(1 - yHat)))
}

crossEntropy = function(gamma, y, x, lambda = 1) {
  # calculate the cross entropy 
  
  # input
  # x: a matrix of predictors, n * p, including the intercept 1
  # gamma: a vector of coefficients, length is p
  
  linearSum = x %*% gamma
  prob = 1 / (1 + exp(-(linearSum)))
  penalty = 0
  n = length(y)
  gamma = as.vector(gamma)
  if (length(gamma) >= 2) {
    penalty = 0.5 * lambda * sum(gamma[2:length(gamma)]^2)
  }
  crossEntropyValue = -sum(y * log(prob) + (1 - y) * log(1 - prob)) + 
    penalty
  return(crossEntropyValue)
}

crossEntropyGradient = function(gamma, y, x, lambda = 1) {
  # gradient of the function crossEntropy
  
  # input
  # see crossEntropy
  muHat = 1 / (1 + exp(-(x %*% gamma)))
  gamma[1] = 0 # the intercept is not penalized
  gradient = -t(x) %*% (y - c(muHat)) + lambda * gamma
  return(c(gradient))
}

crossEntropyHessian = function(gamma, y, x, lambda = 1) {
  # gradient of the function crossEntropy
  
  # input
  # see crossEntropy
  muHat = 1 / (1 + exp(-(x %*% gamma)))
  W = c(muHat * (1 - muHat))
  lambdaDiag = diag(lambda, nrow = dim(x)[2])
  lambdaDiag[1, 1] = 0
  hessian = t(x) %*% (x * W) + lambdaDiag
  return(hessian)
}

crossEntropyNegative = function(gamma, y, x, penalty = 1) {
  return(-crossEntropy(gamma, y, x, penalty))
}
crossEntropyGradientNegative = function(gamma, y, x, penalty = 1) {
  return(-crossEntropyGradient(gamma, y, x, penalty))
}
crossEntropyHessianNegative = function(gamma, y, x, penalty = 1) {
  return(-crossEntropyHessian(gamma, y, x, penalty))
}


optimizeCrossEntropyCV = function(y, x, lambda, fold = 10, bestModel = "min") {
  # use cross validation to decide lambda. The error measure is the 
  # cross entropy measure
  
  # input
  # y: the observed probabilites
  # x: the annotation matrix
  # lambda: the penalites to evaluate, ascendingly sorted
  
  # output
  # the best lambda parameter and the fitted yHat
  
  n = length(y)
  if (fold > n) {
    fold = n
    print("use Leave-one-out cv instead because n < 10")
  }
  stopifnot(fold > 1)
  foldID = sample(rep(1:fold, length = n))
  nLambda = length(lambda)
  cvError = numeric(nLambda)
  cvErrorMatrix = matrix(NA, fold, nLambda)
  for (j in 1:nLambda) {
    lambdaOfThisFold = lambda[j]
    for (i in 1:fold) {
      indiceExcluded = which(foldID == i)
      yTrain = y[-indiceExcluded]
      yTest = y[indiceExcluded]
      xTrain = x[-indiceExcluded, ]
      xTest = x[indiceExcluded, ]
      gammaFitted = optimizeCrossEntropy(yTrain, xTrain, 
                                         lambdaOfThisFold)$gamma
      yPredicted = 1 / (1 + exp(-(xTest %*% gammaFitted)))
      error = crossEntropyError(yTest, yPredicted)
      cvError[j] = cvError[j] + error
      cvErrorMatrix[i, j] = error
    }
  }
  browser()
  IndexOfMinError = which.min(cvError)
  lambdaMinCV = lambda[IndexOfMinError]
  print(paste("lambda with minimum cv error is lambda[", 
              IndexOfMinError, "]", lambdaMinCV))
  if (bestModel == "min") {
    bestIndex = IndexOfMinError
  } else if (bestModel == "1se") {
    stdErrorCV = apply(X = cvErrorMatrix, MARGIN = 2, FUN = sd) / sqrt(fold)
    cvErrorMin1se = cvError[IndexOfMinError] + stdErrorCV[IndexOfMinError]
    IndexOfWithin1se = IndexOfMinError + 1
    while(cvError[IndexOfWithin1se] <= cvErrorMin1se) {
      IndexOfWithin1se = IndexOfWithin1se + 1
    }
    bestIndex = IndexOfWithin1se - 1
  }
  bestLambda = lambda[bestIndex]
  print(paste("best lambda from model selection is lambda[", 
              bestIndex, "]", bestLambda))
  fittedOnFullData = optimizeCrossEntropy(y, x, bestLambda)
  gammaFitted = fittedOnFullData$gamma
  yFitted = fittedOnFullData$yFitted
  
  return(list(yFitted = yFitted, gamma = gammaFitted, 
              cvError = cvError, bestLambda = bestLambda))
}

optimizeCrossEntropy = function(y, x, lambda = 1) {
  # find the minimal cross entropy
  
  # input
  # see crossEntropy
  # browser()
  initGamma = c(log(mean(y) / (1 - mean(y))), rep(0, dim(x)[2] - 1))
  #   optimResult = optim(initGamma, crossEntropy, gr = crossEntropyGradient,
  #                       y, x, lambda,
  #                       control = list(maxit = 1e3), 
  #                       method = "BFGS")
  #   gammaUpdated = optimResult$par
  #   if (optimResult$convergence != 0) {
  #     print("fail to converge")
  #     stop()
  #   }
  
  optimResult2 = maxNR(fn = crossEntropyNegative, 
                       grad = crossEntropyGradientNegative, 
                       hess = NULL, #crossEntropyHessianNegative, 
                       start = initGamma, 
                       lambdatol = 1e-6,  qrtol = 1e-10,
                       y = y, x = x, penalty = lambda)
  # print(paste("max(Q) = ", optimResult2$maximum))
  gammaUpdated = optimResult2$estimate
  if (as.integer(optimResult2$code) != 1L && 
        as.integer(optimResult2$code) != 2L) {
    print(paste("no normal converge with code:", optimResult2$code))
    stop()
  }
  
  # <debug> compare the difference between glm and maxNR
#   model = glm(y ~ x - 1, family = "binomial")
#   gamma = as.matrix(summary(model)$coefficients[, 1])
#   gammaUpdated = gamma
#   fittedPrior = 1 / (1 + exp(-(x %*% gamma)))
#   print(summary(model))
#   W = as.vector(fittedPrior * (1 - fittedPrior))
#   InfoMatrix = t(x) %*% (x * W)
#   print("svd of the information matrix:")
#   print(svd(InfoMatrix)$d)
#   print("rank of the information matrix:")
#   print(qr(InfoMatrix)$rank)
  # </debug>
  
  #print(gammaUpdated)
  #print(paste("current negative cross entropy is:", -optimResult$value))
  linearSum = x %*% gammaUpdated
  yFitted = 1 / (1 + exp(-(linearSum)))
  return(list(yFitted = yFitted, gamma = gammaUpdated))
}

penaltyOfMethods = function(gamma, fittingMethod, methodArgs) {
  if (fittingMethod == "PCEMin" || fittingMethod == "glmnetRidgeMin") {
    penalty = 0.5 * methodArgs[1] * sum(gamma[2:length(gamma)]^2)
  } else if (fittingMethod == "glmnetLASSOMin") {
    penalty = methodArgs[1] * sum(abs(gamma[2:length(gamma)]))
  } else if (fittingMethod == "glmnetENETMin") {
    lambda = methodArgs[1]
    alpha = methodArgs[2]
    penalty = lambda * (0.5 * (1 - alpha) * sum(gamma[2:length(gamma)]^2) + 
              alpha * sum(abs(gamma[2:length(gamma)])))
  } else {
    penalty = 0
    message("no penalty")
  } 
  return(penalty)
}

PIPAndLikelihoodOfGamma = function(gamma, annotation, lociID, BFFileList,
                                   fittingMethod, methodArgs) {
  fittedPrior = 1 / (1 + exp(-(annotation %*% gamma)))
  PIPAndLikelihood = calculatePIPsFromBFs(fittedPrior, lociID, BFFileList,
                                          "numericDev")
  likelihood = PIPAndLikelihood$likelihood
  totalLik = sum(likelihood)
  penaltyTerm = penaltyOfMethods(gamma, fittingMethod, methodArgs)
  totalLik = totalLik - penaltyTerm
  PIP = PIPAndLikelihood$PIPsFromBFs
  return(c(totalLik, PIP))
}

likelihoodOfGamma = function(gamma, annotation, lociID, BFFileList,
                                   fittingMethod, methodArgs) {
  PIPAndLikelihood = PIPAndLikelihoodOfGamma(gamma, 
                                     annotation, lociID, BFFileList,
                                     fittingMethod, methodArgs)
  return(PIPAndLikelihood[[1]])
}

writePIPs = function(PIPs, outputFile) {
  SNPIDPIPs = cbind(1:length(PIPs), PIPs)
  dim(SNPIDPIPs) = c(length(PIPs), 2)
  permutation = order(PIPs, decreasing = T)
  SNPIDPIPs = SNPIDPIPs[permutation, , drop = F]
  # round PIPs 
  SNPIDPIPs[, 2] = signif(SNPIDPIPs[, 2], digits = 7)
  write.table(SNPIDPIPs, outputFile, row.names = F, col.names = F, 
              quote = F)
}

outputResult = function(outputPrefix, gamma, PIPsFromBFsAPP, PIPsFromBFsMAP, 
                        logLikAPP, logLikMAP, lociListFile) {
  gammaFile = paste0(outputPrefix, ".gamma")
  write.table(gamma, gammaFile, row.names = F, col.names = F, 
              quote = F)
#   PIPFile = paste0(outputPrefix, ".app.marginal")
#   writePIPs(PIPsFromBFsAPP, PIPFile)
#   PIPZFile = paste0(PIPFile, "z")
#   connectPIPFileToZFile(lociListFile, PIPFile, PIPZFile)
  PIPFile = paste0(outputPrefix, ".marginal")
  writePIPs(PIPsFromBFsMAP, PIPFile)
  PIPZFile = paste0(PIPFile, "z")
  connectPIPFileToZFile(lociListFile, PIPFile, PIPZFile)
  # likelihoodFile = paste0(outputPrefix, ".app.loglik")
  # write.table(logLikAPP, likelihoodFile, row.names = F, col.names = F, 
  #             quote = F)
  likelihoodFile = paste0(outputPrefix, ".loglik")
  write.table(logLikMAP, likelihoodFile, row.names = F, col.names = F, 
              quote = F)
}

generateLociID = function(lociListFile) {
  lociList = c(as.matrix(read.table(lociListFile, colClasses = "character")))
  lociID = NULL
  for (i in 1:length(lociList)) {
    lines = readLines(paste0(dirname(lociListFile), "/", lociList[i]))
    nLines = length(lines)
    lociID = c(lociID, rep(i, nLines))
  }
  return(lociID)
}

fittedUsingGlmnetENET = function(y, x, lambda, alpha = 0.5, exact = F) {
  if (dim(x)[2] == 1) { # fit only 1s
    meanY = mean(y)
    yFitted = rep(meanY, length(y))
    gamma = log(meanY / (1 - meanY))
    return(list(yFitted = yFitted, gamma = gamma))
  }
  nObservation = dim(x)[1]
  lambda = lambda / nObservation
  lambdaSeries = lambda # * 2^seq(from = 5, by = -0.5, to = -2)
  weights = rep(nObservation, nObservation)
  combinedY = cbind(1-y, y)
  model = glmnet::glmnet(x = x, y = combinedY, #lambda = lambdaSeries, 
                 #weights = weights,
                 family = "binomial", alpha = alpha, standardize = F)
  # <debug>
  # exact = F
  # </debug>
  if (exact) {
    assign("x", x, envir = .GlobalEnv)
    assign("combinedY", combinedY, envir = .GlobalEnv)
    assign("alpha", alpha, envir = .GlobalEnv)
  }
  yFitted = predict(model, newx = x, type = "response", s = lambda, 
                    exact = exact)
  gammaFitted = as.matrix(predict(model, x, 
                                  type = "coefficients", s = lambda, 
                                  exact = exact))
  gamma = c(gammaFitted[1], gammaFitted[3:(dim(x)[2] + 1)])
  return(list(yFitted = yFitted, gamma = gamma))
}

calculatePriorFromAnnotations = function(annotation, PIPsFromBFs, 
                                         fittingMethod = "PCEMin",
                                         parameters = NULL, exact = F) {
  # predict PIPs based on annotation. 
  # Several options:
  # 1. use L1 logistic regression, grouped data with the same number of 
  #    observations in the group
  # 2. use logistic regression
  # The current version use logistic regression, even though there is a 
  # warning for non-integer successes
  
  # browser()
  gamma = NULL
  # t1 = system.time({
  if (fittingMethod == "glmnetLASSOMin" || 
        fittingMethod == "glmnetLASSO1se") {
    fittedResult = fittedUsingGlmnetENET(PIPsFromBFs, annotation, 
                                         lambda = parameters[1], alpha = 1,
                                         exact)
    fittedPrior = fittedResult$yFitted
    gamma = fittedResult$gamma
  } else if (fittingMethod == "glmnetENETMin" || 
               fittingMethod == "glmnetENET1se") {
    fittedResult = fittedUsingGlmnetENET(PIPsFromBFs, annotation, 
                                         lambda = parameters[1], 
                                         alpha = parameters[2],
                                         exact)
    fittedPrior = fittedResult$yFitted
    gamma = fittedResult$gamma
  } else if (fittingMethod == "glmnetRidgeMin" || 
               fittingMethod == "glmnetRidge1se") {
    fittedResult = fittedUsingGlmnetENET(PIPsFromBFs, annotation, 
                                         lambda = parameters[1], alpha = 0,
                                         exact)
    fittedPrior = fittedResult$yFitted
    gamma = fittedResult$gamma
  } else if (fittingMethod == "PCEMin" ||
               fittingMethod == "PCE1se"
  ) {
    CEResult = optimizeCrossEntropy(PIPsFromBFs, annotation, 
                                    lambda = parameters[1])
    fittedPrior = CEResult$yFitted
    gamma = CEResult$gamma
  }  
  # })
  # print("t1:")
  # print(t1)
  # print(summary(glm(PIPsFromBFs ~ annotation - 1, family = "binomial")))
  # coef = summary(glm(PIPsFromBFs ~ annotation - 1, 
  #                    family = "binomial"))$coefficients[, 1]
  # print(crossEntropy(coef, PIPsFromBFs, annotation), digits = 16)
  # browser()
  # t2 = system.time({
  # priorFromAnnotations2 = optimizeCrossEntropy(PIPsFromBFs, annotation)
  # })
  # print("t2:")
  # print(t2)
  return(list(fittedPrior = fittedPrior, gamma = gamma))
}

EMIteration = function(annotation, lociID, BFModelList, 
                       gammaInit, fittingMethod,
                       maxIter, deltaLik, workingDir, 
                       methodArgs, PIPFilePrefix = NULL,
                       exact = F) {
  # Set exact to T improves the fitting precision. This exact parameter is the
  # "exact" parameter for the exact Bayes factor.
  
  #browser()
  iter = 0L
  PIPsFromBFs = NULL
  PIPsFromBFsCurrent = NULL
  likelihoodPre = -Inf
  sumLogLikCurrent = -Inf
  gammaCurrent = gammaInit
  gammaNew = gammaInit
  nSNPAllLoci = length(lociID)
  
  causalPriorProb = 1 / (1 + exp(-(annotation %*% gammaNew)))
  while (iter <= maxIter) {
    # browser()
    # print(paste("== iteration", iter))
    
    timeStart = proc.time()
    ResultFromBFs = calculatePIPsFromBFModels(causalPriorProb, lociID, 
                                              BFModelList, iter)
    time = proc.time() - timeStart
    if (verbose) {
      print("Time used to calculate PIPs:")
      print(time)
    }
    likelihood = ResultFromBFs$likelihood    
    PIPsFromBFs = ResultFromBFs$PIPsFromBFs
    sumLogLik = sum(likelihood)
    
    # TODO: think about how to add penalty here
    #     penaltyTerm = penaltyOfMethods(gammaNew, fittingMethod, methodArgs)
    #     sumLogLik = sumLogLik - penaltyTerm
    
#     if (length(gammaNew) <= 4) {
#       print(paste("gamma:", paste(gammaNew, collapse = " ")))
#     }
#     print(paste("(penalized) log likelihood up to a constant", 
#                 "assuming independence among loci", sumLogLik))
    
    if (sumLogLik < sumLogLikCurrent + deltaLik) {
#       print(paste("Early stop because the increase of the log likelihood < ", 
#                   deltaLik))
      # <debug>
      # browser()
      # </debug>
      break   
    }
    
#     if (iter > 0L) {
#       maxDiff = max(abs(PIPsFromBFs - PIPsFromBFsCurrent))
#       print(paste("max difference in PIP updating =", maxDiff))
#     }
    
    # TODO: think whether this stop is reasonable
    #     if ((iter > 0L) && (maxDiff < 1e-5)) {
    #       print("Stop because miminum differences in PIPs reached")
    #       break
    #     }
    
    PIPsFromBFsCurrent = PIPsFromBFs
    sumLogLikCurrent = sumLogLik
    gammaCurrent = gammaNew
    
    #     # write the result for debug
    #     SNPIDPIPs = cbind(1:nSNPAllLoci, PIPsFromBFs)
    #     permutation = order(PIPsFromBFs, decreasing = T)
    #     SNPIDPIPs = SNPIDPIPs[permutation, ]
    #     if (!is.null(PIPFilePrefix)) {  # controlled output
    #       PIPFile = paste0(PIPFilePrefix, "_", "iter", iter, ".marginal")
    #       write.table(SNPIDPIPs, PIPFile, row.names = F, col.names = F, 
    #                   quote = F)
    #     }
    
    #     if (iter == 0) {
    #       sumLogLikNull = sumLogLik
    #       SNPIDPIPsNull = SNPIDPIPs
    #       gammaNull = gammaNew
    #     }
    time = system.time({
      priorUpdateResult = calculatePriorFromAnnotations(annotation, 
                                                        PIPsFromBFs, fittingMethod, methodArgs,
                                                        exact)
    })
    if (verbose) {
      print("Time used to update parameters:")
      print(time)
    }
    causalPriorProb = priorUpdateResult$fittedPrior
    gammaNew = priorUpdateResult$gamma
    iter = iter + 1
  }
  if (iter > maxIter) {
    print("Stop because max iteration reached")
  }
  ## <debug>
  # print(paste("logLik:", sumLogLikCurrent))
  # print("gamma")
  # print(as.matrix(gammaCurrent))
  ## </debug>
  return(list(PIPsFromBFs = PIPsFromBFsCurrent, gamma = gammaCurrent, 
              logLik = sumLogLikCurrent))
}

effectiveNumberOfParam = function(X, lambda, gamma, fittingMethod, 
                                  alpha = NULL) {
  # calculate the effective number of parameters
  
  # input
  # X: the input matrix
  # lambda: the penalty value for ridge regression and elastic net
  #         Attention: this lambda is different from the lambda in glmnet
  
  # output
  # the effective number of parameter
  
  # effective degree of freedom for logistic regression
  yHat = c(1 / (1 + exp(-(X %*% gamma))))
  w = yHat * (1 - yHat)
  X = X * sqrt(w)
  if (fittingMethod == "PCEMin" || fittingMethod == "glmnetRidgeMin") {
    # center X
    d = svd(X)$d
    # only keep non-zero singular values
    d = d[d > 1e-12]
    effectiveNumber = sum(d^2 / (d^2 + lambda))
  } else if (fittingMethod == "glmnetLASSOMin" || 
               fittingMethod == "glmnetENETMin") {
    # exclude the intercept
    indNonZero = which(abs(gamma) > 1e-12)
    Xactive = X[, indNonZero, drop = F] 
    if (fittingMethod == "glmnetLASSOMin") {
      effectiveNumber = qr(Xactive)$rank
    } else if (fittingMethod == "glmnetENETMin") {
      if (is.null(alpha)) {
        stop("alpha should not be null")
      }
      lambda2 = lambda * (1 - alpha)
      d = svd(Xactive)$d
      # only keep non-zero singular values
      d = d[d > 1e-12]
      effectiveNumber = sum(d^2 / (d^2 + lambda2))
    }
  } 
  
  # # effective degree of freedom for linear regression
  #   if (fittingMethod == "PCEMin" || fittingMethod == "glmnetRidgeMin") {
  #     # center X
  #     mX = colMeans(X)
  #     X = t(t(X) - mX)
  #     d = svd(X)$d
  #     # only keep non-zero singular values
  #     d = d[d > 1e-12]
  #     effectiveNumber = sum(d^2 / (d^2 + lambda)) + 1 # plus the intercept
  #   } else if (fittingMethod == "glmnetLASSOMin" || 
  #              fittingMethod == "glmnetENETMin") {
  #     # exclude the intercept
  #     indNonZero = which(abs(gamma[2:length(gamma)]) > 1e-12)
  #     if (length(indNonZero) == 0L) {
  #       return(1)  # nothing is selected
  #     }
  #     Xactive = X[, 1 + indNonZero, drop = F] 
  #     mX = colMeans(Xactive)
  #     XactiveCentered = t(t(Xactive) - mX)
  #     if (fittingMethod == "glmnetLASSOMin") {
  #       effectiveNumber = qr(XactiveCentered)$rank + 1
  #     } else if (fittingMethod == "glmnetENETMin") {
  #       if (is.null(alpha)) {
  #         stop("alpha should not be null")
  #       }
  #       lambda2 = lambda * (1 - alpha)
  #       d = svd(XactiveCentered)$d
  #       # only keep non-zero singular values
  #       d = d[d > 1e-12]
  #       effectiveNumber = sum(d^2 / (d^2 + lambda2)) + 1 # plus the intercept
  #     }
  #   } 
  
  return(effectiveNumber)
}

AICBICOfEachHyperparamete = function(lambda, alpha,
                                     annotation, lociID, BFModelList,
                                     gammaInit, fittingMethod, 
                                     maxIter, deltaLik, workingDir, 
                                     nSample) {
  # calculate AIC or BIC 
  # calculation of AIC and BIC on the fine mapping model  
  # bic2 use nSNP, bic uses nSample in BIC.   
  
  modelFromFullData = EMIteration(annotation, lociID, BFModelList, 
                                  gammaInit, fittingMethod,
                                  maxIter, deltaLik, workingDir, 
                                  c(lambda, alpha), NULL)
  gamma = modelFromFullData$gamma
  logLik = modelFromFullData$logLik
  PIPsFromBFs = modelFromFullData$PIPsFromBFs
  nEffective = effectiveNumberOfParam(annotation, lambda, 
                                      gamma, fittingMethod, alpha = alpha)
  nSNP = dim(annotation)[1]
  aic = logLik - nEffective
  bic = logLik - log(nSample) / 2 * nEffective  
  bic2 = logLik - log(nSNP) / 2 * nEffective
  return(list(aic = aic, bic = bic, bic2 = bic2, 
              PIPsFromBFs = PIPsFromBFs, gamma = gamma,
              A = NULL, logLik = logLik))
}

evidenceOfEachHyperparameter = function(lambda, alpha,
                                annotation, lociID,  
                                BFModelList,
                                gammaInit, fittingMethod, 
                                maxIter, deltaLik, trueHessian, workingDir) {
  # make sure it works for ENET
  
  nGamma = dim(annotation)[2]
  modelFromFullData = EMIteration(annotation, lociID, BFModelList, 
                                  gammaInit, fittingMethod,
                                  maxIter, deltaLik, workingDir, 
                                  c(lambda, alpha), NULL)
  gamma = modelFromFullData$gamma
  # predict
  # browser()
  if (trueHessian) {
    #runTime = system.time({
    jacobianMatrix = jacobian(PIPAndLikelihoodOfGamma, gamma, 
                              annotation = annotation, 
                              lociID = lociID, BFModelList = BFModelList, 
                              fittingMethod = fittingMethod, 
                              methodArgs = c(lambda, 0.5))
    #})
    #print(runTime)
    print("derivative of the likelihood of gammas")
    print(jacobianMatrix[1, ])
    J = jacobianMatrix[2:dim(jacobianMatrix)[1], ]
  } else {
    J = NULL
  }
  print(paste("gamma:", gamma))
  A = calculateHessianOfNegLik(annotation, gamma, J, lambda, 
                               fittingMethod, alpha = alpha)
  eigenvalues = eigen(A)$values
  logdetA = sum(log(eigenvalues)) 
  
  causalPriorProb = 1 / (1 + exp(-(annotation %*% gamma)))
  ResultFromBFs = calculatePIPsFromBFModels(causalPriorProb, lociID, 
                                            BFModelList, paste0("full"))
  likelihood = ResultFromBFs$likelihood    
  PIPsFromBFs = ResultFromBFs$PIPsFromBFs
  sumLogLik = sum(likelihood)
  print(paste("logLik:", sumLogLik))
  penalizedLogLik = sumLogLik - lambda / 2 * sum(gamma[2:nGamma]^2)
  print(paste("penalized logLik:", penalizedLogLik))
  print(paste("log(det(A))", logdetA))
  print(paste("0.5 * (nGamma - 1) * log(lambda):", 
              1 / 2 * (nGamma - 1) * log(lambda)))
  evidence = penalizedLogLik - 1 / 2 * logdetA + 
    1 / 2 * (nGamma - 1) * log(lambda)
  print(paste("evidence:", evidence))
  return(list(evidence = evidence, PIPsFromBFs = PIPsFromBFs, gamma = gamma,
              A = A, logLik = sumLogLik))
}

CVLikelihoodOfEachHyperparameter = function(lambda, alpha,
                                annotation, lociID,
                                BFModelList,
                                gammaInit, fittingMethod, 
                                maxIter, deltaLik, trueHessian, posteriorUsage,
                                workingDir, foldID) {
  nFold = length(unique(foldID))
  logLike = numeric(nFold)
  logLike[ ] = NA
  for (i in 1:nFold) {
    print(paste("# fold", i))
    lociTrain = which(foldID != i)
    lociTest = which(foldID == i)
    SNPIndiceTrain = which(is.element(lociID, lociTrain))
    SNPIndiceTest = which(is.element(lociID, lociTest))
    lociIDTrain = lociID[SNPIndiceTrain]
    lociIDTest =  lociID[SNPIndiceTest]
    annotationTrain = annotation[SNPIndiceTrain, , drop = F]
    annotationTest = annotation[SNPIndiceTest, , drop = F]
    BFModelListTrain = BFModelList[lociTrain]
    BFModelListTest = BFModelList[lociTest]
    
    # train
    model = EMIteration(annotationTrain, lociIDTrain, BFModelListTrain, 
                        gammaInit, fittingMethod,
                        maxIter, deltaLik, workingDir, 
                        c(lambda, alpha), NULL)
    gamma = model$gamma
    ## <debgu>
    # browser()
    ## </debug>
    
    # test
    if (posteriorUsage == "MAP") {
      priorProb = 1 / (1 + exp(-(annotationTest %*% gamma)))
      # browser()
    } else if (posteriorUsage == "approximate") {
      # compare the Hessian matrix from numerical values and partially from
      # analytical versions
      if (trueHessian) {
        runTime = system.time({
          jacobianMatrix = jacobian(PIPAndLikelihoodOfGamma, gamma, 
                                    annotation = annotationTrain, 
                                    lociID = lociIDTrain, BFModelList = BFModelListTrain, 
                                    fittingMethod = fittingMethod, 
                                    methodArgs = lambda)
        })
        print(runTime)
        J = jacobianMatrix[2:dim(jacobianMatrix)[1], ]
      } else {
        J = NULL
      }
      
      #         runTime = system.time({
      #           hessianMatrix = hessian(likelihoodOfGamma, gamma, 
      #                          annotation = annotationTrain, 
      #                          lociID = lociIDTrain, BFModelList = BFModelListTrain, 
      #                          fittingMethod = fittingMethod, 
      #                          methodArgs = lambda[j])
      #         })
      #         print(runTime)
      #         print("numerical hessian")
      #         print(hessianMatrix)
      
      A = calculateHessianOfNegLik(annotationTrain, gamma, J, lambda,
                                   fittingMethod = fittingMethod, alpha)
      priorProb = predictivePrior(annotationTrain, gamma, A, 
                                  annotationTest)
    }
    
    ResultFromBFs = calculatePIPsFromBFModels(priorProb, lociIDTest, 
                                              BFModelListTest, paste0("fold", i))
    likelihood = ResultFromBFs$likelihood
    logLike[i] = sum(likelihood)
    print(paste("logLik of fold", i, "is", logLike[i]))
  }
  
  cvlogLik = sum(logLike)
  return(list(cvlogLik = cvlogLik, foldLogLik = logLike))
}

scoreBasedHyperparameterSelection = function(annotation, lociID, 
                                 BFModelList, nSample, PIPFilePrefix,
                                 gammaInit, fittingMethod, lambda, alpha,
                                 maxIter, deltaLik, trueHessian, posteriorUsage,
                                 workingDir, scoreCriterion) {
  # use evidence approximation, AIC, or BIC to select hyperparameters
  
  # input
  # scoreCriterion: "ev", "aic" or "bic"
  
  # output
  
  nLambda = length(lambda)
  nAlpha = length(alpha)
  maxModelScore = -Inf
  PIPsFromBFs = NULL
  PIPsFromBFsMAP = NULL
  gamma = NULL
  A = NULL
  logLik = NULL
  logLikMAP = NULL
  for (i in 1:nLambda) {
    for (j in 1:nAlpha) {
      print(paste("## lambda =", lambda[i], "alpha =", alpha[j]))
      if (scoreCriterion == "ev") {
        result = evidenceOfEachHyperparameter(lambda[i], alpha[j],
                                              annotation, lociID, 
                                              BFModelList,
                                              gammaInit, fittingMethod, 
                                              maxIter, deltaLik, trueHessian, workingDir)
        modelScore = result$evidence
      } else if (scoreCriterion == "aic" || scoreCriterion == "bic" || 
                   scoreCriterion == "bic2") {
        result = AICBICOfEachHyperparamete(lambda[i], alpha[j],
                                           annotation, lociID, BFModelList,
                                           gammaInit, fittingMethod, 
                                           maxIter, deltaLik, workingDir, 
                                           nSample)
        if (scoreCriterion == "aic") {
          modelScore = result$aic
        } else if (scoreCriterion == "bic") {
          modelScore = result$bic
        } else if (scoreCriterion == "bic2") {
          modelScore = result$bic2
        }
      }
      print(paste("modelScore:", modelScore))
      if (modelScore > maxModelScore) {
        bestLambda = lambda[i]
        bestAlpha = alpha[j]
        PIPsFromBFsMAP = result$PIPsFromBFs
        gamma = result$gamma
        A = result$A
        logLikMAP = result$logLik
        maxModelScore = modelScore
      }
      print(paste("current best: lambda =", bestLambda, "alpha =", bestAlpha))
    }
  }
  
  # rerun again with the best parameter with the exact predict from glmnet
  # and store intermediate iterations
  modelFromFullData = EMIteration(annotation, lociID, BFModelList, 
                                  gammaInit, fittingMethod,
                                  maxIter, deltaLik, workingDir, 
                                  c(bestLambda, bestAlpha), PIPFilePrefix,
                                  exact = T)
  gamma = modelFromFullData$gamma
  PIPsFromBFsMAP = modelFromFullData$PIPsFromBFs
  logLikMAP = modelFromFullData$logLik
  
  if (scoreCriterion == "aic" || scoreCriterion == "bic" || 
        scoreCriterion == "bic2") {
    if (trueHessian) {
      #runTime = system.time({
      jacobianMatrix = jacobian(PIPAndLikelihoodOfGamma, gamma, 
                                annotation = annotation, 
                                lociID = lociID, BFModelList = BFModelList, 
                                fittingMethod = fittingMethod, 
                                methodArgs = c(bestLambda, bestAlpha))
      #})
      #print(runTime)
      print("derivative of the likelihood of gammas")
      print(jacobianMatrix[1, ])
      J = jacobianMatrix[2:dim(jacobianMatrix)[1], ]
    } else {
      J = NULL
    }
    print(paste("gamma:", gamma))
    A = calculateHessianOfNegLik(annotation, gamma, J, bestLambda,
                                 fittingMethod = fittingMethod, bestAlpha)
  }
  
  if (posteriorUsage == "approximate") {
    priorProb = predictivePrior(annotation, gamma, A, annotation) 
    ResultFromBFs = calculatePIPsFromBFModels(priorProb, lociID, 
                                              BFModelList, paste0("full"))
    likelihood = ResultFromBFs$likelihood    
    PIPsFromBFs = ResultFromBFs$PIPsFromBFs
    logLik = sum(likelihood)
  } else if (posteriorUsage == "sampling") {
    sampledResult = sampleBasedPIPs(mean = gamma, variance = solve(A), 
                                    annotation, lociID, BFModelList,
                                    nSamples = 1e3)
    PIPsFromBFs = sampledResult$PIPs
    logLik = sampledResult$likelihood
  } 
  
  return(list(PIPsFromBFs = PIPsFromBFs, gamma = gamma, logLik = logLik,
              PIPsFromBFsMAP = PIPsFromBFsMAP, logLikMAP = logLikMAP))
}

sampleBasedPIPs = function(mean, variance, annotation, lociID, BFModelList,
                           nSamples = 1000) {
  gammaSamples = mvrnorm(n = nSamples, mu = mean, Sigma = variance)
  likelihoodAll = numeric(nSamples)
  PIPsAll = matrix(NA, nSamples, dim(annotation)[1])
  for (i in 1:nSamples) {
    if (i %% 20 == 0) {
      print(i)
    }
    gamma = gammaSamples[i, ]
    causalPriorProb = 1 / (1 + exp(-(annotation %*% gamma)))
    #causalPriorProb[causalPriorProb < 1e-2] = 1e-2 # avoid 0 prior prob
    ResultFromBFs = calculatePIPsFromBFModels(causalPriorProb, lociID, 
                                              BFModelList, "sampling")
    logLikelihood = sum(ResultFromBFs$likelihood)
    PIPs = ResultFromBFs$PIPsFromBFs
    likelihoodAll[i] = logLikelihood
    PIPsAll[i, ] = PIPs
  }
  likelihoodAveraged = mean(likelihoodAll, na.rm = T)
  PIPsAveraged = colMeans(PIPsAll, na.rm = T)
  nNotNA = sum(!is.na(likelihoodAll))
  stdLikelihoodAveraged = sd(likelihoodAll, na.rm = T) / sqrt(nNotNA - 1)
  stdPIPs = apply(PIPsAll, 2, sd, na.rm = T) / 
    sqrt(colSums(!is.na(PIPsAll)) - 1)
  return(list(likelihood = likelihoodAveraged, PIPs = PIPsAveraged,
              stdLik = stdLikelihoodAveraged,
              stdPIPs = stdPIPs))
}
