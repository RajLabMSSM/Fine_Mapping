#' Fine mapping causal variants in multiple loci
#' 
#' @param lociListFile A file listing all the z score files. 
#' @param maxCausal The assumed maximal number of causal variants in each locus
#' @param nSample The number of individuals used to calculate summary statistics
#' @param priorType01 0 to specify sigmaa. An experimental 1 is used to specify
#'        the proportion of the phenotype variance explained (pve)
#' @param priorValue If priorType01 is 0, it specifies the value of sigmaa. 
#'        Otherwise, it is the pve value. For sigmaa, a vector of values can 
#'        be used and the final Bayes factors are averaged. For GWAS fine
#'        mapping, c(0.1, 0.2, 0.4) might be a good option. For eQTL, 
#'        c(0.1, 0.2, 0.4, 0.8, 1.6) might be a good option to accout for large
#'        eQTL effects. 
#' @param exact Whether to calculate the exact Bayes factors. This is only valid
#'        for quantitative traits. For binary traits, only approximated Bayes
#'        factors are available
#' @param eps The small value to add to the correlation matrix. This is useful
#'        when the LD matrix is from a reference panel
#' @param useIdentityMatrix If true, use the identity matrix as the LD matrix. 
#'        This improves the computing speed when only considering one causal 
#'        variant for each locus
#' @param BFFileListFile Specify the files of precomputed Bayes factors 
#'        corresponding to each locus in the lociListFile. 
#'        If NA, the Bayes factors will be 
#'        calculated first. Bayes factors can be precomputed for each locus
#'        using the C++ program caviarbf. This has the advantage of running 
#'        in parallel, and can be saved for later fine mapping.
#' @param priorProb The intial probability of each variant being causal. If it
#'        is 0, the intial probability will be the number of loci divided by 
#'        the number of all variants
#' @param fittingMethod Specify different penalization models: glmnetLASSOMin 
#'        for L1 penalization, glmnetRidgeMin for L2 penalization, 
#'        glmnetENETMin for elastic net penalization. This is only used when
#'        hyperParamSelection is "cv".
#' @param hyperParamSelection Method to select hyperparameter lambda and alpha.
#'        "cv" for loci based cross validation to select
#'        lambda and alpha, "aic", "bic", "bic2" uses AIC, BIC models. For the
#'        BIC model, "bic2" use the number of SNPs as the number of data points, 
#'        "bic" uses the number of individuals as the number of data points, but
#'        they often show similar results. "topK" produces the top 
#'        significant and relative independent annotations. No penalization 
#'        is used in the fitting. As a by-product, it also outputs the
#'        pvalues for each individual annotation when considered separately
#' @param nFold The fold of cross validation, only used when hyperParamSelection
#'        is "cv". If the number of loci is less than then nFold, the number 
#'        of loci is used
#' @param K Specify the top K annotations to select, used for "topK"
#' @param rThreshold Specify the upper bound of the correlation among
#'        top annotations, used for "topK"
#' @param pvalueThreshold Specify the p value threshold for top annotations, 
#'        used for "topK"
#' @param lambda The lambda values for parameter selection. If lambda is 0, 
#'        then no penalization is used in the fitting. This can be used for 
#'        statistical testing, e.g., comparing two nested annotation models.
#' @param alpha The alpha values for parameter selection, only used when 
#'        the fitting method is glmnetENETMin
#' @param annotationIndices A vector to specify annotations. If it is -1, 
#'        use all annotations, 0 for only the intercept, otherwise it 
#'        specifies the corresponding columns of annotations. 
#'        Default uses all annotations.
#' @param LDSuffix The correlation file suffix, default ".LD"
#' @param annotationSuffix The annotation suffix, default ".annotations"
#' @param outputPrefix The prefix of the output files
#' @param keepBF Whether to keep produced Bayes factor files
#' @param overwriteExistingResults Whether to overwrite existing output files
#' @param maxIter The maximal iteration in the EM method
#' @param deltaLik The minimum improvement of log likelihood to continue 
#'        the EM iteration
#' @param useParallel Whether to use library parallel to use multiple cores
#'        This can increase the speed of calculate PIPs in the EM iteration
#' @param ncores Number of cores to use in parallel when useParallel is T
#' @param verbose If true, more information are printed
#' @keywords caviarbf
#' @export
#' @return The results are saved in several files. 
#'         <outputPrefix>.marginal saves the indices of all variants following 
#'         the order in the list file and the their PIPs sorted in a descring 
#'         order. 
#'        <outputPrefix>.marginalz includes the extra region index, 
#'        the variant ID and the z score. 
#'        <outputPrefix>.loglik stores the relative log likelihood.
#'        <outputPrefix>.gamma stores the estimated annotation effect size 
#'        and the effect size of standardized annotation.
#'        <outputPrefix>.log stores the running log.
#'        <outputPrefix>.time saves the running time.   
#'        
#' 
#' @useDynLib caviarbf
#' @importFrom Rcpp evalCpp
caviarbfFineMapping = function(lociListFile,  
                               maxCausal, nSample,
                               priorType01 = 0, 
                               priorValue = c(0.1, 0.2, 0.4), 
                               exact = T, eps = 0,
                               useIdentityMatrix = F,
                               BFFileListFile = NA,
                               priorProb = 0, 
                               fittingMethod = "glmnetLASSOMin",
                               hyperParamSelection = "cv",
                               nFold = 5,
                               K = 10,
                               rThreshold = 0.2,
                               pvalueThreshold = 0.05,
                               lambda = c(2^seq(from = -15, by = 2, to = 5), 
                                          1e2, 1e3, 1e4, 1e5, 1e6),
                               alpha = c(0, 0.2, 0.3, 0.5, 0.7, 0.8, 1),
                               annotationIndices = -1,   
                               annotationSuffix = ".annotations",
                               LDSuffix = ".LD",
                               outputPrefix,
                               keepBF = F,
                               overwriteExistingResults = F,
                               maxIter = 50, 
                               deltaLik = 0.01,
                               useParallel = F, ncores = 1,
                               verbose = F
) {
  startTime = proc.time()
  
  trueHessian = F
  pairTestThreshold = NULL
  
  # add into the global environment
  
  assign("useParallel", useParallel, .GlobalEnv) 
  assign("ncores", ncores, .GlobalEnv)
  assign("verbose", verbose, .GlobalEnv)
  
  posteriorUsage = "MAP" 
  if (hyperParamSelection == "cv" || hyperParamSelection == "fb") {
    # two options: "MAP" or "approximate"
    posteriorUsage = "MAP"  # decide to use MAP or approximate in each CV fold
    #posteriorUsage = "approximate"
  } else if (hyperParamSelection == "topK") {
    fittingMethod = "glmnetRidgeMin"
    posteriorUsage = "ML"
  }
  
  if (is.na(lambda[1])) {
    lambda = c(2^seq(from = -15, by = 2, to = 5), 1e2, 1e3, 1e4, 1e5, 1e6)  
  }
  if (length(lambda) == 1 && lambda[1] == 0) {
    fittingMethod = "glmnetRidgeMin"
  }
  
  if (fittingMethod == "glmnetENETMin") {
    if (is.na(alpha)) {
      alpha = c(0, 0.2, 0.3, 0.5, 0.7, 0.8, 1)
    } 
  } else {
    alpha = NA 
  }
  
  workingDir = dirname(outputPrefix)
  
  if(!file.exists(workingDir)) {
    dir.create(workingDir)
  } else if (!overwriteExistingResults) {
    startNumber = 1
    workingDirOriginal = workingDir
    while (T) {
      workingDir = paste0(workingDirOriginal, "_s", startNumber)
      if (!file.exists(workingDir)) {
        dir.create(workingDir)
        break
      }
      startNumber = startNumber + 1
    }
  }
  
  logFile = paste0(outputPrefix, ".log")
  sink(logFile, split = T)
    
  # for general running
  tempDir = paste0(workingDir, "/temporary_", basename(lociListFile))
  if(!file.exists(tempDir)) {
    dir.create(tempDir, recursive = T)
  }
  
  lociID = generateLociID(lociListFile) 
  nSNPAllLoci = length(lociID)
  nLoci = length(unique(lociID))
  print(paste0("number of loci: ", nLoci, 
               ", number of SNPs in all loci:", nSNPAllLoci))
  
  if (annotationIndices[1] == 0) {
    annotation = matrix(1, nSNPAllLoci, 1)
  } else {
    annotationResult = aggregateAnnotation(lociListFile, annotationSuffix)
    annotation = annotationResult$annotation
  }
  if (annotationIndices[1] != -1 && annotationIndices[1] != 0) {
    annotation = annotation[, annotationIndices, drop = F]
  }
  if (annotationIndices[1] != 0) {
    # add constant 1s into the annotations
    annotation = cbind(1, annotation)
  }
  
  if (abs(priorProb) < 1e-12) {
    causalPriorProb = nLoci / nSNPAllLoci               
  } else {
    causalPriorProb = priorProb
  }
  gammaInit = c(log(causalPriorProb / (1 - causalPriorProb)),
                rep(0, dim(annotation)[2] - 1))
  
  if (is.na(BFFileListFile)) {
    # generate Bayes factors
    BFFileList = CAVIARBFMultiLoci(tempDir, lociListFile, LDSuffix,
          priorType01, priorValue, nSample, maxCausal, exact, eps, 
          useIdentityMatrix)
    PIPFilePrefix = paste0(tempDir, "/", basename(lociListFile), 
                           "_l", maxCausal)
  } else {
    BFFileList = c(as.matrix(read.table(BFFileListFile, as.is = T)))[1:nLoci]
    PIPFilePrefix = paste0(tempDir, "/", basename(BFFileListFile))
  }
  BFModelList = createBFModels(BFFileList, lociID)
  
  if (hyperParamSelection == "cv") {
    result = CVPartitionOnLoci(annotation, lociID,
                               BFModelList, nSample, PIPFilePrefix,
                               gammaInit, fittingMethod, lambda, alpha,
                               maxIter, deltaLik, trueHessian, 
                               posteriorUsage, tempDir, nFold,
                               pairTestThreshold)
  } else if (hyperParamSelection == "ev" || 
               hyperParamSelection == "aic" || 
               hyperParamSelection == "bic" ||
               hyperParamSelection == "bic2") {
    result = scoreBasedHyperparameterSelection(annotation, 
                                 lociID, 
                                 BFModelList, nSample, PIPFilePrefix,
                                 gammaInit, fittingMethod, lambda, alpha,
                                 maxIter, deltaLik, trueHessian, 
                                 posteriorUsage, tempDir, hyperParamSelection)
  } else if (hyperParamSelection == "topK") {
    result = topKBasedonMaximumLikelihood(annotation, lociID,  
                                          BFModelList, nSample, PIPFilePrefix,
                                          gammaInit, fittingMethod, 
                                          maxIter, deltaLik, trueHessian, 
                                          posteriorUsage, tempDir, 
                                          K, rThreshold, pvalueThreshold)
  } else if (hyperParamSelection == "fb") {
    nAnnotation = dim(annotation)[2] - 1
    stopifnot(nAnnotation >= 1)
    data = list()
    data$annotation = annotation
    data$lociID = lociID
    data$lociListFile = lociListFile
    data$BFModelList = BFModelList
    data$nSample = nSample
    data$gammaInit = gammaInit
    data$maxIter = maxIter
    data$deltaLik = deltaLik
    data$trueHessian = trueHessian
    data$tempDir = tempDir
    result = annotationSelectionFgwasPipeline(caviarbfFgwasPipeline, 
                                              data, nAnnotation, 
                                              PIPFilePrefix,
                                              nFold,
                                              startWithTopK = NULL,
                                              otherArguments = NULL,
                                              otherArgumentsFinal = NULL,
                                              significanceModel = "nested")
    gamma = numeric(nAnnotation + 1)
    gamma[1] = result$resultFinalModel$gamma[1]
    if (length(result$candidateSet) > 0) {
      gamma[1 + result$candidateSet] = result$resultFinalModel$gamma[
        2:length(result$resultFinalModel$gamma)]
    }
    result$resultFinalModel$gamma = gamma
    result = result$resultFinalModel
  }
  
  gamma = result$gamma
  logLikAPP = result$logLik
  logLikMAP = result$logLikMAP
  PIPsFromBFsAPP = result$PIPsFromBFs
  PIPsFromBFsMAP = result$PIPsFromBFsMAP
  
  # standardized gamma
  if (hyperParamSelection != "topK") {
    standardDeviation = apply(annotation, 2, sd)
    gammaStandardized = gamma * standardDeviation
    gamma = cbind(gamma, gammaStandardized)
  } else {
    if (is.null(result$topKIndice)) {
      gamma = cbind(gamma, 0)
    } else {
      standardDeviation = apply(annotation[, result$topKIndice + 1, drop = F], 
                                2, sd)
      gammaStandardized = gamma * c(0, standardDeviation)
      gamma = cbind(gamma, gammaStandardized)
    }
  }
  
  outputResult(outputPrefix, gamma, PIPsFromBFsAPP, PIPsFromBFsMAP, 
               logLikAPP, logLikMAP, lociListFile)
  if (hyperParamSelection == "topK") {
    write.table(result$topKIndice, paste0(outputPrefix, ".topk"),
                row.names = F, col.names = F, quote = F)
    write.table(result$statsEachAnnotation, 
                paste0(outputPrefix, ".stats"),
                row.names = F, col.names = F, quote = F)
  }
  
  if (dim(annotation)[2] == 2) { 
    plotPIPVsAnnotation(PIPsFromBFsMAP, annotation[, 2], outputPrefix)
  }
  
  endTime = proc.time()
  print("total time used")
  print(endTime - startTime)
  timeFile = paste0(outputPrefix, ".time")
  file.create(timeFile)
  usedTime = endTime - startTime
  cat(names(usedTime), "\n", file = timeFile)
  cat(usedTime, "\n", file = timeFile, append = T)
  
  if (!keepBF) {
    unlink(tempDir, recursive=TRUE)
  }
  sink(NULL)
} 

# Calculate Bayes factors for multiple loci
# 
# @param annotaton
# @keywords caviarbf
# @examples
# see CVPartitionOnLoci
CAVIARBFMultiLoci = function(workingDir, lociListFile, LDSuffix = ".LD",
                             priorType01, priorValue,
                             nSample, maxCausal, exact = T, eps = 0,
                             useIdentityMatrix = F) {
  filePrefixes = as.vector(as.matrix(read.table(lociListFile, as.is = T)))
  dataDir = dirname(lociListFile)
  BFFileList = NULL
  for (i in 1:length(filePrefixes)) {
    BFFile = paste0(workingDir, "/", filePrefixes[i], "_l", maxCausal, ".bf")
    zFile = paste0(dataDir, "/", filePrefixes[i])
    rFile = paste0(dataDir, "/", filePrefixes[i], LDSuffix)
#     CAVIARBFWrapper(zFile, rFile, priorType01, priorValue, nSample, maxCausal,
#                     BFFile, exact, eps, useIdentityMatrix)
#     print(zFile)
#     print(rFile)
#     print(priorType01)
#     print(priorValue)
#     print(nSample)
#     print(maxCausal)
#     print(BFFile)
#     print(exact)
#     print(eps)
#     print(useIdentityMatrix)
    calculateBayesFactorsRcpp(zFile, rFile, priorType01, priorValue, nSample, 
                          maxCausal, BFFile, exact, eps, useIdentityMatrix)
    BFFileList = c(BFFileList, BFFile)
  }
  return(BFFileList)
}

# create a BFModelSearch object from each Bayes fator file
# 
# @param BFFileList: a vector of all Bayes factor filenames
# @param lociID: a vector of the loci IDs for each SNP
# @keywords caviarbf
# @examples
# see CVPartitionOnLoci
createBFModels = function(BFFileList, lociID) {
  uniqueLociID = unique(lociID)
  nLoci = length(uniqueLociID)
  BFModelList = list()
  for (i in 1:nLoci) {
    BFFile = BFFileList[i]
    locusIndice = which(lociID == uniqueLociID[i])
    nSNPs = length(locusIndice)
    BFModelList[[i]] = new(BFModelSearch)
    BFModelList[[i]]$loadBayesFactorData(BFFile, nSNPs)
  }
  
  return(BFModelList)
}

# Fine mapping causal variants incorporating functional annotations and using
# loci-based cross-validation to select the best penalty parameters
# 
# @param annotaton
# @keywords caviarbf
CVPartitionOnLoci = function(annotation, lociID,
                             BFModelList, nSample, PIPFilePrefix,
                             gammaInit, fittingMethod, lambda, alpha,
                             maxIter, deltaLik, trueHessian, posteriorUsage,
                             workingDir, nFold, pairTestThreshold = NULL) {
  # Assume lambda is in acsending order and favor large penalization 
  # when the cv logLik is the same
  
  #browser()
  
  nGamma = dim(annotation)[2]
  
  # partition data
  nLoci = length(unique(lociID))
  nFold = min(nLoci, nFold)
  
  foldID = sample(rep(1:nFold, length.out = nLoci))
  
  # cross validation to select the best parameter
  nLambda = length(lambda)
  if (fittingMethod == "glmnetENETMin") {
    nAlpha = length(alpha)
  } else {
    alpha = NA
    nAlpha = 1
  }
  
  if (nFold == 1 && (nLambda != 1 || nAlpha != 1)) {
    stop("only 1 locus! no way to do cross validation by partitioning on loci")
  }
  
  maxLogLik = -Inf
  bestLambda = NULL
  bestAlpha = NULL
  if (nLambda == 1 && nAlpha == 1) {
    bestLambda = lambda
    bestAlpha = alpha
  } else {
    logLike = matrix(NA, nLambda, nAlpha)
    foldLogLik = array(NA, dim = c(nLambda, nAlpha, nFold))
    maxFoldLogLik = NULL
    for (i in 1:nLambda) {
      for (j in 1:nAlpha) {
        print(paste("## lambda =", lambda[i], "alpha =", alpha[j]))
        result = CVLikelihoodOfEachHyperparameter(lambda[i], alpha[j],
                                annotation, lociID, 
                                BFModelList,
                                gammaInit, fittingMethod, 
                                maxIter, deltaLik, trueHessian, posteriorUsage,
                                workingDir, foldID) 
        logLike[i, j] = result$cvlogLik
        foldLogLik[i, j, ] = result$foldLogLik
        if (result$cvlogLik >= maxLogLik) { # favor large regularization
          bestLambda = lambda[i]
          bestAlpha = alpha[j]
          maxLogLik = result$cvlogLik
          maxFoldLogLik = result$foldLogLik
        }
        print(paste("lambda =", lambda[i], "alpha =", alpha[j], 
                    "cvlogLik =", result$cvlogLik))
        print(paste("current best: lambda =", bestLambda, "alpha =", bestAlpha,
                    "cvlogLik =", maxLogLik))
      }
    }
    
    # browser()
    if (!is.null(pairTestThreshold)) {
      # find the best lambda and alpha which has similar cv logLik with 
      # the best but achieve higher regularization. This has similar idea as 
      # that using one standard deviation to search the best
      bestLogLik = NA
      for (i in 1:nLambda) {
        for (j in 1:nAlpha) {
          currentFoldLogLik = foldLogLik[i, j, ]
          pvalue = wilcox.test(maxFoldLogLik, currentFoldLogLik, 
                               paired=T, , alternative="g")$p.value
          if (!is.na(pvalue) && !(pvalue < pairTestThreshold)) {
            bestLambda = lambda[i]
            bestAlpha = alpha[j]
            bestLogLik = sum(currentFoldLogLik)
          }
        }
      }
      print(paste("best after paired test: lambda =", 
                  bestLambda, "alpha =", bestAlpha,
                  "cvlogLik =", bestLogLik))
    }
  }
  
  #  train on the full data with the best hyperparameters with exact predict
  #  in glmnet
  modelFromFullData = EMIteration(annotation, lociID, BFModelList, 
                                  gammaInit, fittingMethod,
                                  maxIter, deltaLik, workingDir, 
                                  c(bestLambda, bestAlpha), PIPFilePrefix,
                                  exact = T)
  gamma = modelFromFullData$gamma
  
  # MAP
  priorProbMAP = c(1 / (1 + exp(-(annotation %*% gamma))))
  ResultFromBFs = calculatePIPsFromBFModels(priorProbMAP, lociID, 
                                            BFModelList, paste0("full"))
  logLikMAP = sum(ResultFromBFs$likelihood)    
  PIPsFromBFsMAP = ResultFromBFs$PIPsFromBFs
  
  # approximate integral
  # "approximate". The "sampling" option is not used due to time cost
  if (trueHessian) {
    runTime = system.time({
      jacobianMatrix = jacobian(PIPAndLikelihoodOfGamma, gamma, 
                                annotation = annotation, 
                                lociID = lociID, BFModelList = BFModelList, 
                                fittingMethod = fittingMethod, 
                                methodArgs = bestLambda)
    })
    print(runTime)
    print("derivative of the likelihood of gammas")
    print(jacobianMatrix[1, ])
    J = jacobianMatrix[2:dim(jacobianMatrix)[1], ]
  } else {
    J = NULL
  }
  logLikAPP = NA
  PIPsFromBFsAPP = NA
  # calculation of app results are commented
#   A = calculateHessianOfNegLik(annotation, gamma, J, bestLambda,
#                                fittingMethod = fittingMethod, bestAlpha)
#   priorProb = predictivePrior(annotation, gamma, A, annotation)
#   ResultFromBFs = calculatePIPsFromBFModels(priorProb, lociID, 
#                                             BFModelList, paste0("full"))
#   logLikAPP = sum(ResultFromBFs$likelihood)    
#   PIPsFromBFsAPP = ResultFromBFs$PIPsFromBFs
  
  return(list(gamma = gamma, PIPsFromBFs = PIPsFromBFsAPP, logLik = logLikAPP,
              PIPsFromBFsMAP = PIPsFromBFsMAP, logLikMAP = logLikMAP,
              cvLogLik = maxLogLik))
}

singleAnnotationTest = function(annotation, lociID,
                              BFModelList, nSample, PIPFilePrefix,
                              gammaInit, fittingMethod, 
                              maxIter, deltaLik, trueHessian, 
                              workingDir, scale = F) {
  # The first column contains all 1s, which is the null model
  nAnnotations = dim(annotation)[2] - 1
  likelihoodAll = numeric(nAnnotations + 1)
  gammaEachAnnotation = matrix(NA, nAnnotations + 1, 2)
  #PIPPerAnnotation = matrix(NA, dim(annotation)[1], nAnnotations + 1)
  for (i in 1:(nAnnotations + 1)) {
    indice = NULL
    if (i == 1) {
      indice = 1
    } else {
      indice = c(1, i)
    }
    if ((i > 1) && (var(annotation[, i]) == 0)) {
      likelihoodAll[i] = likelihoodAll[1]
      gammaEachAnnotation[i, ] = gammaEachAnnotation[1, ]
    } else {
      annotationColumns = annotation[, indice, drop = F]
      if (scale == T && i > 1) {
        annotationColumns[, 2] = scale(annotationColumns[, 2])
      }
      model = EMIteration(annotationColumns, 
                          lociID, BFModelList, 
                          gammaInit[indice], fittingMethod,
                          maxIter, deltaLik, workingDir, 
                          c(0, 0), NULL, exact = T)
      gammaEachAnnotation[i, 1:length(model$gamma)] = model$gamma
      likelihoodAll[i] = model$logLik
      
      # fitted PIP per annotation
      # PIPPerAnnotation[, i] = model$PIPsFromBFs
    }
  }
  twiceLoglikelihoodRatio =  2 * (likelihoodAll[2:dim(annotation)[2]] - 
                                    likelihoodAll[1])
  pvalue = pchisq(twiceLoglikelihoodRatio, df = 1, lower.tail = F)
  statsEachAnnotation = cbind(c(0, twiceLoglikelihoodRatio), 
                              gammaEachAnnotation,
                              c(NA, pvalue))
  return(statsEachAnnotation)
}

topKBasedonMaximumLikelihood = function(annotation, lociID, 
                              BFModelList, nSample, PIPFilePrefix,
                              gammaInit, fittingMethod, 
                              maxIter, deltaLik, trueHessian, 
                              posteriorUsage, workingDir, 
                              K, rThreshold, pvalueThreshold,
                              scale = F) {
  # The first column contains all 1s, which is the null model
  stopifnot(dim(annotation)[2] >= 2)
  statsEachAnnotation = singleAnnotationTest(annotation, lociID,  
                                BFModelList, nSample, PIPFilePrefix,
                                gammaInit, fittingMethod, 
                                maxIter, deltaLik, trueHessian, 
                                workingDir)
  twiceLoglikelihoodRatio = statsEachAnnotation[2:dim(statsEachAnnotation)[1], 
                                                1]
  corAnnotation = cor(annotation[, 2:dim(annotation)[2], drop = F])
  topKIndice = selectTopKBasedonLogLikelihoodDiff(twiceLoglikelihoodRatio, 
                                            corAnnotation,
                                            K, rThreshold, pvalueThreshold)
  annotationSelected = annotation[, c(1, topKIndice + 1), drop = F]
  
  # use the selected annotations
  modelFromSelected = EMIteration(annotationSelected, 
                                  lociID, BFModelList, 
                                  gammaInit[c(1, topKIndice + 1)], 
                                  fittingMethod,
                                  maxIter, deltaLik, workingDir, 
                                  c(0, 0), PIPFilePrefix, exact = T)
  gamma = modelFromSelected$gamma
  
#   # MAP
#   priorProbMAP = c(1 / (1 + exp(-(annotationSelected %*% gamma))))
#   ResultFromBFs = calculatePIPsFromBFModels(priorProbMAP, lociID, 
#                                             BFModelList, "selected")
#   logLikMAP = sum(ResultFromBFs$likelihood)    
#   PIPsFromBFsMAP = ResultFromBFs$PIPsFromBFs
  
  logLikMAP = modelFromSelected$logLik 
  PIPsFromBFsMAP = modelFromSelected$PIPsFromBFs

  return(list(gamma = gamma, PIPsFromBFs = NA, logLik = NA,
              PIPsFromBFsMAP = PIPsFromBFsMAP, logLikMAP = logLikMAP,
              topKIndice = topKIndice, 
              statsEachAnnotation = statsEachAnnotation)
         )
}

plotPIPVsAnnotation = function(PIPPerAnnotation, annotation, 
                                   outputPrefix) {
    pdfFile = paste0(outputPrefix, ".plot.pdf")
    pdf(pdfFile)
    plot(annotation, PIPPerAnnotation, 
         xlab = paste("annotation"), ylab = "fitted PIP")
    dev.off()
}

plotMultPIPVsAnnotation = function(PIPPerAnnotation, annotation, 
                               outputPrefix) {
  nAnnotation = dim(annotation) - 1 # exclude the 1s
  for (i in 1:nAnnotation) {
    pdfFile = paste0(outputPrefix, "_anno", i, ".pdf")
    pdf(pdfFile)
    plot(annotation[, i + 1], PIPPerAnnotation[, i + 1], 
         xlab = paste("annotation", i), ylab = "fitted PIP")
    dev.off()
  }
}

selectTopKBasedonLogLikelihoodDiff = function(twiceLoglikelihoodRatio, 
                                              corAnnotation,
                                              K, rThreshold, pvalueThreshold) {
  # select the top annotations based on the loglikelihood difference
  
  # input
  # twiceLoglikelihoodRatio: twice the log likelihood difference between the 
  #                          model including the annotation and the model 
  #                          with only the intercetp
  # corAnnotation: the correlation matrix between annotations
  
  # output
  # the selected indices staring with 1
  
  rank = order(twiceLoglikelihoodRatio, decreasing = T) 
  
  # select the top K annotations
  tsThreshold = qchisq(pvalueThreshold, df = 1, lower.tail = F)
  
  topKIndice = NULL
  for (i in 1:length(rank)) {
    candidateIndex = rank[i]
    if (twiceLoglikelihoodRatio[candidateIndex] >= tsThreshold) {
      if (is.null(topKIndice)) {
        topKIndice = candidateIndex
      } else { # check correlations
        correlations = corAnnotation[candidateIndex, topKIndice]
        if (all(abs(correlations) <= rThreshold)) {
          topKIndice = c(topKIndice, candidateIndex)
          if (length(topKIndice) == K) {
            break
          }
        }
      }
    }
  }
  
  return(topKIndice)
}

CVPartitionOnSubjects = function() {
  # TODO: partition on subjects. This will be useful if there is only one locus
}
