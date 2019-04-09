require(caviarbf)

main = function() {
  lociListFile = "../example/multiLoci/CAVIARBF/00000001"
  BFFileListFile = NA
  maxCausal = 3
  priorType01 = 0
  priorValue = 0.1
  exact = T
  priorProb = 0
  nSample = 2000
  fittingMethod = "glmnetLASSOMin"
  hyperParamSelection = "cv"
  nFold = 5
  lambda = c(2^seq(from = -15, by = 2, to = 5), 
             1e2, 1e3, 1e4, 1e5, 1e6)
  alpha = c(0, 0.2, 0.3, 0.5, 0.7, 0.8, 1)
  
  dataDir = normalizePath(dirname(lociListFile))
  
  for (eps in c(0, 0.2)) {
    workingDir = paste0(dataDir, "/priorType", priorType01, "_", 
                        paste(priorValue, collapse = "_"),
                        "_exact", "_eps", eps, "_", fittingMethod,
                        "_", hyperParamSelection)
    outputPrefix = paste0(workingDir, "/", basename(lociListFile), 
                          "_l", maxCausal)
    caviarbfFineMapping(lociListFile = lociListFile,  maxCausal = maxCausal, 
                        nSample = nSample, priorType01 = priorType01, 
                        priorValue = priorValue, exact = exact, eps = eps,  
                        BFFileListFile = BFFileListFile,  
                        fittingMethod = fittingMethod, 
                        hyperParamSelection = hyperParamSelection,
                        nFold = nFold, lambda = lambda, alpha = alpha, 
                        outputPrefix = outputPrefix, 
                        overwriteExistingResults = T)
  }
  
  K = 10
  rThreshold = 0.4
  pvalueThreshold = 0.1
  hyperParamSelection = "topK"
  eps = 0
  workingDir = paste0(dataDir, "/priorType", priorType01, "_", 
                      paste(priorValue, collapse = "_"),
                      "_exact", "_eps", eps,
                      "_", hyperParamSelection)
  outputPrefix = paste0(workingDir, "/", basename(lociListFile), 
                        "_l", maxCausal)
  caviarbfFineMapping(lociListFile = lociListFile,  maxCausal = maxCausal, 
                      nSample = nSample, priorType01 = priorType01, 
                      priorValue = priorValue, exact = exact, eps = eps,  
                      BFFileListFile = BFFileListFile,  
                      fittingMethod = fittingMethod, 
                      hyperParamSelection = hyperParamSelection,
                      K = K, rThreshold = rThreshold, 
                      pvalueThreshold = pvalueThreshold,
                      lambda = lambda, alpha = alpha, 
                      outputPrefix = outputPrefix, 
                      overwriteExistingResults = T)
}

main()