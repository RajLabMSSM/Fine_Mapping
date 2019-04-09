#include <Rcpp.h>
#include <RcppEigen.h>
#include <vector>

#include "calculateBayesFactors.h"

// [[Rcpp::depends(RcppEigen)]]

// [[Rcpp::export]]
int calculateBayesFactorsRcpp(std::string zFilename, 
  						   std::string corFilename, 
  						   int priorType, 
  						   std::vector<double> priorValues, 
  						   int nSample, 
  						   int maxCausal, 
  						   std::string outputFile, 
  						   bool exact, 
  						   double eps, 
						     bool useIdentityMatrix) {
  calculateBayesFactors(zFilename, corFilename, 
         priorType, priorValues, 
         nSample,
         maxCausal, outputFile, exact,
         eps, useIdentityMatrix);
  return 0;
} 