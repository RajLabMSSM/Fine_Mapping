#ifndef CALCULATEBAYESFACTORS_H
#define CALCULATEBAYESFACTORS_H

#include <string>
#include <vector>

using namespace std;

void calculateBayesFactors(string zFilename, \
  						   string corFilename, \
  						   int priorType, \
  						   vector<double> priorValues, \
  						   int nSample, \
  						   int maxCausal, \
  						   string outputFile, \
  						   bool exact, \
  						   double eps, \
						   bool useIdentityMatrix);
						   
#endif // CALCULATEBAYESFACTORS_H