#include "InputArgument.h"
#include "calculateBayesFactors.h"

int main(int argc, char** argv) {
  InputArgumentCAVIARBF inputArgument;
  inputArgument.processArguments(argc, argv);

  calculateBayesFactors(inputArgument.zFilename, inputArgument.corFilename, 
         inputArgument.priorType, inputArgument.priorValues, inputArgument.nSample,
         inputArgument.maxCausal, inputArgument.outputFile, inputArgument.exact,
         inputArgument.eps, inputArgument.useIdentityMatrix);
  return 0;
}
