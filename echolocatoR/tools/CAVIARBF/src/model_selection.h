#ifndef MODEL_SELECTION_H
#define MODEL_SELECTION_H

#include <vector>
#include <exception>
#include <fstream>

#include "bf_io.h"
#include "CombinationTable.h"

using namespace std;

struct SNPSetWithProbability {
	vector<int> SNPSet;
	double probability;
};

struct ExhaustiveSelection {
	vector<SNPSetWithProbability> models;
};

struct StepwiseSelection {
	vector<int> selectedOfEachStep;
	vector<double> probabilityOfEachStep;
};

struct BadParameter:public std::exception {
	const char* what() const throw () {
		return "bad input parameters";
	}
};

int locateBIMBAMBFIndexUsingCoLex(unsigned int n, int* SNPVector, \
									  unsigned int k, \
									  CombinationTable& combinationTable);
int locateBIMBAMBFIndexUsingCoLex(unsigned int n, int* SNPVector, \
									  unsigned int k, bool* mask, unsigned int t, \
									  CombinationTable& combinationTable);
void calculatePriors(BayesFactorData& bfData, double prior, \
					double* modelPrior, double& nullModelPrior, \
					int printLevel = 1);
void calculatePriors(BayesFactorData& bfData, vector<double>& prior, \
					double* modelPrior, double& nullModelPrior, \
					int printLevel = 1);
void minusTheLargestBayesFactor(BayesFactorData& bfData, 
					double& largestLogBF);
void calculateModelScoreAndSum(BayesFactorData& bfData, double* modelPrior, \
			double* modelScore, double nullModelScore, double& totalScore, \
			int printLevel = 1);
void marginalProbability(BayesFactorData& bfData, double* modelScore, \
		           double totalScore, StepwiseSelection& marginalResult);
void exhaustiveSearch(BayesFactorData& bfData, double* modelScore, \
		double totalScore, int maxcausalSNP, ExhaustiveSelection& SNPList);
void greedySearch(BayesFactorData& bfData, double* modelScore, \
		           double totalScore, vector<int> selected, double rho, \
                  StepwiseSelection& stepwiseResult);
//void greedySearchMarginal(BayesFactorData& bfData, double* modelScore, \
//					double totalScore, vector<int> selected, double rho, \
//					StepwiseSelection& stepwiseResult);
int sumBoolArray(bool* array, int length);
ostream& operator << (ostream& os, ExhaustiveSelection& es);
ostream& operator << (ostream& os, StepwiseSelection& ss);
void combineExhaustiveStepwise(ExhaustiveSelection& SNPList, \
						     StepwiseSelection& exhaustivestepwiseResult, \
							 ExhaustiveSelection& SNPList2);
void normalizeInputProbPrior(vector<double>&);
double logLikelihood(double totalScore, double& largestLogBF);
double probabilityAtLeast1CausalVariant(double totalScore, \
                           double nullModelPrior);
double globalBF(double totalScore, double nullModelPrior, \
                      double nullModelScore, double& largestLogBF);
void outputStatistics(double totalScore, double nullModelPrior, \
                      double nullModelScore, double& largestLogBF, \
					    ofstream& outputStream);
void preprocessPrior(std::vector<double>& allSNPPriors);

#endif