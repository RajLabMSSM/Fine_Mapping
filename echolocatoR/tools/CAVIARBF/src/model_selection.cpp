#include <assert.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include <ctime>
#include <set>
#include <map>
#include <list>
#include <deque>
#include <iomanip>
#include <stdlib.h>
#include <fstream>

#include "model_selection.h"
#include "CombinationTable.h"
#include "bf_io.h"

using namespace std;

void preprocessPrior(std::vector<double>& allSNPPriors) {
	for (unsigned int i = 0; i < allSNPPriors.size(); i++) {
		if (1 - allSNPPriors[i] < 1e-12) {
			allSNPPriors[i] = 1 - 1e-12;
		}
	}
}

ostream& operator << (ostream& os, ExhaustiveSelection& es) {
	vector<SNPSetWithProbability> models = es.models;
	for (vector<SNPSetWithProbability>::iterator it = models.begin(); \
		  it != models.end(); it++) {
		SNPSetWithProbability& model = *it;
		os << std::scientific << model.probability;
	    for (vector<int>::iterator it2 = model.SNPSet.begin(); \
			   it2 != model.SNPSet.end(); it2++) {
			os << " " << *it2;
		}
		os << endl;
	}
	return os;
}

ostream& operator << (ostream& os, StepwiseSelection& ss) {
	vector<double>::iterator iterProbability = \
	           ss.probabilityOfEachStep.begin();
	for (vector<int>::iterator it = ss.selectedOfEachStep.begin(); \
		  it != ss.selectedOfEachStep.end(); it++) {
		os << *it << " " << std::scientific << *iterProbability << endl;
		iterProbability++;
	}
	return os;
}

void combineExhaustiveStepwise(ExhaustiveSelection& SNPList, \
						     StepwiseSelection& exhaustivestepwiseResult, \
							 ExhaustiveSelection& SNPList2) {
	SNPList2 = SNPList;
	int nExhaustive = SNPList.models.size();
	int nSNP = exhaustivestepwiseResult.selectedOfEachStep.size();
	SNPSetWithProbability SNPSetWithProb = SNPList2.models[nExhaustive - 1];
	for (int i = nExhaustive + 1; i <= nSNP; i++) {
		SNPSetWithProb.SNPSet.push_back( \
		        exhaustivestepwiseResult.selectedOfEachStep[i - 1]);
		SNPSetWithProb.probability = \
				 exhaustivestepwiseResult.probabilityOfEachStep[i - 1];
		SNPList2.models.push_back(SNPSetWithProb);
	}
}

int sumBoolArray(bool* array, int length) {
	int nTrue = 0;
	for (int i = 0; i < length; i++) {
		if (array[i]) {
			nTrue++;
		}
	}
	return(nTrue);
}

void normalizeInputProbPrior(vector<double>& prior) {
	double sumPrior = 0;
	for (unsigned int i = 0; i < prior.size(); i++) {
		sumPrior += prior[i];
	}
	for (unsigned int i = 0; i < prior.size(); i++) {
		prior[i] /= sumPrior;
	}
}

int locateBIMBAMBFIndexUsingCoLex(unsigned int n, int* SNPVector, \
									  unsigned int k, \
									  CombinationTable& combinationTable) {
  // This function uses the relationship between lexicographic ordering and
  // co-lex ordering to calculate the rank of lexicographic ordering

  // input
  // n: the total number of SNPs
  // SNPVector: the causal SNP ID in the model in 1 to n
  // k: the length of SNPVector
  // combinationTable: precomputed combination numbers

  // output
  // The row number of the Bayes factor (header and comments in the file are
  // not counted

  unsigned long long totalRank = combinationTable(n, k) - 1;
  unsigned long long rankColex = 0;
  for (unsigned int i = 1; i <= k; i++) {
    rankColex = rankColex + combinationTable(n - SNPVector[i - 1], k + 1 - i);
  }
  // added 1 because here the rank starts with 1
  unsigned long long rank = totalRank - rankColex + 1;
  for (unsigned int i = 1; i <= k - 1; i++) {
    rank = rank + combinationTable(n, i);
  }
  return(rank);
}

int locateBIMBAMBFIndexUsingCoLex(unsigned int n, int* SNPVector, \
									  unsigned int k, bool* mask, unsigned int t, \
									  CombinationTable& combinationTable) {
  // This function uses the relationship between lexicographic ordering and
  // co-lex ordering to calculate the rank of lexicographic ordering

  // input
  // n: the total number of SNPs
  // SNPVector: the causal SNP ID in the model in 1 to n
  // k: the length of SNPVector
  // mask: the mask indicating which SNP is in the model. This is used to
  //       create subsets of a SNP set
  // t: the total number of SNPs with true mask values
  // combinationTable: precomputed combination numbers

  // output
  // The row number of the Bayes factor (header and comments in the file are
  // not counted

  unsigned long long totalRank = combinationTable(n, t) - 1;
  unsigned long long rankColex = 0;
  int count = 0; // the real count of SNPs with a true mask value
  for (unsigned int i = 1; i <= k; i++) {
	  if (mask[i - 1]) {
		count++;
		rankColex += combinationTable(n - SNPVector[i - 1], t + 1 - count);
	  }
  }
  // added 1 because here the rank starts with 1
  unsigned long long rank = totalRank - rankColex + 1;
  for (unsigned int i = 1; i <= t - 1; i++) {
    rank = rank + combinationTable(n, i);
  }
  return(rank);
}

void calculatePriors(BayesFactorData& bfData, double prior, \
					double* modelPrior, double& nullModelPrior, int printLevel) {
	// calculate priors of models

	// input
	// prior:  prior for each SNP
	// modelPrior: the returned log10 priors
	// nullModelPrior: the prior of the null model,
	//                 i.e., no SNPs associated

	int count = 0;
	int nSNP = bfData.nSNP;
	double priorValue;
	int maxInModel = bfData.maxInModel;
	int numberOfCausal = 1;
	vector<double> accuProbAll(maxInModel + 1, 0);
	vector<double> totalPriorAll(maxInModel + 1, 0);
	vector<double> priorValueAll(maxInModel + 1, 0);
	vector<int> nModelPerSizeAll(maxInModel + 1, 0);

	if (printLevel > 0) {
		cout << "\nmodel_size\ttotal_prior\taccumulated_probability\n";
		cout << std::scientific;
	}
	//std::cout.precision(6);
	double epsilon = 1e-12;
	int nModelPerSize = 0;
	for (numberOfCausal = 0; numberOfCausal <= maxInModel; numberOfCausal++) {
		double totalPrior = 0;
		nModelPerSize = CombinationTable::choose(nSNP, numberOfCausal);
		if (fabs(prior - 0) < epsilon || (prior > 0 && prior < 1)) {
			// p(l) ~ binom(l, prior)
			double prob;
			if (fabs(prior - 0) < epsilon) {
				prob = 1.0 / nSNP;
			} else {
				prob = prior;
			}
			priorValue = log10(prob) * numberOfCausal + \
						 log10(1 - prob) * (nSNP - numberOfCausal);
		} else if (fabs(prior - 1) < epsilon) {
		    // p(l) proportional to (1/2)^l
			if (numberOfCausal > 0) {
				priorValue = log10(0.5) * numberOfCausal - log10(nModelPerSize);
			} else {
				priorValue = log10(0.5) * maxInModel;
			}
		} else {
			throw BadParameter();
		}
		totalPrior += (pow(10, priorValue) * nModelPerSize);

		// save the numbers
		nModelPerSizeAll[numberOfCausal] = nModelPerSize;
		priorValueAll[numberOfCausal] = priorValue;
		totalPriorAll[numberOfCausal] = totalPrior;
		if (numberOfCausal == 0) {
			accuProbAll[numberOfCausal] = totalPrior;
		} else {
			accuProbAll[numberOfCausal] = accuProbAll[numberOfCausal - 1] + \
			                              totalPrior;
		}
	}

	nullModelPrior = std::pow(10, priorValueAll[0]);

	for (numberOfCausal = 0; numberOfCausal <= maxInModel; numberOfCausal++) {
		nModelPerSize = nModelPerSizeAll[numberOfCausal];
		if (printLevel > 0) {
			cout << numberOfCausal << "\t\t" << totalPriorAll[numberOfCausal] \
		     << "\t\t" << accuProbAll[numberOfCausal] << endl;
		}
		if (numberOfCausal > 0) {
			for (int j = count; j < count + nModelPerSize; j++) {
				modelPrior[j] = priorValueAll[numberOfCausal];
			}
			count += nModelPerSize;
		}
	}
}

void calculatePriors(BayesFactorData& bfData, vector<double>& prior, \
					double* modelPrior, double& nullModelPrior, \
					int printLevel) {
	// calculate priors of models with different prior for each SNPs

	// input
	// prior:  prior for each SNP
	// modelPrior: the returned log10 priors
	// nullModelPrior: the prior of the null model,
	//                 i.e., no SNPs associated

	int maxInModel = bfData.maxInModel;
	int numberOfCausal;
	vector<double> accuProbAll(maxInModel + 1, 0);
	vector<double> totalPriorAll(maxInModel + 1, 0);
	
	if (printLevel > 0) {
		cout << "\nmodel_size\ttotal_prior\taccumulated_probability\n";
		cout << std::scientific;
	}

	// nullModelPrior = 1;
	double nullModelPriorLog10 = 0;
	for (unsigned int i = 0; i < prior.size(); i++) {
		nullModelPriorLog10 += log10(1 - prior[i]);
//		// <debug>
//		cout << "1 - prior[i] " << 1 - prior[i] << endl;
//		cout << "nullModelPriorLog10 " << nullModelPriorLog10 << endl;
//		// </debug>
	}
	nullModelPrior = std::pow(10, nullModelPriorLog10);
	totalPriorAll[0] = nullModelPrior;
//	// <debug>
//	cout << "nullModelPriorLog10 " << nullModelPriorLog10 << endl;
//	cout << "nullModelPrior " << nullModelPrior << endl;
//	cout << "log10(nullModelPrior) " << log10(nullModelPrior) << endl;
//	// </debug>

	// non-null model prior
	int nSNPInModel = 0;
	double priorOfAModel = 0;
	for (int i = 0; i < bfData.nModel; i++) {
		nSNPInModel = bfData.nSNPInModel[i];
		priorOfAModel = log10(nullModelPrior);
		for (int j = 0; j < nSNPInModel; j++) {
			int SNPID = bfData.SNPID(i + 1, j + 1);
			priorOfAModel += log10(prior[SNPID - 1]) - \
			                 log10(1 - prior[SNPID - 1]);
//			// <debug>
//			cout << prior[SNPID - 1] << " ";
//			// </debug>
		}
//		// <debug>
//		cout << priorOfAModel << endl;
//		// </debug>
		modelPrior[i] = priorOfAModel;
		totalPriorAll[nSNPInModel] += std::pow(10, priorOfAModel);
	}

	for (numberOfCausal = 0; numberOfCausal <= maxInModel; numberOfCausal++) {
		if (numberOfCausal == 0) {
			accuProbAll[numberOfCausal] = totalPriorAll[numberOfCausal];
		} else {
			accuProbAll[numberOfCausal] = accuProbAll[numberOfCausal - 1] + \
			                              totalPriorAll[numberOfCausal];
		}
		if (printLevel > 0) {
			cout << numberOfCausal << "\t\t" << totalPriorAll[numberOfCausal] \
		     << "\t\t" << accuProbAll[numberOfCausal] << endl;
		}
	}
}

// minus the largest log BF to void overflow in later calculation
void minusTheLargestBayesFactor(BayesFactorData& bfData, \
                                double& largestLogBF) {
	// find the largest Bayes factor
	largestLogBF = bfData.bayesFactor[0];
	for (int i = 1; i < bfData.nModel; i++) {
		if (bfData.bayesFactor[i] > largestLogBF) {
			largestLogBF = bfData.bayesFactor[i];
		}
	}

	// minus the largest Bayes factor
	for (int i = 0; i < bfData.nModel; i++) {
		bfData.bayesFactor[i] -= largestLogBF;
	}
}

void calculateModelScoreAndSum(BayesFactorData& bfData, double* modelPrior, \
				double* modelScore, double nullModelScore, double& totalScore, \
				int printLevel) {
	totalScore = nullModelScore;
	int maxInModel = bfData.maxInModel;
	int nSNP = bfData.nSNP;
	double* posterior = new double[maxInModel + 1];
	double accu_posterior;
	int numberOfCausal = 1;
	int count = 0;
	posterior[0] = nullModelScore;
	for (; numberOfCausal <= maxInModel; numberOfCausal++) {
		int nModelPerSize = CombinationTable::choose(nSNP, numberOfCausal);
		posterior[numberOfCausal] = 0;
		for (int j = count; j < count + nModelPerSize; j++) {
			modelScore[j] = std::pow(10, bfData.bayesFactor[j] + modelPrior[j]);
			totalScore += modelScore[j];
			posterior[numberOfCausal] += modelScore[j];
		}
		count += nModelPerSize;
	}

	if (printLevel > 0) {
		cout << "\nmodel_size\ttotal_posterior\taccumulated_probability\n";
	}
	numberOfCausal = 0;
	accu_posterior =  0;
	for (; numberOfCausal <= maxInModel; numberOfCausal++) {
		double  posterior_prob = posterior[numberOfCausal] / totalScore;
		accu_posterior += posterior_prob;
		if (printLevel > 0) {
			cout << numberOfCausal << "\t\t" \
		     << posterior_prob << "\t\t" << accu_posterior << endl;
		}
	}

	delete [] posterior;
}

double logLikelihood(double totalScore, double& largestLogBF) {
	double logLikelihoodUpToAConstant = log(totalScore) + \
	                                    largestLogBF * log(10);
	return logLikelihoodUpToAConstant;
}

double probabilityAtLeast1CausalVariant(double totalScore, \
                           double nullModelScore) {
	return 1 - nullModelScore / totalScore;
}

double globalBF(double totalScore, double nullModelPrior, \
                      double nullModelScore, double& largestLogBF) {
	double globalBayesFactor = log10(totalScore - nullModelScore) -  \
	                           log10(1 - nullModelPrior) + largestLogBF;
	return globalBayesFactor;
}

void outputStatistics(double totalScore, double nullModelPrior, \
                      double nullModelScore, double& largestLogBF, \
						ofstream& outputStream) {
	double logLikelihoodUpToAConstant = \
	              logLikelihood(totalScore, largestLogBF);
	double probabilityAtLeast1Causal =  \
				   probabilityAtLeast1CausalVariant(totalScore, nullModelScore);
	double globalBayesFactor = globalBF(totalScore, nullModelPrior, \
	                                    nullModelScore, largestLogBF);
	outputStream << std::scientific << std::setprecision(15) \
	         << logLikelihoodUpToAConstant << \
			  " # log(p(Data)) - log(p(Data|M0))" << endl;
	cout << logLikelihoodUpToAConstant << \
	             " # log(p(Data)) - log(p(Data|M0))" << endl;
	outputStream << probabilityAtLeast1Causal << \
	   " # Probability of at least 1 causal variant P(M_alternative)" << endl;
	cout << probabilityAtLeast1Causal << \
	   " # Probability of at least 1 causal variant P(M_alternative)" << endl;
	// renormalize the prior probabability of alternative models
	outputStream << globalBayesFactor << \
	       " # log10 Bayes factor of the region (M_alternative vs M0)" << endl;
	cout << globalBayesFactor << \
	       " # log10 Bayes factor of the region (M_alternative vs M0)" << endl;
}

void marginalProbability(BayesFactorData& bfData, double* modelScore, \
		           double totalScore, StepwiseSelection& marginalResult) {
	int nSNP = bfData.nSNP;
	vector<double> marginalProb(nSNP, 0);
	int nSNPInModel = 0;
	double score = 0;
	for (int i = 0; i < bfData.nModel; i++) {
		score = modelScore[i];
		nSNPInModel = bfData.nSNPInModel[i];
		for (int j = 0; j < nSNPInModel; j++) {
			int SNPID = bfData.SNPID(i + 1, j + 1);
			marginalProb[SNPID - 1] += score;
		}
	}

	// get the probability
	for (int i = 0; i < nSNP; i++) {
		marginalProb[i] /= totalScore;
	}

	// sort the SNPs based on marginal probability
	multimap<double, int> SNPIDProbMap;
	for (int i = 0; i < nSNP; i++) {
		SNPIDProbMap.insert(std::pair<double, int>(marginalProb[i], i + 1));
	}
	multimap<double, int>::reverse_iterator rit;
	for (rit = SNPIDProbMap.rbegin(); rit != SNPIDProbMap.rend(); ++rit) {
		marginalResult.selectedOfEachStep.push_back(rit->second);
		marginalResult.probabilityOfEachStep.push_back(rit->first);
	}
}

void exhaustiveSearch(BayesFactorData& bfData, double* modelScore, \
		double totalScore, int maxcausalSNP, ExhaustiveSelection& SNPList) {
	// exhaustive search to find the best combinations of SNP sets
	// up to maxcausalSNP

	// input
	// modelScore: bayes factor * prior
	// maxcausalSNP: search combinations from 1 to maxcausalSNP

	// output
	//

	unsigned long long nModel = bfData.nModel;
	int maxInModel = bfData.maxInModel;

	// calculate SNP set scores
	int count = 0;
	double* SNPSetScore = new double[nModel];
	int* causalSNPs = new int[maxcausalSNP];
	int indexWithoutTheFirstSNP;
	int indexWithTheFirstSNP;
	int* SNPIDInModel;
	int nSNP = bfData.nSNP;
	unsigned int nSubset = 0;
	bool* mask = NULL;
	int* nCausalSNP = NULL;
	bool* mask_row;
	int max_n = nSNP;
	int max_k = maxcausalSNP;
	CombinationTable combinationTable(max_n,  max_k);
	for (int causalSize = 1; causalSize <= maxcausalSNP; causalSize++) {
		int nModelPerSize = CombinationTable::choose(nSNP, causalSize);
		if (causalSize > 1) {
			nSubset = std::pow(2.0, causalSize - 1);
			mask = new bool[nSubset * causalSize];
			nCausalSNP = new int[nSubset];
			// initialize the mask
			for (unsigned int row = 0; row < nSubset; row++) {
				mask_row = mask + row * causalSize;
				mask_row[0] = true;
				for (int k = 1; k < causalSize; k++) {
					mask_row[k] = (((1U << (k - 1)) & row) != 0);
				}
				nCausalSNP[row] = sumBoolArray(mask_row, causalSize);
			}
		}

		for (int index = count; index < count + nModelPerSize; index++) {
			if (causalSize == 1) { // only one causal SNP
				SNPSetScore[index] = modelScore[index];
			} else {
				// use the sum of all subset of the casual SNP set
				SNPIDInModel = bfData.SNPIDInModel + index * maxInModel;
				memcpy(causalSNPs, SNPIDInModel, sizeof(int) * causalSize);
//				// <debug>
//				for (int i = 0; i < causalSize; i++) {
//					cout << " " << causalSNPs[i];
//				}
//				cout << endl;
//				// </debug>
				indexWithoutTheFirstSNP =  \
				   locateBIMBAMBFIndexUsingCoLex(nSNP, causalSNPs + 1, \
				                                 causalSize - 1, \
				                                 combinationTable);
//				// <debug>
//				cout << "indexWithoutTheFirstSNP: " << indexWithoutTheFirstSNP \
//				      << endl;
//				// </debug>
				SNPSetScore[index] = SNPSetScore[indexWithoutTheFirstSNP - 1];
				// add the score of rest combinations including the first SNP
				for (unsigned int k = 0; k < nSubset; k++) {
					mask_row = mask + k * causalSize;
					indexWithTheFirstSNP = locateBIMBAMBFIndexUsingCoLex(nSNP, \
					                       causalSNPs, causalSize, \
										      mask_row, nCausalSNP[k], \
											  combinationTable);
//					// <debug>
//				   cout << "indexWithTheFirstSNP: " << indexWithTheFirstSNP \
//				      << endl;
//				    // </debug>
					SNPSetScore[index] += modelScore[indexWithTheFirstSNP - 1];
				}
			}
		}
		count += nModelPerSize;

		if (causalSize > 1) {
			delete [] mask;
			delete [] nCausalSNP;
		}
	}

	// <debug>
//	cout << std::hex << \
//      *reinterpret_cast<unsigned __int64 *>(SNPSetScore + 17951) << endl;
//	cout << std::hex << \
//      *reinterpret_cast<unsigned __int64 *>(SNPSetScore + 18299)<< endl;
//	if (SNPSetScore[18299] < SNPSetScore[17951]) {
//		cout << "SNPSetScore[18299] < SNPSetScore[17951]" << endl;
//	} else if (SNPSetScore[18299] == SNPSetScore[17951]) {
//		cout << "SNPSetScore[18299] == SNPSetScore[17951]" << endl;
//	} else {
//		cout << "SNPSetScore[18299] > SNPSetScore[17951]" << endl;
//	}
	// </debug>

	// search the best SNP combinations for each size of models and report the
	// corresponding probabilities
	cout << "\n --exhaustive search-- \n";
	count = 0;
	for (int causalSize = 1; causalSize <= maxcausalSNP; causalSize++) {
		int nModelPerSize = CombinationTable::choose(nSNP, causalSize);
		double maxScore = SNPSetScore[count];
		unsigned long long indexMaxScore = count;
		for (int index = count + 1; index < count + nModelPerSize; index++) {
			if (SNPSetScore[index] > maxScore) {
				maxScore = SNPSetScore[index];
				indexMaxScore = index;
			}
		}
		double probability = maxScore / totalScore;
		SNPSetWithProbability bestSet;
		bestSet.probability = probability;
		SNPIDInModel = bfData.SNPIDInModel + indexMaxScore * maxInModel;
		for (int i = 0; i < causalSize; i++) {
			bestSet.SNPSet.push_back(SNPIDInModel[i]);
		}
		SNPList.models.push_back(bestSet);

		// output to console
		cout << probability;
		for (int i = 0; i < causalSize; i++) {
			cout << " " << bestSet.SNPSet[i];
		}
		cout << endl;

		count += nModelPerSize;
	}

	delete [] causalSNPs;
	delete [] SNPSetScore;
}

int next_combination(int* c, int R, int N) {
	// generate the next combination of size R from size N

	// input:
	// N: the total number of objects index by 0 ... N - 1
	// R: the total number of chosen objects
	// c: the selected R object indice

	//output:
	// 0: the end of the enumeration
	// a nonzero number: c stores the next combination

    int i = 0;
    while ( i < R - 1 && c[i + 1] == c[i] + 1) { // for each bump
       c[i] = i;                // fall back
		i++;
	}
    return N - ++c[i];             // push forward and verify
}

void greedySearch(BayesFactorData& bfData, double* modelScore, \
		           double totalScore, vector<int> selected, double rho, \
                  StepwiseSelection& stepwiseResult) {
	int nSNP = bfData.nSNP;
	int maxInModel = bfData.maxInModel;
	double* SNPInclusionScore = new double[nSNP];
	double* scoreOfEachStep = new double[nSNP];
	vector<int> selectedOfEachStep;
	vector<double> probabilityOfEachStep;
	set<int> selectedSet;   // keep an ordered set of selected
	double currentScore = 0;
	// initialize with model socres of individual SNPs
	for (int i = 0; i < nSNP; i++) {
		SNPInclusionScore[i] = modelScore[i];
	}
	int newSelected;
	double maxScore;
	int maxIndex;
	bool* selectedFlag = new bool[nSNP];
	memset(selectedFlag, false, nSNP * sizeof(bool));
	int max_n = nSNP;
	int max_k = maxInModel;
	CombinationTable combinationTable(max_n,  max_k);
	// stepwise selection
	if (selected.size() > 0) {
		cout << "\n --stepwise selection from exhaustive results-- \n";
		cout << "order of the first " << selected.size() \
		     << " SNPs are ignored!\n";
	} else {
		cout << "\n --stepwise selection-- \n";
	}
	cout << "step\tSNP\tprobability" << endl;
	for (int step = 1; step <= nSNP; step++) {
		if (step <= (int)selected.size()) {
			newSelected = selected[step - 1]; // force the selection
		} else {  // select the one with the highest SNPInclusionScore
			maxScore = -1;
			maxIndex = -1;
			for (int SNPID = 1; SNPID <= nSNP; SNPID++) {
				if (!selectedFlag[SNPID - 1]) {
				   if (SNPInclusionScore[SNPID - 1] > maxScore) {
					   maxIndex = SNPID;
					   maxScore = SNPInclusionScore[SNPID - 1];
				   }
				}
			}
			newSelected = maxIndex;
		}
		selectedSet.insert(newSelected);
		selectedOfEachStep.push_back(newSelected);
		selectedFlag[newSelected - 1] = true;
		currentScore += SNPInclusionScore[newSelected - 1];
		scoreOfEachStep[step - 1] = currentScore;
		double currentProbability = currentScore / totalScore;
		probabilityOfEachStep.push_back(currentProbability);
		cout << step << "\t" << newSelected << "\t" << \
		        currentProbability << endl;
		if (currentProbability > rho) {
			break;
		}
		// update the inclusion score for the next comparison
		int* combinationPosition = NULL;
		if (step > 1) {
			combinationPosition = new int[step - 1];
		}
		int* causalSNP = new int[step + 1];
		bool* initMask = new bool[step + 1];
		bool* mask = new bool[step + 1];
		unsigned long long index;
		for (int SNPID = 1; SNPID <= nSNP; SNPID++) {
			if(selectedFlag[SNPID -1]) {
				continue;
			}
			// cout << "update SNPID: " << SNPID << endl;
			memset(initMask, false, (step + 1) * sizeof(bool));
			memset(mask, false, (step + 1) * sizeof(bool));
			set<int> tempSelected = selectedSet;
			tempSelected.insert(SNPID);
			int i = 0;
			int p = 0;
			for (set<int>::iterator it = tempSelected.begin(); \
				                     it != tempSelected.end(); ++it) {
				 causalSNP[i] = *it;
				 if (causalSNP[i] != SNPID && causalSNP[i] != newSelected) {
					 if (step > 1) {
						 combinationPosition[p] = i;
						 p++;
					 }
				 } else {
					 initMask[i] = true;
				 }
				 i++;
			}
			// enumerate all combinations of setSize elements
			int maxSetSize = min(maxInModel - 2, int(step - 1));
			for (int setSize = 0; setSize <= maxSetSize; setSize++) {
				if (setSize == 0) {
					 index = locateBIMBAMBFIndexUsingCoLex(nSNP, \
								causalSNP, step + 1, \
								initMask, 2,
								combinationTable);
					SNPInclusionScore[SNPID - 1] += modelScore[index - 1];
				} else {
					int* combination = new int[setSize];
					// initial value
					for (int element = 0; element < setSize; element++) {
						combination[element] = element;
					}
					int result;
					int ID;
					while(true) {
						memcpy(mask, initMask, sizeof(bool) * (step + 1));
//					   cout << "combination of size: " << setSize << endl;
						for (int element = 0; element < setSize; element++) {
							ID = combinationPosition[combination[element]];
							mask[ID] = true;
//							cout << " " << causalSNP[ID];
						}
//						cout << endl;
//						cout << "causalSNP and mask" << endl;
//						for (int temp = 0; temp < step + 1; temp++) {
//							cout << " " << causalSNP[temp];
//						}
//						cout << endl;
//						for (int temp = 0; temp < step + 1; temp++) {
//							cout << " " << mask[temp];
//						}
//						cout << endl;
						index = locateBIMBAMBFIndexUsingCoLex(nSNP, \
										   causalSNP, step + 1, \
											  mask, setSize + 2, \
											  combinationTable);
						SNPInclusionScore[SNPID - 1] += modelScore[index - 1];
						result = next_combination(combination, setSize, step - 1);
						if (result == 0) {
							break;
						}
					}
					delete [] combination;
				}
			}
		}

		if (step > 1) {
			delete [] combinationPosition;
		}
		delete [] causalSNP;
		delete [] initMask;
		delete [] mask;
	}
	stepwiseResult.selectedOfEachStep = selectedOfEachStep;
	stepwiseResult.probabilityOfEachStep = probabilityOfEachStep;
	delete [] SNPInclusionScore;
	delete [] scoreOfEachStep;
	delete [] selectedFlag;
}


// This function does not achieve what I initially thought
void greedySearchMarginal(BayesFactorData& bfData, double* modelScore, \
					double totalScore, vector<int> selected, double rho, \
					StepwiseSelection& stepwiseResult) {
	int nSNP = bfData.nSNP;
	unsigned long long nModel = bfData.nModel;
	vector< list<int> > modelIndexPerSNP(nSNP);
	vector<double> contributedScore(nSNP, 0);
	deque< vector< list<int>::iterator > > indexLocationInList(nModel);
	vector<int> selectedOfEachStep;
	vector<double> probabilityOfEachStep;
	double currentScore = 0;
	int nSNPInModel;
	int SNPID;
	double score;

	cout << "\n --marginal stepwise selection-- \n";
	// scan all the models to fill in model indice, its location and score
	for (unsigned int i = 0; i < nModel; i++) {
		score = modelScore[i];
		nSNPInModel = bfData.nSNPInModel[i];
		for (int j = 0; j < nSNPInModel; j++) {
			SNPID = bfData.SNPID(i + 1, j + 1);
			// store model indice
			list<int>& modelIndex = modelIndexPerSNP[SNPID - 1];
			modelIndex.push_back(i);
			// store indice location
			list<int>::iterator& listIterator = (--modelIndex.end());
			indexLocationInList[i].push_back(listIterator);
			// add up scores
			contributedScore[SNPID - 1] += score;
		}
	}

	// <debug>
	for (int i = 0; i < nSNP; i++) {
		cout << setprecision(16) << log10(contributedScore[i]) << endl;
	}
	// </debug>

	vector<bool> selectedFlag(nSNP, false);
	for (int i = 0; i < nSNP; i++) {
		// pick the SNP with the largest score

		// <debug>
//		for (int k = 0; k < nSNP; k++) {
//			cout << setprecision(16) << contributedScore[k] << endl;
//		}
		// </debug>
		double maxScore = -1e300;
		int selectedSNPID = -1;
		for (int j = 0; j < nSNP; j++) {
			if (!selectedFlag[j] && contributedScore[j] > maxScore) {
				selectedSNPID = j + 1;
				maxScore = contributedScore[j];
			}
		}
		if (selectedSNPID == -1) {
			cout << "error in select the next SNP" << endl;
			exit(-1);
		}

		// update selection and probabilities
		selectedFlag[selectedSNPID - 1] = true;
		selectedOfEachStep.push_back(selectedSNPID);
		currentScore += maxScore;
		double currentProbability = currentScore / totalScore;
		probabilityOfEachStep.push_back(currentProbability);
		cout << i + 1 << "\t" << selectedSNPID << "\t" << \
		        currentProbability << endl;
		if (currentProbability > rho + 1e-6) {
			break;
		}
		// update model score and indices for each unselected SNPs
		int index;
		list<int>& modelIndex = modelIndexPerSNP[selectedSNPID - 1];
		for (list<int>::iterator it = modelIndex.begin();
				it != modelIndex.end(); it++) {
			// update indices and scores
			index = *it;
			score = modelScore[index];
			nSNPInModel = bfData.nSNPInModel[index];
			for (int j = 0; j < nSNPInModel; j++) {
				SNPID = bfData.SNPID(index + 1, j + 1);
				if (!selectedFlag[SNPID - 1]) {
					// remove model indice
					list<int>& modelIndexUnSelected =
					                       modelIndexPerSNP[SNPID - 1];
					list<int>::iterator& position =
					                       (indexLocationInList[index])[j];
					// cout << "remove " << *position << " from " << SNPID << endl;
					modelIndexUnSelected.erase(position);
					// subtract scores
					contributedScore[SNPID - 1] -= score;
				}
			}
		}
	}

	stepwiseResult.selectedOfEachStep = selectedOfEachStep;
	stepwiseResult.probabilityOfEachStep = probabilityOfEachStep;
}