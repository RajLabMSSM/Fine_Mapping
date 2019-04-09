#include <Rcpp.h>

#include <cmath>

#include "model_selection.h"
#include "bf_io.h"

class BFModelSearch {
public:
    void loadBayesFactorData(Rcpp::String filename, int nSNPs);
	void calculateLikelihoodAndPIPsFromVector(vector<double> allSNPPriors);
	void calculateLikelihoodAndPIPs(Rcpp::NumericVector allSNPPriors);

	~BFModelSearch() {
		delete pbfData;
		delete [] modelPrior;
		delete [] modelScore;
	}

private:
    BayesFactorData* pbfData;
	double* modelPrior;
	double* modelScore;

public:
	double totalScore;
	double largestLogBF;
	double nullModelPrior;
	double nullModelScore;
	double logLikelihoodUpToAConstant;
	double probabilityAtLeast1Causal;
	double globalBayesFactor;
	vector<double> PIPs;
	vector<int> SNPIDs;
};

void BFModelSearch::loadBayesFactorData(Rcpp::String filename, int nSNPs) {
	pbfData = new BayesFactorData(filename.get_cstring(), nSNPs);
	minusTheLargestBayesFactor(*pbfData, largestLogBF);
	modelPrior = new double[pbfData->nModel];
	modelScore = new double[pbfData->nModel];
}

void BFModelSearch::calculateLikelihoodAndPIPsFromVector(\
                                  vector<double> allSNPPriors) {
	int printLevel = 0;
	preprocessPrior(allSNPPriors);
	calculatePriors(*pbfData, allSNPPriors, modelPrior, nullModelPrior, \
                    printLevel);
	nullModelScore = nullModelPrior * std::pow(10, -largestLogBF);
	calculateModelScoreAndSum(*pbfData, modelPrior, modelScore, \
									 nullModelScore, totalScore, printLevel);
	// <debug>
	// cout << "nullModelPrior:" << nullModelPrior << endl;
	// cout << "largestLogBF:" << largestLogBF << endl;
	// cout << "nullModelScore:" << nullModelScore << endl;
	// cout << "totalScore:" << totalScore << endl;
	// </debug>
	logLikelihoodUpToAConstant = logLikelihood(totalScore, largestLogBF);
	probabilityAtLeast1Causal =  \
				   probabilityAtLeast1CausalVariant(totalScore, nullModelScore);
	globalBayesFactor = globalBF(totalScore, nullModelPrior, \
	                                    nullModelScore, largestLogBF);

	// calculate PIPs
	StepwiseSelection marginalProbabilitySelection;
	marginalProbability(*pbfData, modelScore, totalScore, \
		                    marginalProbabilitySelection);
	PIPs = marginalProbabilitySelection.probabilityOfEachStep;
	SNPIDs = marginalProbabilitySelection.selectedOfEachStep;
}

void BFModelSearch::calculateLikelihoodAndPIPs(\
                               Rcpp::NumericVector allSNPPriors) {
	calculateLikelihoodAndPIPsFromVector(\
							Rcpp::as<vector<double> >(allSNPPriors));
}


RCPP_MODULE(BFModelSearchModule){
    using namespace Rcpp ;

    class_<BFModelSearch>( "BFModelSearch" )

    .constructor()

	.method( "loadBayesFactorData", &BFModelSearch::loadBayesFactorData)
    .method( "calculateLikelihoodAndPIPs", \
               &BFModelSearch::calculateLikelihoodAndPIPs)

	.field_readonly("totalScore", &BFModelSearch::totalScore)
	.field_readonly("largestLogBF", &BFModelSearch::largestLogBF)
	.field_readonly("nullModelPrior", &BFModelSearch::nullModelPrior)
	.field_readonly("nullModelScore", &BFModelSearch::nullModelScore)
	.field_readonly("logLikelihoodUpToAConstant", \
                         &BFModelSearch::logLikelihoodUpToAConstant)
	.field_readonly("probabilityAtLeast1Causal", \
                         &BFModelSearch::probabilityAtLeast1Causal)
	.field_readonly("globalBayesFactor", &BFModelSearch::globalBayesFactor)
	.field_readonly("PIPs", &BFModelSearch::PIPs)
	.field_readonly("SNPIDs", &BFModelSearch::SNPIDs)

    ;
}
