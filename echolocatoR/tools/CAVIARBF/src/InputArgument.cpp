#include <unistd.h>
#include <algorithm>

#include "InputArgument.h"
#include "version.h"

using namespace std;

const string version = CAVIARBF_VERSION;

InputArgumentModelSearch::InputArgumentModelSearch() {
}

InputArgumentModelSearch::~InputArgumentModelSearch()
{
}

void InputArgumentModelSearch::processArguments(int argc, char** argv) {
	try {
		TCLAP::CmdLine cmd("Search models and output probabilities "
		                   "based on Bayes factors", ' ', version);
		TCLAP::ValueArg<std::string> outputFileArg("o","output", \
		                        "output file prefix", \
								    true, "", "string");
		TCLAP::ValueArg<std::string> inputFileArg("i","input", \
								"input file storing Bayes factors", \
								true, "", "string");
		TCLAP::ValueArg<int> nSNPsArg("m","snp-number", \
		                        "the total number of variants in the data", \
								    true, 0, "integer");
		TCLAP::ValueArg<double> priorArg("p","prior", \
		"the prior probability of each SNP being causal, in the range [0, 1)\n"
		 "if it is 0, the prior probability will be set to 1 / m, \n"
		 "where m is the number of variants in the region", \
		 true, 0, "numeric");
		TCLAP::ValueArg<std::string> priorFileArg("f","prior-file", \
			  "the file name specifying the prior probabilities of variants", \
								    true, "", "string");
		TCLAP::SwitchArg outputStepwiseResultArg("s","stepwise", \
			"output stepwise result with rho confidence level", false);
		TCLAP::SwitchArg outputExhaustiveResultArg("e","exhaustive", \
	"output exhaustive search result with rho confidence level", false);
		TCLAP::SwitchArg outputExhaustiveStepwiseResultArg("x","mixed", \
        "output result first using exhaustive and then stepwise search "
		 "with rho confidence level", false);

		cmd.xorAdd(priorArg, priorFileArg);
		cmd.add(outputFileArg);
		cmd.add(outputExhaustiveStepwiseResultArg);
		cmd.add(outputExhaustiveResultArg);
		cmd.add(outputStepwiseResultArg);
		cmd.add(nSNPsArg);
		cmd.add(inputFileArg);

		cmd.parse(argc, argv);

		inputFile = inputFileArg.getValue();
		outputFile = outputFileArg.getValue();
		nSNPs = nSNPsArg.getValue();
		if (priorArg.isSet()) {
			prior = priorArg.getValue();
			priorInAFile = false;
		} else if (priorFileArg.isSet()) {
			priorFile = priorFileArg.getValue();
			priorInAFile = true;
		} else {
			throw BadCommandLine();
		}
		outputStepwiseResult = outputStepwiseResultArg.getValue();
		outputExhaustiveResult = outputExhaustiveResultArg.getValue();
		outputExhaustiveStepwiseResult = \
		            outputExhaustiveStepwiseResultArg.getValue();

	} catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " \
		<< e.argId() << std::endl;
	}
}

InputArgumentCAVIARBF::InputArgumentCAVIARBF() {
}

InputArgumentCAVIARBF::~InputArgumentCAVIARBF()
{
}

void InputArgumentCAVIARBF::processArguments(int argc, char** argv) {
	// Process the command line arguments to assign the following values
	//	  string zFilename;
	//	  string corFilename;
	//	  int priorType;
	//	  vector<double> priorValues;
	//	  int nSample;
	//	  int maxCausal;
	//	  string outputFile;

	try {
		TCLAP::CmdLine cmd("Calculate Bayes factors based on "
							 "summary statistics", ' ', version);
		TCLAP::ValueArg<std::string> zFilenameArg("z","zfile", \
		                        "input file for marginal test statistics", \
								    true, "", "string");
		TCLAP::ValueArg<std::string> corFilenameArg("r","rfile", \
								"input file for correlation matrix", \
								false, "", "string");
		TCLAP::ValueArg<int> priorTypeArg("t","prior-type", \
			"prior type for variant effect size\n"
			"\t0: specify sigmaa\n"
			"\t1: specify the proportion of variance explained (pve)", \
			true, 0, "integer");
		TCLAP::ValueArg<std::string> priorValueArg("a", "prior-values", \
							"the prior values associated with the prior type", \
								    true, "0.1 0.2 0.4", "numeric");	
		TCLAP::ValueArg<double> nSampleArg("n","sample-number", \
		                        "the total number of samples in the data", \
								    true, 0, "numeric");
		TCLAP::ValueArg<double> epsArg("e","epsilon", \
		    "the value added to the diagonal of the correlation matrix. "
			 "The default value is 0 if using exact Bayes factors, "
			 "and 1e-3 if using approximated Bayes factors. "
			 "A recommendded value is 0.1 "
	"when using estimated correlation, e.g., from the 1000 genomes project", \
								    false, 0, "numeric");
		TCLAP::ValueArg<int> maxCausalArg("c","max-causal", \
						"the maximal number of causal variants in the model", \
								    true, 3, "integer");
		TCLAP::ValueArg<std::string> outputFileArg("o","output", \
			  "the output file name for Bayes factors", \
								    true, "", "string");
                TCLAP::SwitchArg useApproximationArg("","appr", \
		   "calculate the approximated Bayes factors instead of the exact", \
		      false);									
		TCLAP::SwitchArg identityCorrelationArg("i","identity", \
                 "use the identity matrix for the correlation matrix"\
                  ", useful for c = 1", false);		
		cmd.add(epsArg);
		cmd.add(useApproximationArg);
		cmd.add(outputFileArg);
		cmd.add(maxCausalArg);
		cmd.add(nSampleArg);
		cmd.add(priorValueArg);
		cmd.add(priorTypeArg);
		cmd.add(corFilenameArg);
		cmd.add(zFilenameArg);
		cmd.add(identityCorrelationArg);

		cmd.parse(argc, argv);

		zFilename = zFilenameArg.getValue();
		corFilename = corFilenameArg.getValue();
		priorType = priorTypeArg.getValue();
		nSample = (int)nSampleArg.getValue();
		maxCausal = maxCausalArg.getValue();
		outputFile = outputFileArg.getValue();
		useIdentityMatrix = identityCorrelationArg.getValue();
		
                exact = !useApproximationArg.getValue();
                
		
		if (!exact && !epsArg.isSet()) {
			eps = 1e-3; // default 1e-3 for approximated BFs
		} else {
			eps = epsArg.getValue();
		}
		
		string priorValueString = priorValueArg.getValue();
		std::replace(priorValueString.begin(), priorValueString.end(), \
		              ',', ' ');
		std::istringstream priorValuesStream(priorValueString);
		double priorValue;
		while(priorValuesStream >> priorValue) {
			priorValues.push_back(priorValue);
		}
		
	} catch (TCLAP::ArgException &e)  {// catch any exceptions
		std::cerr << "error: " << e.error() << " for arg " \
		<< e.argId() << std::endl;
	}

	if (priorType != 0 && priorType != 1) {
		std::cerr << "error: prior-type should be 0 or 1" << std::endl;
		exit(1);
	}

	if (corFilename.length() == 0 && !useIdentityMatrix) {
	  std::cerr << "error: either input a correlation matrix file or " \
	            "use the identity matrix" << std::endl;
	  exit(1);
	}
}
