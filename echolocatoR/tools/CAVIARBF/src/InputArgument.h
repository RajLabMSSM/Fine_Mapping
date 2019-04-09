#ifndef INPUTARGUMENT_H
#define INPUTARGUMENT_H

#include <string>
#include <tclap/CmdLine.h>

using namespace std;

class InputArgumentModelSearch
{
public:
	string inputFile;
	string outputFile;
	int nSNPs;
	double prior;
    bool priorInAFile;
	string priorFile;
	bool outputStepwiseResult;
	bool outputExhaustiveResult;
	bool outputExhaustiveStepwiseResult;
public:
	InputArgumentModelSearch();
	~InputArgumentModelSearch();
	void processArguments(int argc, char** argv);
};

class InputArgumentCAVIARBF
{
public:
  string zFilename;
  string corFilename;
  int priorType;
  vector<double> priorValues;
  int nSample;
  int maxCausal;
  string outputFile;
  bool exact;
  double eps;
  bool useIdentityMatrix;

public:
	InputArgumentCAVIARBF();
	~InputArgumentCAVIARBF();
	void processArguments(int argc, char** argv);
};



struct BadCommandLine:std::exception {
	char const* what() const throw() {
		return "Bad command line arguments";
	}
};

#endif // INPUTARGUMENT_H
