#ifndef BF_IO_H
#define BF_IO_H

#include <string>
#include <vector>
#include <exception>

struct BayesFactorData {
	int nModel;  // total number of models
	int nSNP;   // total number of SNPs analyzed
	int maxInModel; // maxmimal number of SNPs in a model
	double* bayesFactor; // an array of calculated Bayes factors
	int* SNPIDInModel; // a two dimensional matrix for SNP IDs in each model
	int* nSNPInModel; // an array recording the number of SNPs in each model

	BayesFactorData(const char* filename, int _nSNPs);
	~BayesFactorData() {
		delete [] bayesFactor;
		delete [] SNPIDInModel;
		delete [] nSNPInModel;
	}

	std::vector<int> SNPIDInRow(int row);
	int SNPID(int row, int col);
};

struct ReadError:public std::exception {
	const char* what() const throw () {
		return "Reading input file error";
	}
};

struct FileNotExist:public ReadError {
	const char* what() const throw () {
		return "File does not exist";
	}
};

struct LineNumberError:public ReadError {
	const char* what() const throw () {
		return "line number is not correct";
	}
};

struct PriorProbRangeError:public ReadError {
	const char* what() const throw () {
		return "some prior probabilities are negative";
	}
};

void readPriorFile(const char* filename, int nSNP, \
					 std::vector<double>& allSNPPriors);

#endif


