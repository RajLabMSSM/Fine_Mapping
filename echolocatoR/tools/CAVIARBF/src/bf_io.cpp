#include <fstream>
#include <iostream>
#include <sstream>
#include <ctype.h>
#include <stdlib.h>

#include "bf_io.h"
#include "CombinationTable.h"

using namespace std;

int countField(std::string line) {
	int nField = 0;
	bool whiteSpaceBefore = true;
	for (unsigned int i = 0; i < line.length(); i++) {
		if (!isspace(line[i]) && whiteSpaceBefore) {
			nField++;
			whiteSpaceBefore = false;
		} else if (isspace(line[i])) {
			whiteSpaceBefore = true;
		}
	}
	return nField;
}

BayesFactorData::BayesFactorData(const char* filename, int _nSNP) {
	std::fstream infile(filename);
	int count = 0;
	if (infile.is_open()) {
		std::cout << "reading file: " << filename << std::endl;
		std::string line;
		std::getline(infile, line); // comment line
		std::getline(infile, line); // header line
		// count the number of fields
		int nField = countField(line);

		// allocate memory
		nSNP = _nSNP;
		maxInModel = nField - 2;
		nModel = 0;
		for (int i = 1; i <= maxInModel; i++) {
			nModel += CombinationTable::choose(nSNP, i);
		}
		bayesFactor = new double[nModel];
		nSNPInModel = new int[nModel];
		SNPIDInModel = new int[nModel * maxInModel];

		// readin numbers
		int k = 1; // number of SNPs in the current model
		int nModelPerSNPNumber = CombinationTable::choose(nSNP, k);
		std::string item;
		while (std::getline(infile, line)) {
			count++;
			if (count > nModelPerSNPNumber) {
				k++;
				nModelPerSNPNumber += CombinationTable::choose(nSNP, k);
			}
			nSNPInModel[count - 1] = k;

			std::istringstream iss(line);
			iss >> item;
			bayesFactor[count - 1] = atof(item.c_str());
			iss >> item; // se ignored

			for (int i = 1; i <= k; i++) {
				iss >> item;
				SNPIDInModel[(count - 1) * maxInModel + i - 1] =  \
				                       atoi(item.c_str());
			}
		}
		infile.close();
	} else {
		std::cout << filename << ": ";
		throw FileNotExist();
	}
	if (count != nModel) {
		std::cout << filename << ": ";
		throw LineNumberError();
	}
	std::cout << "finish reading input file" << std::endl;
//	std::cout << "data from the last line" << std::endl;
//	std::cout << bayesFactor[nModel - 1];
//	std::vector<int> SNPID = SNPIDInRow(nModel);
//	for (unsigned int i = 0; i < SNPID.size(); i++) {
//		std::cout << " " << SNPID[i];
//	}
//	std::cout << std::endl;
}

std::vector<int> BayesFactorData::SNPIDInRow(int row) {
	// access SNP IDs in a row

	// input:
	// row: it starts from 1
	std::vector<int> SNPID;
	for (int i = 1; i <= nSNPInModel[row - 1]; i++) {
		SNPID.push_back(SNPIDInModel[(row - 1) * maxInModel + i - 1]);
	}
	return SNPID;
}

int BayesFactorData::SNPID(int row, int col) {
	// access SNP IDs by specifying the row and col

	// input:
	// row: it starts from 1
	// col: it starts from 1
	return(SNPIDInModel[(row - 1) * maxInModel + col - 1]);
}

void readPriorFile(const char* filename, int nSNP, \
					  std::vector<double>& allSNPPriors) {
	std::fstream infile(filename);
	int count = 0;
	allSNPPriors.clear();
	if (infile.is_open()) {
		std::ifstream infile(filename);
		std::string line;
		double prior;
		cout << endl;
		while (std::getline(infile, line)) {
			std::istringstream iss(line);
			iss >> prior;
			// make sure prior is not too small or too large
			if (prior < 0 || prior > 1) {
				throw PriorProbRangeError();
			}
			allSNPPriors.push_back(prior);
			count++;
		}
	} else {
		throw FileNotExist();
	}
	if (nSNP != count) {
		throw LineNumberError();
	}
	std::cout << "finish reading prior file" << std::endl;
}