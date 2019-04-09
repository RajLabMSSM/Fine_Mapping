#include <iostream>
#include <string>
#include <cassert>
#include <cstdlib>
#include <fstream>
#include <vector>

#include <Eigen/Dense>
using namespace Eigen;

#include "calculateBayesFactors.h"


using namespace std;

double caviarSingleModelAppr(VectorXd z, MatrixXd corMatrix, VectorXd u, \
                             double eps = 1e-3) {
  // This function calculate the log10BF of different models
  // using an approximation that assumes the variance explained by each 
  // model is small

  // input:
  // let m be the number of selected SNPs in the model
  // z: a vector of length m, the marginal t test statistic
  // corMatrix: a matrix of m * m,
  // the correlation matrix of SNPs, use 1/n in the correlation calculation
  // u: a vector of length m,
  // the variance of beta's prior normal distribution multipled by n,
  // n is the number of samples
  // The effect size assumes the genotypes have been normalized,
  // i.e., mean 0 and var 1

  // output:
  // log10(BF)

  int m = z.size();
  // sigmaX = diag(1 / u, m) + corMatrix
  // a small value, 1e-3, is added to ensure stable results when u is large
  double inverseU;
  double inverseUMax = -1;
  for (int i = 0; i < m; i++) {
	  inverseU = 1 / u(i);
	  corMatrix(i, i) += inverseU;
	  if (inverseU > inverseUMax) {
		  inverseUMax = inverseU;
	  }
  }

  if (inverseUMax < eps) {
	  for (int i = 0; i < m; i++) {
		  corMatrix(i, i) += eps;
	  }
  }

  // output sigmaX
  //cout << "sigmaX" << endl << corMatrix << endl;

//  // A = PLU  LU decomposition
//  PartialPivLU<MatrixXd> luResult = corMatrix.lu();
//  MatrixXd luMatrix = luResult.matrixLU();

  // A = LL^T  Cholesky decomposition
  LLT<MatrixXd> lltResult = corMatrix.llt();
  MatrixXd lMatrix = lltResult.matrixL();

//  // A = PTLDL*P Robust Cholesky decomposition
//  LDLT<MatrixXd> ldltResult = corMatrix.ldlt();
//  VectorXd d = ldltResult.vectorD();

  // output LU
  //cout << "LU" << endl << luMatrix << endl;

  //  logBF = -0.5 * (sum(log(diag(abs(L))) + log(diag(abs(R)))) +
  //                  + sum(log(u)))
  //  0.5 * (t(z) %*% solve(sigmaX, z))
  double logDet = 0;
  for (int i = 0; i < m; i++) {
	  // logDet += log(abs(luMatrix(i, i)));  // A = PLU
	  logDet += (2 * log(abs(lMatrix(i, i))));  // A = LL^T
	  // logDet += log(abs(d(i)));  // A = PTLDL*P
	  logDet += log(u(i));
  }
  MatrixXd zTSigmaXz = z.transpose() * lltResult.solve(z);
  double logBF = -0.5 * logDet + 0.5 * zTSigmaXz(0);

  return(logBF);
}

double caviarSingleModelExact(VectorXd& z, MatrixXd& corMatrix, VectorXd& u, \
							     int nSample, double eps = 0) {
	// calculate the eaxct Bayes factor as in BIMBAM
		
	// input:
	// nSample: number of individuals in the data
	// others see caviarSingleModel
	
	// output:
	// see caviarSingleModel
	
	int m = z.size();
	int n = nSample;
		
	double inverseU;
	for (int i = 0; i < m; i++) {
	  inverseU = 1 / u(i);
	  corMatrix(i, i) = 1 + inverseU + eps;  
	}
	
	LLT<MatrixXd> lltResult = corMatrix.llt();
	MatrixXd lMatrix = lltResult.matrixL();
	double logDet = 0;
	for (int i = 0; i < m; i++) {
	  logDet += (2 * log(abs(lMatrix(i, i))));  // A = LL^T
	  logDet += log(u(i));
	}
	MatrixXd zTSigmaXz = z.transpose() * lltResult.solve(z);
	double logBF = -0.5 * logDet -  \
	              0.5 * n * log(1 - zTSigmaXz(0));	
	
	return(logBF);
}

double caviarSingleModelMultiPrior(VectorXd& z, MatrixXd& corMatrix, \
                              VectorXd& u, VectorXd& priorValues, \
							     int nSample, bool exact, double eps) {
   VectorXd logBFs;
   int nPrior = priorValues.size();
   logBFs.resize(nPrior, 1);
   VectorXd u2;
   if (!exact) {
		for (int i = 0; i < priorValues.size(); i++) {
			u2 = u * priorValues(i);
			logBFs(i) =  \
				caviarSingleModelAppr(z, corMatrix, u2, eps);
		}
	} else {
		for (int i = 0; i < priorValues.size(); i++) {
			u2 = u * priorValues(i);
			logBFs(i) = caviarSingleModelExact(z, corMatrix, \
				u2, nSample, eps);
		}
	}
	if (nPrior >= 1) {
		// average all BFs
		double maxLogBFs = logBFs.maxCoeff();
		double meanLogBFs = \
		          log((logBFs.array() - maxLogBFs).exp().sum() / nPrior);
		meanLogBFs += maxLogBFs;
		return meanLogBFs / log(10);
	} else {
		return logBFs(0) / log(10);
	}
}

bool nextCombinationLex(int* t, int k, int n) {
	// generate the next combination in lexicographic order

	// input:
	// n: the total number of objects index by 1 ... n
	// k: the total number of chosen objects
	// t: the selected R object indice

	// output:
	// false: the end of the sequence
	// true: the next combination is in t

    int i = k;  // start from the end
	int ti = 0;
	// find the first place to increase the index
	while( i >= 1 && t[i - 1] == (n - k + i)) {
		i--;
	}
	if (i == 0) { // no next sequence
		return false;
	} else {
		// increase the ith place and set the following indice
		ti = t[i - 1];
		for (int j = i; j <= k; j++) {
			t[j - 1] = ti + 1 + j - i;
		}
		return true;
	}
}

void outputBFAndSNPIDs(VectorXd& BFOutputBufffer, MatrixXi& causalSNPID, \
                       VectorXd& nCausalBuffer, ofstream& outputStream, \
					     int count) {
	int maxCausalSNP = causalSNPID.cols();
	for (int i = 0; i < count; i++) {
		outputStream << fixed << showpos   \
		             << BFOutputBufffer(i) << "\tNA";
		for (int j = 0; j < maxCausalSNP; j++) {
			if (j < nCausalBuffer[i]) {
				outputStream << "\t" << noshowpos << causalSNPID(i, j);
			} else {
				outputStream << "\tNA";
			}
		}
	    outputStream << endl;
	}
}

void caviarBF(VectorXd& z, MatrixXd& corMatrix, vector<double>& weights, \
             bool PVEForPrior, VectorXd& priorValues, int nSample, \
			   int maxCausal, const char* outputFile, bool exact, double eps, \
			    bool useIdentityMatrix) {
  // calcualte the Bayesian factor for all subset of SNPs of size nCausal
  // and write the result in the same format as BIMBAM multiple SNP BF

  // input:
  // maxCausal: maximal causal SNPs
  // exact: whether to use exact BF or approximate BF
  // others see function caviarSingleModel

  unsigned int m = z.size();
  int nOutputBuffer = 1000000;
  VectorXd BFOutputBufffer(nOutputBuffer);
  MatrixXi causalSNPID(nOutputBuffer, maxCausal);
  VectorXd nCausalBuffer(nOutputBuffer);
  int nCausal = 1;
  int count = 0;
  bool nextOK;

  VectorXd u;
  if (!PVEForPrior) {
	  u = MatrixXd::Constant(m, 1, 1 * nSample);
	  priorValues.array() = priorValues.array().square(); // sigmaa^2
  } else {
	  u = MatrixXd::Constant(m, 1, 1 * nSample);
	  // h / (1 - h)
	  priorValues.array() = priorValues.array() / (1 - priorValues.array());
  }
  if (weights.size() > 0) {
	  for (unsigned int i = 0; i < m; i++) {
		  u(i) *= weights[i];
	  }
  }

  // preprocess z and u to obtain exact BF
  if (exact) {
	double scaling1 = (nSample - 1.0) / nSample;
	u *= scaling1;
	ArrayXd temporaryResult;
	temporaryResult = (z.array().square() + nSample - 2).pow(0.5).inverse();
	z.array() *= temporaryResult;
  }
  
  // write the first two lines
  ofstream outputStream(outputFile);
  outputStream << "## note:bf=log10(Bayes Factor), SNP IDs are from 1 ... m" \
               << endl;
  outputStream << "bf\t\tse";
  for (int i = 0; i < maxCausal; i++) {
	  outputStream << "\tsnp" << i + 1;
  }
  outputStream << endl;

  for (nCausal = 1; nCausal <= maxCausal; nCausal++) {
    int* t = new int[nCausal]; // selected SNP IDs in 1 ... m
	// initialize t
	for (int i = 0; i < nCausal; i++) {
		t[i] = i + 1;
	}
	VectorXd zSelected(nCausal);
	MatrixXd corMatrixSelected(nCausal, nCausal);
	corMatrixSelected.setIdentity();
	VectorXd uSelected(nCausal);
	while (true) {
		if (PVEForPrior) {  // sigmaa based on PVE
			if (weights.size() > 0) {
				double sumOfVarianceSelected = 0;
				for (int i = 0; i < nCausal; i++) {
					int ind1 = t[i] - 1;
					sumOfVarianceSelected += weights[ind1];
				}
				for (int i = 0; i < nCausal; i++) {
					int ind1 = t[i] - 1;
					uSelected(i) = u(ind1) / sumOfVarianceSelected;
				}
			} else {
				for (int i = 0; i < nCausal; i++) {
					uSelected(i) = u(0) / nCausal;
				}
			}
		} else { // sigmaa fixed no matter how many causal SNPs in the model
			for (int i = 0; i < nCausal; i++) {
				int ind1 = t[i] - 1;
				uSelected(i) = u(ind1);
			}
		}
		for (int i = 0; i < nCausal; i++) {
			int ind1 = t[i] - 1;
			causalSNPID(count, i) = t[i];
			zSelected(i) = z(ind1);
			if (!useIdentityMatrix) {
  			    for (int j = 0; j < nCausal; j++) {
  				    int ind2 = t[j] - 1;
  				    corMatrixSelected(i, j) = corMatrix(ind1, ind2);
  			    }
			}
		}
		nCausalBuffer[count] = nCausal;
		BFOutputBufffer[count] = \
		       caviarSingleModelMultiPrior(zSelected, corMatrixSelected, \
					uSelected, priorValues, nSample, exact, eps);
		count++;
		if (count == nOutputBuffer) {
			// output BF and SNP IDs
			outputBFAndSNPIDs(BFOutputBufffer, causalSNPID, nCausalBuffer, \
					outputStream, count);
			count = 0;
		}
		nextOK = nextCombinationLex(t, nCausal, m);
		if (!nextOK) {
			break;
		}
	}
	delete [] t;
  }
  outputBFAndSNPIDs(BFOutputBufffer, causalSNPID, nCausalBuffer, \
					outputStream, count);
}

void readZFile(const char* zFilename, VectorXd& z, vector<string>& SNPNames, \
			     vector<double>& weights) {
	// read z file. It is a 2 or 3 columns file. The first is the SNP name,
	// the second is the z score from a single SNP association test, such as
	// a t test statistic. The optional 3rd column is the weight for each SNP

	ifstream infile(zFilename);
	string line;
	vector<double> zValues;
	double zValue;
	string SNPName;
	SNPNames.clear();
	double weight;
	weights.clear();
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		iss >> SNPName >> zValue;
		if (iss >> weight) {
			weights.push_back(weight);
		}
		zValues.push_back(zValue);
		SNPNames.push_back(SNPName);
	}
	int m = zValues.size();
	z.resize(m, 1);
	for (int i = 0; i < m; i++) {
		z(i) = zValues[i];
	}
}

void readCorrelationFile(const char* corFilename, MatrixXd& corMatrix) {
	// read correlation file.
	// It is an m * m matrix written in a file

	ifstream infile(corFilename);
	string line;
	unsigned int m = 0;
	vector<double> corLine;
	double corValue;
	unsigned int nLine = 0;
	while (std::getline(infile, line)) {
		std::istringstream iss(line);
		corLine.clear();
		while(iss >> corValue) {
			corLine.push_back(corValue);
		}
		if (m == 0) {  // reading the first line
			m = corLine.size();
			corMatrix.resize(m, m);
		} else {
			assert(m == corLine.size());
		}
		for (unsigned int i = 0; i < m; i++) {
			corMatrix(nLine, i) = corLine[i];
		}
		nLine++;
	}
	assert(nLine == m);
}

void calculateBayesFactors(string zFilename, \
  						   string corFilename, \
  						   int priorType, \
  						   vector<double> priorValues, \
  						   int nSample, \
  						   int maxCausal, \
  						   string outputFile, \
  						   bool exact, \
  						   double eps, \
						   bool useIdentityMatrix) {
  bool PVEForPrior = false;
  if (priorType == 1) {
	  PVEForPrior = true;
  }
  vector<double> weights;
  VectorXd z;
  MatrixXd corMatrix;
  vector<string> SNPNames;
  Map<VectorXd> priorValuesMapped(priorValues.data(), \
                            priorValues.size());
  VectorXd priorValuesXd(priorValuesMapped);

  readZFile(zFilename.c_str(), z, SNPNames, weights);
  if (!useIdentityMatrix) {
    readCorrelationFile(corFilename.c_str(), corMatrix);
  }
  caviarBF(z, corMatrix, weights, PVEForPrior, \
          priorValuesXd, nSample, \
          maxCausal, outputFile.c_str(), \
		      exact, eps, \
          useIdentityMatrix);
}