#ifndef __LINALG__H
#define __LINALG__H

#include <vector>
#include <math.h>

using namespace std;

/**
* This Class is a container of functions related to linear algebra.
*/
class LinAlg
{
public:

	//extracts a column from a matrix
	static const vector<double> getCol(const vector<vector<double>>& mat, const size_t j);

	//calculates the absolut value of a vector
	static const double getAbs(const vector<double>& vec);

	//transposes a matrix
	static const vector<vector<double>> transpose(const vector<vector<double>>& mat);

	//multiplies a matrix and a vector
	static const vector<double> multMatVec(const vector<vector<double>>& mat, const vector<double>& vec);

	//normalizes a vector
	static const vector<double> normalizeVec(const vector<double>& vec);

	//calculates the difference between two vectors
	static const double getDistanceOfVectors(vector<double> a, const vector<double>& b);

	//calculates the dotproduct of two vectors
	static const double dotProduct(const vector<double>& a, const vector<double>& b);

	//solves a problem of the type matrix * vector = vector
	static const vector<double> getBestFitParameters(const vector<vector<double>>& a, const vector<double>& x, vector<double> b);

	//checks if diagonal terms of @r are too small and regulates the resulting fit
	static void collinearity(vector<double>& a, const vector<vector<double>>& r, const vector<double>& b, const double threshold, const double kappa);
	
private:

	//Constructor
	LinAlg(){}
};

#endif
