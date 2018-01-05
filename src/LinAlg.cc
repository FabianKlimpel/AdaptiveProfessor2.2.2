#include "Professor/LinAlg.h"

using namespace std;

/**
 * This function returns a column of a matrix
 * @mat: matrix, out of which the column will be extracted
 * @j: number of column
 * @tmp: column that will be returned
 */
const vector<double> LinAlg::getCol(const vector<vector<double>>& mat, const size_t j){
	vector<double> tmp;
	for(size_t i = 0; i < mat.size(); i++)
		tmp.push_back(mat[i][j]);
	return tmp;
}

/**
 * This function delivers the absolut value of a vector
 * @vec: vector of interest
 * @abs: absolut value
 */
const double LinAlg::getAbs(const vector<double>& vec){
	double abs = 0.;
	//summing up the sqaure of all components of the vector
	for(size_t i = 0; i < vec.size(); i++)
		abs += vec[i] * vec[i];
	//return the sqrt of the sum
	return sqrt(abs);
}

/**
 * This function transposes a given matrix
 * @mat: matrix, that will be transposed
 * @tmp: temporary storage of the transposed matrix
 */
const vector<vector<double>> LinAlg::transpose(const vector<vector<double>>& mat){
	vector<vector<double>> tmp;
	
	//setting @tmp's size to the size of @mat
	tmp.resize(mat[0].size());
	for(size_t i = 0; i < tmp.size(); i++)
		tmp[i].resize(mat.size());
	//transposing by switching the according indices
	for(size_t i = 0; i < mat.size(); i++)
		for(size_t j = 0; j < mat[0].size(); j++)
			tmp[j][i] = mat[i][j];
	return tmp;
}

/**
 * This function multiplies a matrix and a vector
 * @mat: matrix of the product
 * @vec: vector of the product
 * @result: returned vector, represents the product of @mat and @vec
 * @tmp: helper of summing up a row of @mat and the components of @vec
 */
const vector<double> LinAlg::multMatVec(const vector<vector<double>>& mat, const vector<double>& vec){

		vector<double> result;
		
		//setting @result's size to the number of rows in @mat
		result.resize(mat.size());
		double tmp = 0;

		//walking over every row of @mat
		for(size_t i = 0; i < mat.size(); i++)
		{
			//Walking over every column of @mat at a given row and walking over every component of @vec. The product is added to @tmp
			for(size_t j = 0; j < vec.size(); j++)
				tmp += mat[i][j] * vec[j];

			//the sum is assigned to one component of the result and tmp is resetted
			result[i] = tmp;
			tmp = 0;
		}
		return result;
}

/**
 * This function normalizes a vector. If the length of the vector is 0, the norm will be the vector itself in order to prevent nan's.
 * @vec: vector, that will be normalized
 * @abs: absolute value of the vector
 * @result: if the length is != 0, the normalized vector will be stored in this variable
 */
const vector<double> LinAlg::normalizeVec(const vector<double>& vec){

	//calculate the absolute value of @vec
	const double abs = getAbs(vec);
	
	//if the length of @vec is != 0, it will be normalized, else the vector itself will be returned 
	if(abs != 0)
	{
		vector<double> result;
		result.resize(vec.size());

		//calculating the normalized components of the vector
		for(size_t i = 0; i < vec.size(); i++)
			result[i] = vec[i] / abs; 
			
		return result;
	}
	else
		return vec;
}

/**
 * This function calculates the distance between two vectors
 * @a, @b: vectors of interest
 */
const double LinAlg::getDistanceOfVectors(vector<double> a, const vector<double>& b){

	//calculate the difference componentwise
	for(size_t i = 0; i < a.size(); i++)
		a[i] = a[i] - b[i];
		
	//return its absolut value
	return getAbs(a);
}

/**
 * This function delivers the dotproduct of two vectors.
 * The special cases in this function are dotproduct with at least one zerovector. This causes a problem
 * because two vectors are exactly the same if they are zerovectors. This is used as a part of the convergence criteria.
 * In order to avoid this, the dotproduct of two zerovectors is defined as 1.
 * Furthermore, if only one of those vectors is a zerovector, this may causes a similar problem.
 * @a, @b: vectors for the dotproduct
 * @result: result of the dotproduct
 */
const double LinAlg::dotProduct(const vector<double>& a, const vector<double>& b){
	
	//the case of both vectors are zerovectors is handled by explicit definition that the dotproduct of those (normalized) vectors is 1
	if(getAbs(a) == 0. && getAbs(b) == 0.)
		return 1.;
		
	double result = 0;
	
	//summing up the product componentwise
	for(size_t i = 0; i < a.size(); i++)
		result += a[i] * b[i];
	return result;
}


/**
 * This function solves a problem of the type A*x=b and finds the values of the vector x.
 * x is passed to the function since previously set fit parameters are skipped in this function.
 * Every parameter that is left to the estimation of its value is given as a nan.
 * @a: matrix A
 * @x: vector of parameters to fit
 * @b: vector b
 * @result: resulting parameters
 */
const vector<double> LinAlg::getBestFitParameters(const vector<vector<double>>& a, const vector<double>& x, vector<double> b){

	//setting the size of @result
	vector<double> result;
	result.resize(x.size());

	//starting at the bottom line, walking upwards and calculating the parameters
	for(size_t i = x.size(); i > 0; i--)
	{
		//nan is an indicator in a component of @x for a component, that needs to be calculated
		if(isnan(x[i - 1])){
			//calculating a fit parameter
			result[i - 1] = b[i - 1] / a[i - 1][i - 1];
			
			//forward the solution to every row above the regarded one
			for(size_t j = 0; j < i - 1; j++)
				b[j] -= a[j][i - 1] * result[i - 1];
		}
		else{
			//If a values isn't nan, it was set by @LinAlg::collinearity(). It's value will be forwarded as above.
			result[i - 1] = x[i - 1];
			for(size_t j = 0; j < i - 1; j++)
				b[j] -= a[j][i - 1] * result[i - 1];
		}
	}
	return result;
}

/**
 * This function checks, if the columns of the matrix M are collinear. If this is the case, then the diagonal terms of the matrix R would be very small.
 * As a result of this, the corresponding fit parameters would become very big. Numerically, this leads to the situation that the fit function would only 
 * depend on a few fit parameters. Therefore, those parameters are set in this function to avoid this behaiour.
 * @a: Vector that is used to stores set fit parameters
 * @r: Matrix R
 * @b: Vector bprime
 * @threshold: threshold for judging small values of @r on the diagonal
 * @kappa: shift parameter in order to decrease the value of a fit parameter
 * @riihat, @bihat: helper; Used for calculating the value of @a.
 */
void LinAlg::collinearity(vector<double>& a, const vector<vector<double>>& r, const vector<double>& b, const double threshold, const double kappa){
	
	//checking for RR constraints
	for(size_t i = 0; i < a.size(); i++)
	{
		a[i] = nan("1");
		//check if the value in @r is small
		if(fabs(r[i][i]) < threshold)
		{
			//Setting the helpervariables for the creation of a smaller @a than one would get if the fit would run normally.
			//The shift has it's source in @kappa and the according term of @_d.
			const double riihat = sqrt(r[i][i] * r[i][i] + kappa);
			const double bihat = (r[i][i] / riihat) * b[i];
			a[i] = bihat / riihat;
		}
	}
}
