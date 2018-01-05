#ifndef __GRADCALC__H
#define __GRADCALC__H

#include <vector>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"

using namespace std;

namespace Professor{class ParamPoints;}

/**
 * This objects calculates gradient vectors of a polynomial function
 */
class GradCalc
{
public:
	//Constructor
	GradCalc(){}

	//Resize @_structure
	void initStructure(const Professor::ParamPoints& pts);
	
	//Getter of the the @i-th gradient vector
	const vector<double> getGradVector(const size_t i, const std::vector<double>& bfp, Professor::ParamPoints& pts, const int order);
	
	//Getter of all gradient vectors
	const vector<vector<double>> getAllGradVectors(const std::vector<double>& bfp, Professor::ParamPoints& pts, const int order);

private:

	//Add terms to @_structure
	void extendStructure(const size_t i, const size_t j, const size_t numFitParams, Professor::ParamPoints& pts, const int order);
	
	//Adds a monomial to @_structure
	const double addMonomial(const size_t i, const size_t j, const size_t k, Professor::ParamPoints& pts, const int order);
	
	//Setter of @_structure
	void setStructure(const size_t i, const size_t j, const vector<double> vec){_structure[i][j] = vec;}

	//List of the fit parameter independent part of the gradient vectors. Structure: [number of anchor points][dimension][monomial]
	vector<vector<vector<double>>> _structure;
};
#endif
