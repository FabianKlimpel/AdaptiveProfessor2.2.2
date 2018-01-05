#include "Professor/GradCalc.h"

/**
 * This functions sets up @_structure.
 * @pts: Storage of all anchor points
 */
void GradCalc::initStructure(const Professor::ParamPoints& pts){
	if(_structure.empty())
	{
			_structure.resize(pts.numPoints());
			for(size_t j = 0; j < pts.numPoints(); j++)
				_structure[j].resize(pts.dim());
	}	
}

/**
 * This function calculates the fit parameter independent part of a monomial of a gradient vector
 * @i: bin number
 * @j: dimension
 * @k: monomial
 * @pts: storage of the anchor points and powers
 * @order: order of the polynomial function
 * @tmp: storage of the value of the fit parameter independent part of the monomial
 */
const double GradCalc::addMonomial(const size_t i, const size_t j, const size_t k, Professor::ParamPoints& pts, const int order){

	double tmp = 1;
	//extracting the respective parameter values
	for(size_t l = 0; l < pts.dim(); l++)
		//if the parameter is the one that is derived in the component of the gradient, then its derivative is used
		if(j == l)
			//if the value of the parameter is 0 & the power of the parameter to derive is != 1, the whole monomial will be 0 but the special case in the derivative of 0^0 = 1
			//if the power is 0, the whole monomial will be 0 after derived 
			if((pts.pointScaled(i)[l] == 0 && pts.getPower(order)[k][l] != 1) || pts.getPower(order)[k][l] == 0)
			{	
				tmp = 0;
				continue;
			}							
			else			
				//multiply the derived factor
				tmp *= pow(pts.pointScaled(i)[l], pts.getPower(order)[k][l] - 1) * pts.getPower(order)[k][l];	
		else
			//if at least one of the not derived parameters is 0 and its power is !=0, the whole monomial become 0
			if(pts.pointScaled(i)[l] == 0 && pts.getPower(order)[k][l] != 0)
			{
				tmp = 0;
				continue;
			}
			else
				tmp *= pow(pts.pointScaled(i)[l], pts.getPower(order)[k][l]);
	return tmp;
}

/**
 * This function adds new terms to @_structre
 * @i: bin number
 * @j: dimension
 * @numFitParams: number of fit parameters in the polynomial function
 * @pts: storage of the anchor points and powers
 * @order: order of the polynomial function
 * @tmpstruc: temporary storage of the new parts of @_structure
 */
void GradCalc::extendStructure(const size_t i, const size_t j, const size_t numFitParams, Professor::ParamPoints& pts, const int order){

	vector<double> tmpstruc = _structure[i][j];
	tmpstruc.resize(numFitParams);

	//start the calculation after the last calculated component and walk up to the new maximum needed
	for(size_t k = _structure[i][j].size(); k < tmpstruc.size(); k++)
		//put the new monomial part to the @tmpstruc at the specific point in the list
		tmpstruc[k] = addMonomial(i, j, k, pts, order);

	//store the result in @_structure
	setStructure(i, j, tmpstruc);	
}

/**
 * This function serves as a getter of a gradient vector. If the vector is not available yet it will be calculated.
 * @i: bin number
 * @bfp: fit parameters
 * @pts: storage of the anchor points and powers
 * @order: order of the polynomial function
 * @functiongradient: vector that contains the value of the gradient evaluated at the @i-th anchor point
 */
const vector<double> GradCalc::getGradVector(const size_t i, const vector<double>& bfp, Professor::ParamPoints& pts, const int order){
	
	//Initialize
	initStructure(pts);
	vector<double> functiongradient;
	functiongradient.assign(_structure[0].size(), 0);

	//calculating the components of the gradient
	for(size_t j = 0; j < functiongradient.size(); j++)
		//if enough monomials were calculated already, they can be taken directly
		if(_structure[i][j].size() > bfp.size()) 
			//walking over every monomial and calculate the contribution to the gradient vector
			for(size_t k = 1; k < bfp.size(); k++)
				functiongradient[j] += bfp[k] * _structure[i][j][k];
		else
		{
			//extend @_structure and calculate the gradient
			extendStructure(i, j, bfp.size(), pts, order);
			for(size_t k = 1; k < _structure[i][j].size(); k++)
				functiongradient[j] += bfp[k] * _structure[i][j][k];
		}
	//normalizing the gradient
	double tmp = LinAlg::getAbs(functiongradient);
	
	//in order to avoid NaN's, the gradient is only normalized if it isn't a zerovector, else it's normalized vector is the vector itself
	if(tmp != 0)
		for(size_t i = 0; i < functiongradient.size(); i++)
			functiongradient[i] /= tmp;
		
	return functiongradient;
}

/**
 * This function calculates all gradient vectors
 * @bfp: fit parameters
 * @pts: storage of the anchor points and powers
 * @order: order of the polynomial function
 * @result: list of gradient vectors evaluated at each anchor point
 */
const vector<vector<double>> GradCalc::getAllGradVectors(const std::vector<double>& bfp, Professor::ParamPoints& pts, const int order){
	vector<vector<double>> result;
	//calculate each gradient vector
	for(size_t i = 0; i < _structure.size(); i++)
		result.push_back(getGradVector(i, bfp, pts, order));
	return result;	
}
