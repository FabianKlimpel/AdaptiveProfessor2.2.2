#ifndef __HYPERCUBEIPOL__H
#define __HYPERCUBEIPOL__H

#include <vector>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"
#include "Professor/QRHandler.h"

using namespace std;

/**
 * This class handles the selection of anchor points for a hypercube around a given point and calculates a polynomial fit of 1st order through all points.
 */
class HyperCubeIpol
{
public:
	//Constructor
	HyperCubeIpol();

	//Getter of the indices referring to a hypercube surrounding a point with index @i
	const vector<size_t>& gethypercube(const size_t i, const vector<vector<double>>& pts);

	//Getter of @_hypercubes
	const vector<vector<size_t>>& hypercubes() {return _hypercubes;}
	
	//Getter of the fit parameters of all anchor points
	const vector<vector<double>>& getAllFitParams(const vector<vector<double>>& pts, const vector<double>& ptvals);
	
	//Getter of the fit parameter for the @i-th anchor point
	const vector<double>& getFitParams(const size_t i, const vector<vector<double>>& pts, const vector<double>& ptvals);
	
private:

	//Construction of a hypercube around the @i-th anchor point
	void buildHyperCube(const size_t i, const vector<vector<double>>& pts);
	
	//Calculator of the fit parameters for the @i-th anchor point
	void calcFitParams(const size_t i, const vector<vector<double>>& pts, const vector<double>& ptvals);

	//list of indices for hypercubes of the anchor points
	vector<vector<size_t>> _hypercubes;
	
	//list of fit parameters of the anchro points
	vector<vector<double>> _fitparams;
};

#endif

