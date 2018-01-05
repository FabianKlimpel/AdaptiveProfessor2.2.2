#ifndef __FITHANDLER__H
#define __FITHANDLER__H

#include <iostream>
#include <vector>
#include <math.h>
#include "Professor/LinAlg.h"
#include "Professor/OutputHandler.h"
#include "Professor/ParamPoints.h"
#include <Eigen/Dense>
#include "Professor/QRHandler.h"

using namespace std;
using namespace Eigen;

/**
* This class handles the fitting procedure.
*/
class FitHandler
{

public:
	//Default Constructor for fitting
	FitHandler();

	//Constructor used for uncertainty calculation
	FitHandler(const vector<double>& fitparams, const vector<double>& pterrs, Professor::ParamPoints& pts, const int order, const int num_ipol, const string histname);

	//calculates the next iterationstep
	void nextStep(Professor::ParamPoints& pts, const double threshold, const double kappa);

	//calculator of fit uncertainties
	void setFitErrors(const size_t num_ipol, const string histname);

	//getter of the smoothness
	const double getDsmooth(Professor::ParamPoints& pts);

	//getter of Chi2
	const double getChi2() const;

	//init of a new bin
	void startBin(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, const double threshold, const double kappa);

	//getter of the number of fitparameters
	const size_t getNumFitParams() const;

	//setter of a specific bin in a specific iterationstep
	void setToIteration(Professor::ParamPoints& pts, const vector<double>& ptvals, const vector<double>& pterrs, const int num_ipol, const double threshold, const double kappa, const int bestiteration);

	//calculates the dot product between every gradient vector
	const vector<double> getAllGradDotProducts(Professor::ParamPoints& pts);

	//adds 0's as fit parameters and changes the order of the fit parameters in order to become usable in Professor 2.2.1
	void sortFitParams(const vector<int>& match_index);

	//Getter
	const int getMaxPower() const {return _qrh.getMaxPower();}
	const vector<double>& getFitErrors() const {return _bfperr;}
	const vector<vector<int>>& getPowers() const {return _qrh.getPower();}
	const size_t getIterationCounter() const {return _qrh.getIterationCounter();}
	const vector<double>& get_a() const {return _a;}
	const vector<double>& getFitParams() const {return _bfp;}

private:
	/**
	* @_a: vector for the fitparameters; is used as indicator of RR constrained parameters
	* @_d: identity matrix
	* @_bfp, @_bfperr: best fit parameters and corresponding errors
	* @_sigma: error of the reference values
	* @_dsmooth_current: Storage of the dsmooth value in a certain step in order to avoid recalculation. The default value indicates an unset value.
	* @_qrh: Object that handles the linear algebra tasks in the iterations
	* @_gc: Object that handles the gradient calculation
	*/
	vector<double> _a, _d, _bfp, _bfperr, _sigma;
	double _dsmooth_current = 2;
	QRHandler _qrh;
	GradCalc _gc;

	//rescaling @_bfp because of the normalization
	void rescaleBestFitParameters();

	//getter of a single element of the precision matrix
	const double getInvCovMatElement(const size_t i, const size_t j) const;
	
};

#endif

