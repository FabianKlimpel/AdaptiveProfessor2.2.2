#ifndef __OUTPUTHANDLER__H
#define __OUTPUTHANDLER__H

#include <iostream>
#include <vector>
#include <fstream>
#include "Professor/LinAlg.h"
#include "Professor/ParamPoints.h"
#include "Professor/FitHandler.h"
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

class FitHandler;

/**
* This Class is a container for output functions
*/
class OutputHandler
{
public:

	//Default constructor 
	OutputHandler();
	
	//Constructor that sets up some necessary parameters
	OutputHandler(const Professor::ParamPoints& pts, const bool outdotflag, const bool summaryflag, const size_t& num_ipol);
	
	//writes a summary of the fit to the terminal
	void writeBinResult(const size_t num_ipol, Professor::ParamPoints& pts, FitHandler& fh) const;

	//writes a dotproduct-summary
	void writeDotProduct(const size_t num_ipol, Professor::ParamPoints& pts, FitHandler& fh) const;

	//writes to the summaryfile
	void writeSummary(FitHandler& fh, Professor::ParamPoints& pts) const;
	
	//setup of the summary file
	void setupSummary() const;

	//writes a covariance matrix to file
	void writeCovMat(const MatrixXd& mat, const size_t num_ipol, const string histname) const;

private:

	/**
	 * @distances: distances of every anchor point to another
	 */
	vector<double> _distances;
	
	//calculates the distances between the anchor points
	void setDistances(const Professor::ParamPoints& pts);
};

#endif
