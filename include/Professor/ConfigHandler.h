#ifndef __CONFIGHANDLER__H
#define __CONFIGHANDLER__H

#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/**
* This class handles a config file provided by the user.
* It handles the reading of the file and stores all necessary information in the object.
*/
class ConfigHandler
{
public:
	//Constructor
	ConfigHandler(const string configfile);

	//getter of member variables
	const bool getSummaryFlag(){return _summaryflag;}
	const bool getOutDotFlag(){return _outdotflag;}
	const bool getCovmatFlag(){return _covmat;}
	const double getThresholdFit(){return _thresholdfit;}
	const double getThresholdData(){return _thresholddata;}
	const double getThresholdErr(){return _thresholderr;}
	const double getChi2Mean(){return _chi2mean;}
	const double getKappa(){return _kappa;}

private:
	//setter for the member variables
	void readThresholdFit(const string line);
	void readThresholdData(const string line);
	void readThresholdErr(const string line);
	void readChi2Mean(const string line);
	void readKappa(const string line);
	void readSummaryFlag(const string line);
	void readOutDotFlag(const string line);
	void readCovMat(const string line);

	/**
	* @_summaryflag, @_outdotflag: flags for writing additional summaries
	* @_covmat: flag for writing covariance matrices
	* @_thresholdfit: threshold of the RR constrain in the fitting
	* @_thresholddata: threshold of the RR constrain in the hypercube-fitting of the data
	* @_thresholderr: threshold of the RR constrain in getter of the fitting errors
	* @_chi2mean: number of modified chi2's to store in order to state a best value regarding the modified chi2
	* @_kappa: shifting parameter if the RR constrain is applied
	*/
	bool _summaryflag = true, _outdotflag = false, _covmat = true;
	double _thresholdfit = 1e-10, _thresholddata = 1e-10, _thresholderr = 1e-10, _chi2mean = 100, _kappa = 1e-10;
};

#endif
