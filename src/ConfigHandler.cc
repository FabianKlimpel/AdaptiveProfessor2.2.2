#include "Professor/ConfigHandler.h"

/**
 * Constructor that reads a config file and sets all member variables according to the file
 * @configfile: name of the config file
 * @ifile: input stream for reading the config file
 * @line: string that gets the content of a line in the configfile
 */
ConfigHandler::ConfigHandler(const string configfile){
	
	if(!configfile.empty()){
		cout << "reading config file ...";
		ifstream ifile;
		ifile.open(configfile);
		string line;
		
		if(ifile.is_open())
		{
			//linewise file reading
			while(getline(ifile, line))
			{
				//checking for signal words and calling the respective function
				if(line.substr(0, 12) == "thresholdfit")
					readThresholdFit(line);
					
				if(line.substr(0, 13) == "thresholddata")
					readThresholdData(line);
					
				if(line.substr(0, 12) == "thresholderr")
					readThresholdErr(line);
					
				if(line.substr(0, 8) == "chi2mean")
					readChi2Mean(line);
					
				if(line.substr(0, 5) == "kappa")
					readKappa(line);
					
				if(line.substr(0, 7) == "summary")
					readSummaryFlag(line);
					
				if(line.substr(0, 6) == "outdot")
					readOutDotFlag(line);
					
				if(line.substr(0, 6) == "covmat")
					readCovMat(line);
			}		
		}
		cout << "complete" << endl;
	}
}

/**
 * This function reads the threshold of the RR constrain of the fit
 * @line: contains the threshold value
 */
void ConfigHandler::readThresholdFit(const string line){
		_thresholdfit = atof(line.substr(13).c_str());
}

/**
 * This function reads the threshold of the RR constrain of the hypercube fitting
 * @line: contains the threshold value
 */
void ConfigHandler::readThresholdData(const string line){
		_thresholddata = atof(line.substr(14).c_str());	
}

/**
 * This function reads the threshold of the RR constrain of the error of the fit
 * @line: contains the threshold value
 */
void ConfigHandler::readThresholdErr(const string line){
		_thresholderr = atof(line.substr(13).c_str());	
}

/**
 * This function reads the number of Chi2 values used for the mean calculation
 * @line: contains the number of points
 */
void ConfigHandler::readChi2Mean(const string line){
		_chi2mean = atof(line.substr(9).c_str());	
}

/**
 * This function reads the kappa of the RR constrain
 * @line: contains the value
 */
void ConfigHandler::readKappa(const string line){
		_kappa = atof(line.substr(6).c_str());
}

/**
 * This function reads the flag for writing the summary output
 * @line: contains the value
 */
void ConfigHandler::readSummaryFlag(const string line){
	_summaryflag = atoi(line.substr(8).c_str());	
}

/**
 * This function reads the flag for writing the dot product summary
 * @line: contains the value
 */
void ConfigHandler::readOutDotFlag(const string line){
	_outdotflag = atoi(line.substr(7).c_str());
}

/**
 * This function reads the flag for writing the covariance matrices
 * @line: string that contains the flag
 */
void ConfigHandler::readCovMat(const string line){
	_covmat = atoi(line.substr(7).c_str());
}
