#include "Professor/OutputHandler.h"

/**
 * Constructor
 */
OutputHandler::OutputHandler(){}

/**
 * Constructor
 * Here, the sizes of the vectors that will later contain fitparameters etc. will be set. Additionally, the distances and centers will be calculated
 * @pts: Storage of the anchor points
 * @outdotflag: Flag for writing dotproducts to file
 * @summaryflag: Flag for writing additional summaries to file
 * @num_ipol: Bin number
 */
OutputHandler::OutputHandler(const Professor::ParamPoints& pts, const bool outdotflag, const bool summaryflag, const size_t& num_ipol){	
	//if the dotproduct-summary should be written, the distances are calculated
	if(outdotflag)
		setDistances(pts);
	//if the summary should be written, the file will be created
	if(summaryflag && (num_ipol == 0))
		setupSummary();
}	

/**
 * This function calculates the distances between the anchors points
 * @pts: Storage of the anchor points
 * @skips: No distances between an anchor point and itself is calculated, so the index shift needs is stored in this variable
 */
void OutputHandler::setDistances(const Professor::ParamPoints& pts){

	//Resizing the @_distances vectors. No distances between an anchor point and itself will be calculated, therefore the size is modified accordingly.
	_distances.resize(pts.numPoints() * pts.numPoints() - pts.numPoints());
	
	for(size_t i = 0; i < pts.numPoints(); i++)
	{		
		size_t skips = 0;
		//compare every point with every other
		for(size_t j = 0; j < pts.numPoints(); j++)
			//if both points are the same, the iteration is skipped
			if(i == j)
				skips--;				
			else
				//calculate the distance between the vectors and store them
				_distances[i * pts.numPoints() - i - 1 + j + skips] = LinAlg::getDistanceOfVectors(pts.pointScaled(i), pts.pointScaled(j));
	}
}

/**
 * This function creates the summary file and writes its header
 * @outsummary: Provides the output to file functionality
 */
void OutputHandler::setupSummary() const{	
	//set up the file and write the header
	ofstream outsummary;
	outsummary.open("summary", ofstream::out | ofstream::app);
	outsummary << "Chi^2\tChi2^2,red\tIterations\tDsmooth" << endl;
	outsummary.close();
}

/**
 * This function writes the result of the fit of a bin to the terminal
 * @num_ipol: Bin number
 * @pts: Storage of the anchor points
 * @fh: Handler of the fit
 */
void OutputHandler::writeBinResult(const size_t num_ipol, Professor::ParamPoints& pts, FitHandler& fh) const{
	cout << endl << "Result for bin " << num_ipol << ":" << endl;

	//the indices of the constrained monomials are printed to the terminal
	cout << "RR constraint:\t";
	for(size_t i = 0; i < fh.getNumFitParams(); i++)
		if(!std::isnan(fh.get_a()[i]))
			cout << i << "\t";
	cout << endl;
	
	//write further summary variables
	cout << "Dsmooth:\t\t" << fh.getDsmooth(pts) << endl;
	cout << "chi2:\t\t\t" << fh.getChi2() << endl;
	cout << "iterationcounter:\t" << fh.getIterationCounter() << endl;
	cout << "max. power:\t\t" << fh.getMaxPower() << endl;
	cout << "-------------------------" << endl;	
}

/**
 * This function writes the dotproduct-summary
 * @num_ipol: Bin number
 * @pts: Storage of the anchor points
 * @fh: Handles the fit
 * @outdot: Provides the output to file functionality
 * @rhdot, @fhdot: Stores the dot products of the reference data and the fit at every anchor point
 */
void OutputHandler::writeDotProduct(const size_t num_ipol, Professor::ParamPoints& pts, FitHandler& fh) const{
	//set up the output
	ofstream outdot;
	outdot.open(("dotproduct" + to_string(num_ipol)).c_str());
	
	//calculate all dot products
	vector<double> rhdot = pts.getAllGradDotProducts(), fhdot = fh.getAllGradDotProducts(pts);
	
	//write some summary parameters
	outdot << fh.getDsmooth(pts) << "\t" << fh.getChi2() << "\t" << fh.getIterationCounter() << endl;
		
	//write the dot products of the reference data, the fit and the distances between the anchor points used for the dot product
	for(size_t k = 0; k < rhdot.size(); k++)
		outdot << rhdot[k] << "\t" << fhdot[k] << "\t" << _distances[k] << "\t";
	outdot << endl;
						
	outdot.close();
}
	
/**
 * This function writes a summary of a fit to the summary file
 * @fh: Handles the fit
 * @pts: Storage of the anchor points
 * @outsummary: Provides the output to file functionality
 */
void OutputHandler::writeSummary(FitHandler& fh, Professor::ParamPoints& pts) const{
	//open the file and continue writing at the end
	ofstream outsummary;
	outsummary.open("summary", ofstream::out | ofstream::app);
	//write out the resulting Chi^2, the number of iterations needed and the smoothness
	outsummary << fh.getChi2() << "\t" << fh.getIterationCounter() << "\t" << fh.getDsmooth(pts) << endl;
	outsummary.close();
}

/**
 * This function writes the covariance matrix of the fit parameters to file.
 * @mat: Covariance matrix
 * @num_ipol: Bin number
 * @outcovmat: Provides the output to file functionality 
 */
void OutputHandler::writeCovMat(const MatrixXd& mat, const size_t num_ipol, const string histname) const{
	//open the file
	ofstream outcovmat;
	if(histname[0] == '/')
		outcovmat.open((histname.substr(1, histname.find("/", 1) - 1) + "_" + histname.substr(histname.find("/", 1) + 1) + "_" + to_string(num_ipol)).c_str());
	else
		outcovmat.open((histname + to_string(num_ipol)).c_str());
	
	//write the matrix
	for(size_t row = 0; row < (size_t) mat.rows(); row++)
	{
		for(size_t col = 0; col < (size_t) mat.cols(); col++)
			//if a matrix element would become +/- inf, the value is replaced by +/- maximum of double
			if(mat(row, col) == std::numeric_limits<double>::infinity())
				outcovmat << std::numeric_limits<double>::max() << " ";
			else
				if(-mat(row, col) == std::numeric_limits<double>::infinity())
					outcovmat << -std::numeric_limits<double>::max() << " ";
				else
					outcovmat << mat(row, col) << " ";
		outcovmat << "\n";
	}
	outcovmat.close();
}
