#ifndef __QRHANDLER__H
#define __QRHANDLER__H

#include <iostream>
#include <vector>
#include "Professor/LinAlg.h"
#include <limits>
#include "Professor/ParamPoints.h"
#include <cmath>

using namespace std;
namespace Professor{class ParamPoints;};

/**
 * This class handles the matrix M and the vector b of a problem of the type M*x=b by performing a Gram-Schmidt-QR-decomposition.
 */
class QRHandler
{

public:

	//Constructor
	QRHandler(){}

	//Setup and store information
	void init(Professor::ParamPoints& pts, const vector<double>& ptvals);
	
	//Performs an iteration
	void iterate(Professor::ParamPoints& pts, bool walkthrough = false);
	
	//Sets the object to given iteration number
	void load(int order, Professor::ParamPoints& pts, const size_t numFitParams);
	
	//Deletes stored information
	void reset();

	//Getter
	const vector<vector<double>>& getM() const {return _m;}
	const vector<vector<double>>& getR() const {return _r;}
	const vector<double>& getBprime() const {return _bprime;}
	const vector<double>& getB() const {return _b;}
	const int getMaxPower() const {return _max;}
	const size_t getIterationCounter() const {return _iterationcounter;}
	const vector<vector<int>>& getPower() const {return _power;}

private:
	/**
	 * @_m: storage of the matrix M
	 * @_q, @r: storage of the matrices Q & R of the QR-decomposition
	 * @_mprime: columnwise normalized M
	 * @_b: reference values
	 * @_bprime: normalized @_b
	 * @_power: list that contains the powers of the variables at the terms of the fitting function
	 * @_max: maximum of powers
	 * @_iterationcounter: number of iterations in the fitting process
	 */
	vector<vector<double>> _m, _q, _r, _mprime;
	vector<double> _b, _bprime;
	vector<vector<int>> _power;
	int _max; 
	size_t _iterationcounter;

	//Initializer of @_q and @_r
	const bool initQR();
	
	//Performs an iteration of @_m
	void iterateM(Professor::ParamPoints& pts);
	
	//Performs and iteration of @_b
	void iterateb();
	
	//Updates the sizes of @_q and @_r
	void resizeQR();

	//Adds components to @_q and @_r
	void expandQR();

	//Calculates @_mprime from @_m
	void makeMprime();

	//Adds elements to @_m
	void increaseM(Professor::ParamPoints& pts);	
};

#endif

