#include "Professor/ParamPoints.h"
#include <algorithm>

namespace Professor {

  using namespace std;


  ParamPoints::ParamPoints(const vector< vector<double> >& ppoints) {
    /// @todo Throw a ParamPointsError or similar rather than this assert
    assert(!ppoints.empty());
    // _parampoints.clear();
    // _locked = false;
    // for (size_t i = 0; i < p.size(); ++i) {
    //   _parampoints.push_back(p[i]);
    // }
    _parampoints = ppoints;
    _locked = true;
    _pow.setDim(dim());
  }


  vector<double> ParamPoints::ptcenters() const {
    vector<double> temp_max, temp_min;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_max.push_back(*max_element(temp.begin(), temp.end()));
      temp_min.push_back(*min_element(temp.begin(), temp.end()));
    }
    vector<double> center;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      center.push_back(temp_min[i] + 0.5* (temp_max[i] - temp_min[i]));
    }

    return center;
  }


  vector<double> ParamPoints::ptmins() const {
    vector<double> temp_min;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_min.push_back(*min_element(temp.begin(), temp.end()));
    }
    return temp_min;
  }


  vector<double> ParamPoints::ptmaxs() const {
    vector<double> temp_max;
    for (size_t i = 0; i < dim(); i++) { // iteration over coordinates
      vector<double> temp;
      for (size_t j = 0; j < numPoints(); j++) { // iteration over anchors
        temp.push_back(_parampoints[j][i]);
      }
      temp_max.push_back(*max_element(temp.begin(), temp.end()));
    }
    return temp_max;
  }


  vector< pair<double, double> > ParamPoints::ptedges() const {
    const vector<double> mins = ptmins();
    const vector<double> maxs = ptmaxs();
    vector< pair<double, double> > edge;
    edge.reserve(dim());
    for (size_t i = 0; i < dim(); i++) {
      edge.push_back( pair<double, double>(mins[i], maxs[i]));
    }
    return edge;
  }

  void ParamPoints::setNames(std::vector<std::string > names) {
    if (_names.size() == 0) { // No names set so far
      if (dim() == names.size()) { // Sanity check
        for (size_t i = 0; i < names.size(); i++) {
          _names.push_back(names[i]);
        }
      }
      else {
        stringstream ss;
        ss << "ParamPoints::setNames: dimension mismatch (" << dim() << "dimensions vs. " << names.size() << " names)  ";
        throw ParamPointsError(ss.str());
      }
    }
    else { // Names already set
      stringstream ss;
      ss << "ParamPoints::setNames: Names already set!";
      throw ParamPointsError(ss.str());
    }
  }

  void ParamPoints::printMeta() const {
    cout << "Nr. of points: " << numPoints() << endl;
    cout << "Dimension:     " << dim() << endl;
    cout << "Center:       ";
    for (size_t i = 0; i < dim(); i++) {
      cout << " " << ptcenters()[i];
    }
    cout << endl;
    cout << "Edges:" << endl;
    const vector< pair<double, double> > edges = ptedges();
    for (size_t i = 0; i < dim(); i++) {
      cout << edges[i].first << " < " << edges[i].second << endl;
    }
    cout << endl;
  }


  void ParamPoints::printPoints() const {
    for (size_t i = 0; i < numPoints(); ++i) {
      cout << "Point " << i << ":" << endl;
      for (size_t j = 0; j < dim(); ++j) {
        cout << _parampoints[i][j] << " ";
      }
      cout << endl;
    }
  }

  /////////////////////////////////////////////////////////////////////////////////////////
	
/**
 * This function maps the parameter points @_parampoints on a [0,1]-hypercube and stores them in @_parampoints_scaled
 * @diff: Temporary storage of the difference between the maximum and minimum in every dimension
 */
void ParamPoints::rescale(){
	
	//if the points were not calculated yet, do it
	if(pointsScaled().empty())
	{
		//resize @_parampoints_scaled so that its values can be set easy
		_parampoints_scaled.resize(numPoints());
		for(size_t i = 0; i < numPoints(); i++)
			_parampoints_scaled[i].resize(dim());
		
		//walk over every dimension and calculate the difference between the maximum and the minimum value
		vector<double> diff;
		for(size_t i = 0; i < dim(); i++)
			diff.push_back(ptmaxs()[i] - ptmins()[i]);
		
		//walk over every parameter set dimension-wise and rescale it on a [0,1]-hypercube
		for(size_t i = 0; i < numPoints(); i++)
			for(size_t j = 0; j < dim(); j++)
				_parampoints_scaled[i][j] = (points()[i][j] - ptmins()[j]) / diff[j];	
	}		
}

/**
 * This function empties all lists in the object
 */
void ParamPoints::clearAll(){
	_parampoints.clear();
	_parampoints_scaled.clear();
	_gradients.clear();
	_names.clear();
	_pow.clearAll();
}

/**
 * This function sets the gradient vectors of each anchor point
 * @ptvals: Bin values
 */
void ParamPoints::setGradients(const vector<double> ptvals){
	GradCalc gc;
	HyperCubeIpol hci;
	//walk over every anchor point and calculate the gradient vector
	for(size_t i = 0; i < numPoints(); i++)
		_gradients.push_back(gc.getGradVector(i, hci.getFitParams(i, points(), ptvals), *this, 1));
}

/**
 * This function calculates the dot product of every pair of gradient vectors of the anchor points
 * @result: storage of the resulting values of the dotproducts
 */
const vector<double> ParamPoints::getAllGradDotProducts() const{
	
	vector<double> result;
	
	//walk over every possible combination of anchor points
	for(size_t i = 0; i < _gradients.size(); i++)
		for(size_t j = 0; j < _gradients.size(); j++)
			//do not calculate the dotproduct if the points are the same
			if(i != j)
				result.push_back(LinAlg::dotProduct(_gradients[i], _gradients[j]));
	return result;
}
}
