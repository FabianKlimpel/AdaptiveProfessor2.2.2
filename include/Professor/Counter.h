#ifndef __COUNTER__H
#define __COUNTER__H

#include <vector>

//This code is copied from Professor 2.1.3
//Source: prof2 -> counter.h / counter.cc / Ipol.cc->mkStructure()
using namespace std;

/**
 * This class creates the powers for the parameters in the monomials.
 */
class Counter {
public:

    /**
     * This constructor sets up the member variables
     * @dim: this parameter represents the dimension of the parameters
     * @maxval: this is the order of the monomial
     */
    Counter(const size_t dim, const int maxval) {
      for (unsigned int i = 0; i < dim; ++i) _data.push_back(0);
      _maxval = maxval;
    };
    
    //Destructor
    ~Counter(void);

	//calculating another step in order to get the powers of a monomial
    const bool next(const int index);

	//sums up the powers
    const int sum() const;

	//getter of the powers
    const vector<int> data() const {return _data;}

private:
	//@_maxval stores the information about the order
	//@_data stores the powers
    int _maxval;
    vector<int> _data;
};
#endif
