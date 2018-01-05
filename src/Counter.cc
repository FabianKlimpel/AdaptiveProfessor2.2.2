#include "Professor/Counter.h"

//Destructor
Counter::~Counter(void) { }

/**
 * This function calculates a new step in order to calculate the powers for a monomial.
 * The idea is a list of powers, starting at [0, 0, ... ,0], that will be increased in every step
 * until it fits the needed order. Starting with the last index, the value of an index gets increased by 1.
 * If the index is at the maximum, it will be resetted to 0 and the next index will be increased
 * in a recursive way.
 * That way, every possible combination of powers in a certain order can be calculated.
 * @index: this is the parameter, that will be modified
 */
const bool Counter::next(const int index) {
	//if the index is at the maximum...
	if (_data[index] == _maxval) {
		//if the first index is at his maximum, the function returns false in order to break an outer for-loop
		//this condition means that every combination was created
		if(index == 0) return false;
		//...the index will be resetted and the next index will be calculated
		_data[index] = 0;
		return next(index - 1);
	}
	else {
		//If the entry is smaller than the maximum, it will be increased.
		//the return value is meant as continuitation of an outer for-loop
		_data[index]++;
		return true;
	}
}

/**
 * This function calculates the sum of all powers in @_data
 * @sum_v: storage for the sum
 */
const int Counter::sum() const{
	int sum_v = 0;
	for (size_t n = 0; n < _data.size(); n++)
	  sum_v += _data[n];
	return sum_v;
}
