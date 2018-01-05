#include "Professor/Power.h"

/**
 * This function calculates the list of powers up to a given order
 * and adds them to the overall list @_powerlist
 * @order: order up to which the powers will be calculated
 * @tmp_power: temporary storage for the new list of powers until it will be added to @_powerlist
 * @zero: list of 0's as the 0th order of powers
 * @c: object that delivers the powers of the monomials
 */
void Power::setPowerOfOrder(const int order){
	vector<vector<int>> tmp_power;

	//setting the 0th order
	const vector<int> zero(n, 0);
	tmp_power.resize(1);
	tmp_power[0] = zero;

	//loop up to the maximum order
	for (int i = 0; i <= order; ++i) 
	{
		//create a @Counter object of the order @i
		Counter c(n, i);
		
		//calculate every combination of powers
		//if they fit the order @i, they will be stored in @tmp_power
		while (c.next(n-1)) 
			if (c.sum() == i) 
				tmp_power.push_back(c.data());     
	}		
	//saving the calculated list in @_powerlist at the position of the order		
	_powerlist[order] = tmp_power;
}

/**
 * This function is a getter for the list of powers of a given order.
 * If the order isn't calculated yet, it will be calculated in this function.
 * @order: order of the polynomial
 */
const vector<vector<int>> Power::getPowerOfOrder(const int order){

	//the length of @_powerlist is equal to the highest order calculated
	//if the requested order is higher, the size of @_powerlist will be adjusted
	if(order > ((int) _powerlist.size()) - 1)
		_powerlist.resize(order + 1);
	
	//if the order is not calculated yet, it will be calculated
	if(_powerlist[order].empty())
		setPowerOfOrder(order);

	return _powerlist[order];
}
