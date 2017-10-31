#ifndef __OBJFUNC_H__
#define __OBJFUNC_H__
#include "pricer.h"

// using FFT for each opton
// eta = 0.5, N = 2^6

class objFunc {
public:


	objFunc();
	double evalObjFunc(double* v, const objFuncInput& input);
	void getPrice(double* v, const objFuncInput& input, std::vector<double>& price);
	~objFunc();
};

#endif