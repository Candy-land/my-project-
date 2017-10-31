#pragma once
#ifndef __OBJFUNCINPUT_H__
#define __OBJFUNCINPUT_H__

#include <iostream>
#include <vector>

// input for fixed maturity, fixed option type
// we can use class instead of struct as well as no member functions are included
struct objFuncInput {
	int numOfInstruments;
	double *Strike;// = new double[numOfInstruments];       // numOfInstruments strikes	
	double *pMarket;// = new double[numOfInstruments]; // numOfInstruments option prices


					// param for FFT
	double alpha;
	double eta;
	double N;  // 2 to the power of N

			   // param for market
	double s0;
	double r;
	double q;
	double t;
	int putCall;    // 1 for call, -1 for put
};



#endif // !__OBJFUNCINPUT_H__

