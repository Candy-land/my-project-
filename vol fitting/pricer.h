#ifndef __PRICER_H__
#define __PRICER_H__

#include <iostream>
#include <string>
#include <vector>
#include <iterator>
#include <ctime>
#include <map>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <complex>

#include "objFuncInput.h"
//#include "characterFunc.h"


//typedef double(*penaltyFunction)(double, double, double);

#define PI  3.1415926535897932

typedef std::complex<double> dcmplx;
typedef dcmplx(*charFuncOfLogOfStock)(dcmplx&, double, double, double, double, double*);     // pointer to function


dcmplx* FFT_simple(dcmplx*, int);
dcmplx* IFFT_simple(dcmplx*, int);
void FFT_calculate(dcmplx*, int, int, dcmplx*, dcmplx*, dcmplx*);
dcmplx* FFT_get_twiddle_factors(int);

// characteristic related functions
dcmplx charFuncOfLogOfStock_VG(dcmplx&, double, double, double, double, double*);

// pricing function
dcmplx optModifiedCFunc(dcmplx&, double, double, double, double, double, double*, charFuncOfLogOfStock);

double optPriceFFT(double, double, double, double, double, double, double*, int, charFuncOfLogOfStock, double);

//---------------------------------------------------------------
// interpolation
double interpolate(double, double, double, double, double);
void getPrice(double*, double*, std::vector<double>&, std::vector<double>&, int);
//void printprice(double*, double*, double*, int);
//void optPriceFFT(double, double, double, double, double, double, double*, int, charFuncOfLogOfStock, double*, double*);

#endif