#include "objFunc.h"

using namespace std;

objFunc::objFunc() {}
objFunc::~objFunc(){}

double objFunc::evalObjFunc(double *v, const objFuncInput& input) {
	double alpha = input.alpha;
	double eta = input.eta;
	int pow2 = input.N;
	double N = pow(2, pow2);

	double s0 = input.s0;
	double r = input.r;
	double q = input.q;
	double t = input.t;

	int putCall = input.putCall;

	double* param = v;
		
	int numOfInstruments = input.numOfInstruments;

	double *Strike  = input.Strike;
	double *pMarket = input.pMarket;// new double[numOfInstruments]; // numOfInstruments option prices
	
	double pModel;

	double wi;
	
	double rms = 0.0;
	for (int i=0; i < numOfInstruments; ++i){
		pModel = optPriceFFT(putCall*alpha, s0, r, q, t, eta, v, N, &charFuncOfLogOfStock_VG, Strike[i]);
		//cout << putCall << " " << alpha << " " << s0 << " " << r << " " <<  q << " " <<  t << " " << Strike[i] << " " << N << " " << pModel << " " << v << " " << charFuncOfLogOfStock_VG << endl;
		if (i == 0){
			wi = 1.0;
		}
		else{
			wi = 1.0;
		}
		rms += wi*(pModel - pMarket[i])*(pModel - pMarket[i]);
		//rms += fabs(pModel - pMarket[i]);
	}

	rms = sqrt(rms);

	std::cout << "mae: " << rms << std::endl;

	return rms;

}


void objFunc::getPrice(double* v, const objFuncInput& input, std::vector<double>& price) {
	double alpha = input.alpha;
	double eta = input.eta;
	int pow2 = input.N;
	double N = pow(2, pow2);

	double s0 = input.s0;
	double r  = input.r;
	double q  = input.q;
	double t  = input.t;

	int putCall = input.putCall;


	double rms = 0.0;
	
	//cout << "getPrice " << v[0] << " " << v[1] << " " << v[2] << endl;
	double* param = v;

	int numOfInstruments = input.numOfInstruments;

	double *Strike = input.Strike;

	for (int i = 0; i < numOfInstruments; i++) {
		double pModel = optPriceFFT(putCall*alpha, s0, r, q, t, eta, v, N, &charFuncOfLogOfStock_VG, Strike[i]);
		price.push_back(pModel);
	}



}