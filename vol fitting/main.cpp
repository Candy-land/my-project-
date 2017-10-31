#include "_NelderMead.h"

#include <iostream>
#include <vector>

int main() {


	int numOfInstruments = 11;
	double *k = new double[numOfInstruments]; 
	double *pMarket = new double[numOfInstruments]; // numOfInstruments option prices
	
	
	k[0] = 175.0;
	k[1] = 180.0;
	k[2] = 185.0;
	k[3] = 190.0;
	k[4] = 195.0;
	k[5] = 200.0;
	k[6] = 205.0;
	k[7] = 210.0;
	k[8] = 215.0;
	k[9] = 220.0;
	k[10] = 225.0;

	pMarket[0] = 30.8217;
	pMarket[1] = 27.4800;
	pMarket[2] = 24.3799;
	pMarket[3] = 21.5245;
	pMarket[4] = 18.9132;
	pMarket[5] = 16.5415;
	pMarket[6] = 14.4018;
	pMarket[7] = 12.4841;
	pMarket[8] = 10.7759;
	pMarket[9] =  9.2634;
	pMarket[10] = 7.9320;
	
	vector<double> start;
	
	start.push_back( 0.5);
	start.push_back( 0.2);
	start.push_back(-0.2);
	
	objFuncInput input;
	input.numOfInstruments = numOfInstruments;
	input.Strike = k;
	input.pMarket = pMarket;

	input.alpha = 1.5;
	input.eta = 0.25;
	input.N = 10;

	input.s0 = 200;
	input.r = 0.01;
	input.q = 0.015;
	input.t = 0.5;

	input.putCall = 1;

	
	objFunc objf = objFunc();
	optimSet optimSet = {1000, 1000, 1e-08};
	vector<double> prices;
	cout << _NelderMead(start, optimSet, objf, input) << endl;
	
	for (int i=0; i < start.size(); ++i) {
		cout << start[i] << " ";
	}
	cout << endl;
		

	delete[] k;
	delete[] pMarket;
}