#pragma once

#include <vector>
#include <string>
#include "optimSet.h"
#include "objFuncInput.h"
#include "objFunc.h"

#include <unistd.h>

using namespace std;

#define rho 1.0
#define psi 0.5
#define chi 2.0
#define sigma 0.5


double _NelderMead(vector<double>&start, optimSet& optimSet, objFunc& objf, objFuncInput& input){
	
	int MAX_IT = optimSet.MAX_IT;
	int MAX_FUNC_EVALS = optimSet.MAX_FUNC_EVALS;
	double EPSILON = optimSet.ESPSILON;
	string whichAction = "";
	int nN = start.size();

	double **v = new double*[nN + 1];
	for (int i = 0; i < nN + 1; i++) {
		v[i] = new double[nN];
	}
	double *f = new double[nN + 1];

	double *y = new double[nN];
	double *xbar = new double[nN];
	double *xr = new double[nN];
	double *xe = new double[nN];
	double *xc = new double[nN];
	double *xcc = new double[nN];

	double fxe;     // value of function at expansion point */
	double fxr;     // value of function at shrinkage point */
	double fxc;     // value of function at contraction point */
	double fxcc;    // value of function at reflection point */

	double centroid;
	int vertexMin;             // vertex with smallest value
	int vertexMax;             // vertex with largest value
	int vertexSecondLargest;   // vertex with second largest value


	double dfMax, df;
	double dvMax, dv;

	int i, j;

	int funcEvals;    // track the number of function evaluations */
	int iter;         // track the number of iterations */


					  // create the initial simplex
					  // assume one of the vertices is 0,0

	for (j = 0; j<nN; ++j) {
		
		v[0][j] = start[j];
	}
	double usual_delta = 0.05;             // 5 percent deltas for non-zero terms
	double zero_term_delta = 0.00025;      // Even smaller delta for zero elements of x
	for (i = 0; i<nN; ++i) {

		for (j = 0; j<nN; ++j) {
			y[j] = start[j];
		}

		if (y[i] != 0)
			y[i] = (1 + usual_delta) * y[i];
		else
			y[i] = zero_term_delta;

		for (j = 0; j<nN; ++j) {
			v[i + 1][j] = y[j];
		}
	}
	// find the initial function values
	for (i = 0; i < nN + 1; i++) {
		f[i] = objf.evalObjFunc(v[i], input);
	}

	funcEvals = nN + 1;
	// ==========================================================
	//
	// begin the main loop of the minimization
	//
	// ==========================================================
	for (iter = 1; iter < MAX_IT + 1; ++iter) {

		// find the index of the largest value
		vertexMax = 0;
		for (i = 0; i < 1 + nN; i++) {
			if (f[i] > f[vertexMax]) {
				vertexMax = i;
			}
		}
		// find the index of the smallest value
		vertexMin = 0;
		for (i = 0; i < 1 + nN; i++) {
			if (f[i] < f[vertexMin]) {
				vertexMin = i;
			}
		}
		// find the index of the second largest value
		vertexSecondLargest = vertexMin;
		for (i = 0; i < 1 + nN; i++) {
			if (f[i] > f[vertexSecondLargest] && f[i] < f[vertexMax]) {
				vertexSecondLargest = i;
			}
		}

		// Compute the reflection point
		// calculate the centroid
		for (j = 0; j < nN; ++j) {
			centroid = 0.0;
			for (i = 0; i < 1 + nN; ++i) {
				if (i != vertexMax) {
					centroid += v[i][j];
				}
			}
			xbar[j] = centroid / nN;
		}


		// reflect vertexMax to new vertex xr
		for (j = 0; j < nN; ++j) {
			xr[j] = (1 + rho)*xbar[j] - rho*v[vertexMax][j];
		}
		fxr = objf.evalObjFunc(xr, input);
		funcEvals++;

		if (fxr < f[vertexMin]) {

			// Calculate the expansion point
			for (j = 0; j < nN; ++j) {
				xe[j] = (1 + rho*chi)*xbar[j] - rho*chi*v[vertexMax][j];
			}
			fxe = objf.evalObjFunc(xe, input);

			funcEvals++;

			if (fxe < fxr) {
				whichAction = "expand";
				for (j = 0; j < nN; ++j) {
					v[vertexMax][j] = xe[j];
				}
				f[vertexMax] = fxe;

			}
			else {
				whichAction = "reflect";
				for (j = 0; j < nN; ++j) {
					v[vertexMax][j] = xr[j];
				}
				f[vertexMax] = fxr;
			}
		}
		else { // (fxr >= f[vertexMin])
			if (fxr < f[vertexSecondLargest]) {
				whichAction = "reflect";
				for (j = 0; j < nN; ++j) {
					v[vertexMax][j] = xr[j];
				}
				f[vertexMax] = fxr;
			}
			else { // fxr >= f[vertexSecondLargest]
				   // Perform contraction
				if (fxr < f[vertexMax]) {
					// Perform an outside contraction
					whichAction = "outside contraction";
					for (j = 0; j < nN; ++j) {
						xc[j] = (1 + psi*rho)*xbar[j] - psi*rho*v[vertexMax][j];
					}
					fxc = objf.evalObjFunc(xc, input);
					funcEvals++;
					if (fxc <= fxr) {
						whichAction = "contract outside";
						for (j = 0; j < nN; ++j) {
							v[vertexMax][j] = xc[j];
						}
						f[vertexMax] = fxc;
					}
					else {
						// Perform a shrink
						whichAction = "shrink";
					}
				}
				else {
					// Perform an inside contraction
					for (j = 0; j < nN; ++j) {
						xcc[j] = (1 - psi)*xbar[j] + psi*v[vertexMax][j];
					}
					fxcc = objf.evalObjFunc(xcc, input);
					funcEvals++;
					if (fxcc < f[vertexMax]) {
						whichAction = "contract inside";
						for (j = 0; j < nN; ++j) {
							v[vertexMax][j] = xcc[j];
						}
						f[vertexMax] = fxcc;

					}
					else {
						// Perform a shrink
						whichAction = "shrink";
					}
				}
			}
			if (whichAction == "shrink") {

				for (i = 0; i <= nN; ++i) {
					if (i != vertexMin) {
						for (j = 0; j < nN; ++j) {
							v[i][j] = v[vertexMin][j] + sigma * (v[i][j] - v[vertexMin][j]);
						}
						f[i] = objf.evalObjFunc(v[i], input);
					}
				}
				funcEvals = funcEvals + nN;
			}
		}

	

		// test for function evaluations
		if (funcEvals > MAX_FUNC_EVALS) {
			break;
		}


		// test for convergence
		dfMax = abs(f[vertexMax] - f[vertexMin]);
		dvMax = 0;
		for (i = 0; i<nN; ++i) {
			for (j = 0; j < nN; ++j) {
				dv = abs(v[i][j] - v[i + 1][j]);
				if (dv > dvMax) {
					dvMax = dv;
				}
			}
		}
		if ((dvMax < EPSILON) && (dfMax < EPSILON)) {
			break;
		}

	}
	// end main loop of the minimization

	// find the index of the smallest value
	vertexMin = 0;
	for (i = 0; i <= nN; ++i) {
		if (f[i] < f[vertexMin]) {
			vertexMin = i;
		}
	}
	
	for (j = 0; j<nN; ++j) {
		start[j] = v[vertexMin][j];
	}
	cout << endl;
	funcEvals++;

	// delete  
	for (int i = 0; i < nN + 1; i++) {
		delete[] v[i];
	}
	
	double mae = f[vertexMin];
	delete[] v;

	delete[] f;
	delete[] y;
	delete[] xbar;
	delete[] xr;
	delete[] xe;
	delete[] xc;
	delete[] xcc;
	
	return mae;


}
