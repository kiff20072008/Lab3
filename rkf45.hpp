#ifndef RKF45_HPP_
#define RKF45_HPP_

#include <math.h>
#include <iostream>
#include <iomanip>


#define  EPSILON  2.2e-16

int rkf45(int(*F)(int n, double t, double y[], double yp[]), int NEQN,
	double Y[], double YP[],
	double *T, double TOUT,
	double *RELERR, double ABSERR,
	double *H,
	int *NFE, int MAXNFE, int *IFLAG);

int fehl45(int(*F)(int n, double t, double y[], double yp[]),
	double T, double H,
	double Y[], double YP[], double F1[],
	double F2[], double F3[], double F4[],
	double F5[], int NEQN);

int rkfinit(int NEQN, int *fail);
int rkfend();

#endif