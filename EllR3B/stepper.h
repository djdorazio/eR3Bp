/*
 *  stepper.h
 *  
 *
 *  Created by Daniel D'Orazio on 6/23/14.
 *  Copyright 2014 Columbia University. All rights reserved.
 *
 */
#include"nr3.h"

struct StepperBase {
	Doub &x;
	Doub xold;
	VecDoub &y, &dydx;
	Doub atol, rtol;
	bool dense;
	Doub hdid;   // actual stepsize used
	Doub hnext; // step size predicted by controller for next step
	double EPS;
	Int n, neqn;
	VecDoub yout, yerr;
	StepperBase(VecDoub_IO &yy, VecDoub_IO &dydxx, Doub &xx, const Doub atoll, const Doub rtoll, bool dens) : x(xx), y(yy), dydx(dydxx), atol(atoll), rtol(rtoll), dense(dens), n(y.size()), neqn(n), yout(n), yerr(n) {}
};
