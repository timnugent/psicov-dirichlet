// LogGamma.h
// 31 Oct 1995
// Kevin Karplus

// Compute log Gamma(x), using the math library lgamma function,
// but cache the results to avoid unnecessary recomputation



#ifndef LogGamma_H
#define LogGamma_H

#include <iostream>

using namespace std;

extern double LogGamma(double x);

extern double LogGamma_1(double x);	// first derivative
inline double Psi(double x) {return LogGamma_1(x);}

// log Gamma and its first two derivatives.
extern void LogGamma_derivs(double x, double &log_gamma, 
		double &log_gamma_1, double &log_gamma_2);

extern void LogGamma_print_summary(ostream &out);
	// print some statistics about the LogGamma cache
#endif
