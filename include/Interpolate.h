// Interpolate.h
// Thu Dec 20 05:54:36 PST 2007 Kevin Karplus

#ifndef INTERPOLATE_H
#define INTERPOLATE_H

using namespace std;
#include <vector>



// Do linear interpolation based on a table
// The table consists of two vectors, x_tab and y_tab, of the same
//	length, sorted in increasing order of x.
// If x is outside the range of x_tab, then extrapolate.
// Clipping can be done instead of extrapolation by duplicating the
//	first or last point of the table
double interpolate(double x, 
	const vector<double> &x_tab,
	const vector<double> &y_tab);
	
// Because vectors are a bit of a pain to initialize, 
//	we also provide an implementation using ordinary arrays.
double interpolate(double x, 
	unsigned int size,
	const double * x_tab,
	const double * y_tab);
	


#endif
