// Multinomial.h
// Kevin Karplus
// 23 September 1995

#ifndef Multinomial_H
#define Multinomial_H

#include "Prob.h"

// compute the number of ways that a sample with
//	sum(counts) elements can have the distribution
//	given by counts.
// (Use Gamma functions to get continuous version of function)

Prob Multinomial(const float *counts, int num_dim);
Prob Multinomial(const float *counts, int num_dim, double SumCounts);

// SumCounts is an optional argument that is the sum of counts---
// already knowing it saves recomputing it

#endif
