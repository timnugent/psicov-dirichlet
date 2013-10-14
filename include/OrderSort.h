// OrderSort.h
// Kevin Karplus
// 9 June 1995
// (Added int, double, and Prob versions 19 Aug 1999)

// ABSTRACT
// OrderSort takes an array of scalars (value) and
//	returns a newly allocated array of integers (index)
//	such that value[index[i+1]] >= value[index[i]]

#ifndef ORDERSORT_H
#define ORDERSORT_H

class Prob;
int * OrderSort(const float * value, int num_to_sort);
int * OrderSort(const int * value, int num_to_sort);
int * OrderSort(const double * value, int num_to_sort);
int * OrderSort(const Prob * value, int num_to_sort);
#endif
