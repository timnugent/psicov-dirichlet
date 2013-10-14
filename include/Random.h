// Random.h
//	utility routines for interface to random()
//	By using only these interface routines, changing random number generators
//	should require changing only this header file

#ifndef RANDOM_H
#define RANDOM_H

using namespace std;
#include <vector>
#include <stdlib.h>	// for random()
#include <assert.h>


// return value uniformly in (low,high)
inline double drandom(double low=0., double high=1.)
{   // return (random()+0.5)/(RAND_MAX+1.0) * (high-low) + low;
    // return drand48()*(high-low)+low;
    return (rand()+0.5)/(RAND_MAX+1.0) * (high-low) + low;
}

// return value uniformly in [0..modulus-1]
inline unsigned long int irandom(unsigned long int modulus)
{//    unsigned long int too_big = (RAND_MAX/modulus)*modulus;
//    assert(too_big>0);
//    long int r = random();
//    while (r>=too_big) r=random();
//    return r%modulus;
    return static_cast<unsigned long int> (drandom()*modulus);
}

// return value uniformly in [low..high]
inline long int irandom(long int  low, long int high)
{    return low+irandom(static_cast<unsigned long int>(high-low+1));
}

inline void set_random(unsigned long int seed)
{   // srandom(seed);
    // srand48(seed);
    srand(seed);
}

inline int get_random_given_cum_weights(int numr, const double *cum_weight)
{
    if (numr<=1) return 0;
    double w= cum_weight[numr-1];
    if (w<=0) return irandom(numr);
    double x = drandom()*w;
    
    int clo=-1, chi=numr-1;
    int cmid;
    // do binary search to determine component
    while (clo<chi-1)
    {   // invariant:     cum_weight[clo] < x <=  cum_weight[chi]
	assert (x <= cum_weight[chi]);
	cmid = (clo+chi+1)/2;
	if (x > cum_weight[cmid]) clo=cmid;
	else chi=cmid;
    }
    return chi;
}

inline int get_random_given_weights(int numr, const double *weight)
{
    if (numr<=1) return 0;
    double *cum_weight= new double[numr];
    double sum=0;
    int r;
    for (r=0; r<numr; r++)
    {	double w = weight[r];
        if (w<0) w=0;
        sum += w;
	cum_weight[r] = sum;
    }
    r = sum>0? get_random_given_cum_weights(numr, cum_weight)
    	: irandom(numr);
    delete [] cum_weight;
    return r;
}

inline int get_random_given_weights(const vector<double> &weight)
{
    int numr = weight.size();
    if (numr<=1) return 0;
    double *cum_weight= new double[numr];
    double sum=0;
    int r;
    for (r=0; r<numr; r++)
    {	double w = weight[r];
        if (w<0) w=0;
        sum += w;
	cum_weight[r] = sum;
    }
    r = sum>0? get_random_given_cum_weights(numr, cum_weight)
    	: irandom(numr);
    delete [] cum_weight;
    return r;
}


inline int get_random_given_weights(int numr, const float *weight)
{
    if (numr<=1) return 0;
    double *cum_weight= new double[numr];
    double sum=0;
    int r;
    for (r=0; r<numr; r++)
    {	double w = weight[r];
        if (w<0) w=0;
        sum += w;
	cum_weight[r] = sum;
    }
    r = sum>0? get_random_given_cum_weights(numr, cum_weight)
    	: irandom(numr);
    delete [] cum_weight;
    return r;
}

// return a random name of the specified length
inline char* random_name(int length=10)
{
    const char vowel[6] = "aeiou";
    const char consonant[19] = "bcdfghjklmnprstvwz";
    char *word = new char[length+1];
    for (int i=0; i<length; i++)
    {   word[i] = i%2?	vowel[irandom(5)]
    		    : consonant[irandom(18)];
    }
    word[length] = 0;
    return word;
}


// CHANGE LOG:
// 4 August 2003 Kevin Karplus
//	added get_random_given_cum_weights and get_random_given_weights
// 27 Feb 2004 Kevin Karplus
//	Converted old-style casts.
// 9 April 2004 Kevin Karplus
//	merged various Random.h files into Utilities/
// Wed Feb  9 09:30:27 PST 2005 Kevin Karplus
//	Made get_random_given_weghts and get_random_given_cum_weights 
//	a little safer (handling all-zero weights, by dropping back to uniform)
// Mon Jun 13 04:08:23 PDT 2005 Kevin Karplus
//	Made  get_random_given_weights and get_random_given_cum_weights 
//	take const arguments
// Mon Jun 13 04:43:48 PDT 2005 Kevin Karplus
//	Added get_random_given_weights with const float* argument.
// Wed Jul 20 16:28:30 PDT 2005 Kevin Karplus
//	Changed get_random_given_cum_weights to do binary search
// Sat Dec  8 14:13:32 PST 2007 Kevin Karplus
//	Added random_name()
// Sun Dec  9 15:45:29 PST 2007 Kevin Karplus
//	Added vector version of get_random_given_weights for double
#endif
