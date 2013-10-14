// LogGamma.cc
// 31 Oct 1995
// Kevin Karplus

// Compute log Gamma(x), using the math library lgamma function,
// but cache the results to avoid unnecessary recomputation
// Only intended for x>0.

// LogGamma is slightly faster than lgamma, particularly if repeated values
//	close to each other are used.
// LogGamma_derivs is about four times? as fast as using 3 calls to lgamma to
//	compute the first two derivatives.

// To Do:
//	check to see if Stirling approx good enough for LogGamma,
//	as well as derivatives.
//
//	Maybe use an extra term of Stirling approx?

#include "LogGamma.h"
#include <math.h>
#include <iostream>
#include "minmax.h"
#include <assert.h>

static const double epsilon= 0.005;	// tolerable error in LogGamma
//	Note: the actual error is far less, since this is
//	the error one would get by moving the full interval length
//	away from the point at the highest derivative in the cache.
//	Since we actually interpolate, the expected max error is
//	far less.

//	Worst observed for epsilon=0.005, CacheFrom=2.0, CacheUpTo=128
// x= 2.00233 lgamma_0(x)= 0.000987688 LogGamma_0(x)= 0.000987774 err= 8.55934e-08
// x= 2.01006 lgamma_1(x)= 0.429249 LogGamma_1(x)= 0.429249 err= 1.40765e-07
// x= 9.71096 lgamma_2(x)= 0.108441 LogGamma_2(x)= 0.10846 err= 1.91138e-05

//	Worst observed for epsilon=0.005, CacheFrom=0.5, CacheUpTo=128
// x= 0.500715 lgamma_0(x)= 0.570962 LogGamma_0(x)= 0.570962 err= 6.54432e-07
// x= 1.46345 lgamma_1(x)= 0.00175763 LogGamma_1(x)= 0.00175794 err= 3.09863e-07
// x= 9.71096 lgamma_2(x)= 0.108441 LogGamma_2(x)= 0.10846 err= 1.91138e-05

//	Worst observed for epsilon=0.005, CacheFrom=0.8, CacheUpTo=128
// x= 0.800827 lgamma_0(x)= 0.151263 LogGamma_0(x)= 0.151263 err= 3.05239e-07
// x= 1.46345 lgamma_1(x)= 0.00175763 LogGamma_1(x)= 0.00175794 err= 3.09863e-07
// x= 9.71096 lgamma_2(x)= 0.108441 LogGamma_2(x)= 0.10846 err= 1.91138e-05


//	Worst observed for epsilon=0.005, CacheFrom=1.0, CacheUpTo=128
// x= 1 lgamma_0(x)= 0 LogGamma_0(x)= 2.05276e-07 err= 2.05276e-07
// x= 1.46345 lgamma_1(x)= 0.00175292 LogGamma_1(x)= 0.00175322 err= 3.05435e-07
// x= 9.71096 lgamma_2(x)= 0.108441 LogGamma_2(x)= 0.10846 err= 1.91138e-05

//	Worst observed for epsilon=0.01, CacheFrom=1.0, CacheUpTo=128
// x= 0.329383 lgamma_0(x)= 0.997873 LogGamma_0(x)= 0.997873 err= 5.85024e-07
// x= 1.47686 lgamma_1(x)= 0.014634 LogGamma_1(x)= 0.0146352 err= 1.22461e-06
// x= 8.05173 lgamma_2(x)= 0.132245 LogGamma_2(x)= 0.132228 err= 1.72378e-05

//	Worst observed for epsilon=0.001, CacheFrom=1.0, CacheUpTo=32
// x= 0.32912 lgamma_0(x)= 0.998707 LogGamma_0(x)= 0.998707 err= 1.15549e-08
// x= 32.0046 lgamma_1(x)= 3.45018 LogGamma_1(x)= 3.45018 err= 5.89796e-08
// x= 8.05173 lgamma_2(x)= 0.132245 LogGamma_2(x)= 0.132228 err= 1.72439e-05

//	Worst observed for epsilon=0.0005, CacheFrom=1.0, CacheUpTo=16
//	(error less that 1e-8 for LogGamma)
// x= 16.0002 lgamma_1(x)= 2.74103 LogGamma_1(x)= 2.74103 err= 5.69013e-07
// x= 8.05173 lgamma_2(x)= 0.132245 LogGamma_2(x)= 0.132228 err= 1.71799e-05


static int NumLogGammaLookup=0;	// how many calls to LogGamma
static int NumLogGammaLog=0;	// how many calls to log
static int NumLogGammaMiss=0;	// how many calls to lgamma
static int NumLogGammaUsed=0;	// how many slots used

static int NumInverse=0;	// how many 1/x in scaling for derivatives
static int NumInverse2=0;	// how many 1/x^2 in scaling for derivatives

static const double EulerConst = 0.57721566490;
	// -EulerConst = Gamma'(1)/Gamma(1) = Gamma'(1)

// some constants for psi(x) taken from Gradsteyn and Ryzhik, others
//	computed by program.
static const double CacheFrom = 0.333333333;
static const double abs_psi_CacheFrom = 3.132034;

	// upper bound on derivative of lgamma at CacheFrom

// older values for CacheFrom
  //static const double CacheFrom = 2.0;
  //static const double abs_psi_CacheFrom = 0.433784;
  
  //static const double CacheFrom = 1.0;
  //static const double abs_psi_CacheFrom = EulerConst;
  
  //static const double CacheFrom = 0.8;
  //static const double abs_psi_CacheFrom = 0.96501;

  //static const double CacheFrom = 0.5;
  //static const double abs_psi_CacheFrom = 1.963510026;
  

static const double CacheUpTo = 128;
static const double abs_psi_CacheTo = 4.84812;
	// upper bound on derivative of lgamma at CacheTo
	// note: for x > approx 3, psi(x) < ln(x)

// older values for CacheUpTo
 //static const double CacheUpTo = 16;
 //static const double abs_psi_CacheTo = 2.74101;

  //static const double CacheUpTo = 32;
  //static const double abs_psi_CacheTo = 3.45003;

			
static const double MaxSlopeLogGamma = 
	(abs_psi_CacheFrom > abs_psi_CacheTo)? 
		abs_psi_CacheFrom: abs_psi_CacheTo;

static const double InverseInterval = MaxSlopeLogGamma/epsilon;
	// 1/ size of the interval between data points in the cache.
static const double InverseInterval2 = (InverseInterval*InverseInterval);

// create cache in a lazy manner. Use -1.0 to indicate an empty entry, since
// we know that the min of lgamma is around lgamma(1.46163) and is
// certainly larger than -0.121487
static const double LogGammaNull = -1.0;
static int LogGammaCacheSize = 0;
static double *LogGammaCache = NULL;

static void BuildCache()
{
    LogGammaCacheSize = static_cast<int>(CacheUpTo*InverseInterval) +3;
    LogGammaCache = new double[LogGammaCacheSize];
    for (int i = 0; i < LogGammaCacheSize; i++)
    {    LogGammaCache[i] = LogGammaNull;
    }
}

inline void EnsureCache()
{    if (LogGammaCache == NULL)
     {   BuildCache();
     }
}

inline static int CacheIndex(double x, double &frac)
{   int i= static_cast<int>(x*InverseInterval);
    frac = x*InverseInterval - i;
    return i;  
}

inline static double CacheLookup(int i)
{   assert(i>=0);
    assert(i<LogGammaCacheSize);
    if (LogGammaCache[i] <= LogGammaNull)
    {   NumLogGammaUsed ++;
        NumLogGammaMiss++;
	LogGammaCache[i] = lgamma(i* (epsilon/MaxSlopeLogGamma));
    }
    return LogGammaCache[i];
}

double LogGamma(double x)
{   EnsureCache();
    assert(x>0);
    NumLogGammaLookup++;
    
    double base=0.;
    while (x<CacheFrom)
    {   base -= log(x);
	x++;
	NumLogGammaLog++;
    }
    
    if (x>CacheUpTo)
    {   NumLogGammaMiss++;
	return lgamma(x)+base;
    }
    
    double frac;
    int i = CacheIndex(x, frac);
    
    double lowy = CacheLookup(i);
    double hiy = CacheLookup(i+1);

    return lowy + (hiy-lowy)*frac + base;
}

void LogGamma_print_summary(ostream & out)
{   EnsureCache();
    out << "\n";
    out << "LogGamma caching from " <<  CacheFrom << " to " << CacheUpTo
	<< " with tolerance " << epsilon << "\n";
    out << "LogGamma cache has " << LogGammaCacheSize
		 << " slots "  << 1./InverseInterval 
		 << " wide, of which " << NumLogGammaUsed 
		 << " were used\n";
    if (NumLogGammaLookup)
    {	out << "LogGamma cache had " << NumLogGammaLookup 
		<< " lookups,\n"
		<< "of which " << NumLogGammaMiss
		<< " (" << (100.*NumLogGammaMiss)/NumLogGammaLookup
		<< "%) required lgamma,"
		<< "and " << NumLogGammaLog
		<< " (" << (100.*NumLogGammaLog)/NumLogGammaLookup
		<< "%) required log.\n";
	out << "Derivative calculations used " << NumInverse 
		<< " 1/x computations"
		<< " and " << NumInverse2 << " 1/x^2 computations\n";

    }
}

// LogGamma_1(x)
//	is an approximation to the derivative of lgamma(x)
double LogGamma_1(double x)
{   EnsureCache(); 
    assert(x>0);
    assert (CacheUpTo- CacheFrom > 1.);	
	// need to have sufficient table so that pre-scaling doesn't
	// miss the table.
    NumLogGammaLookup++;
    
    double lg1=0.;
    while (x<CacheFrom)
    {   lg1 -= 1./x;
	NumInverse++;
	x++;
    }
    
    if (x > CacheUpTo)
    {	// use Stirling approx if off the top of the cache:
	double denom = 1./ (x*(12*x+1));
	lg1 +=  -0.5/(x) + log(x) - denom;
	return lg1;
    }
    
    double frac;
    int i = CacheIndex(x, frac);
    double y_0 = CacheLookup(i);
    double y_1 = CacheLookup(i+1);
    double y_minus = CacheLookup(i-1);

    return (0.5*(y_1-y_minus) + (y_1-2*y_0+y_minus)*frac)*InverseInterval
	+ lg1;
}

void LogGamma_derivs(double x, double &lg0, double &lg1, double &lg2)
{   EnsureCache(); 
    assert(x>0);
    assert (CacheUpTo- CacheFrom > 1.);	
	// need to have sufficient table so that pre-scaling doesn't
	// miss the table.
    NumLogGammaLookup++;
    
    lg0=0; lg1=0; lg2=0;
    while (x<CacheFrom)
    {   lg0 -= log(x);
	NumLogGammaLog++;
	lg1 -= 1./x;
	NumInverse++;
	lg2 += 1./(x*x);
	NumInverse2++;
	x++;
    }
    
    int InterpolateForLogGamma=1;
    if (x > CacheUpTo)
    {	lg0 += lgamma(x);
	NumLogGammaMiss++;
	
    // use Stirling approx if off the top of the cache:
	double denom = 1./ (x*(12*x+1));
	lg1 +=  -0.5/(x) + log(x) - denom;
	lg2 +=  0.5/(x*x) + 1./x  + (24*x + 1) * (denom *denom);
	return;
    }
    
    double frac;
    int i = CacheIndex(x, frac);
    double y_0 = CacheLookup(i);
    double y_1 = CacheLookup(i+1);
    double y_2 = CacheLookup(i+2);
    double y_minus = CacheLookup(i-1);

    if (InterpolateForLogGamma)
	lg0 += y_0 + (y_1-y_0)*frac;
    
    double mid = y_1-2*y_0+y_minus;
    lg1 += (0.5*(y_1-y_minus) + mid*frac)*InverseInterval;
    lg2 += (mid + (y_2-y_minus + 3*(y_0 - y_1))*frac)*InverseInterval2;

}

// CHANGE LOG:
// 1 Nov 1995 Kevin Karplus
//	Implemented cache with fixed spacing.
// 2 Nov 1995 Kevin Karplus
//	Implemented interpolating cache.
// 27 Dec 1995 Kevin Karplus
//	Implemented LogGamma_1 and LogGamma_derivs using cache and
//	Stirling approx.
// 15 March 2004 Kevin Karplus
//	Changed old-style cast to static_cast
// Fri Dec 14 11:37:40 PST 2007 Kevin Karplus
//	Mark Diekhans fixed the cache to be dynamically allocated
