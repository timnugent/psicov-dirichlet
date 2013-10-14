// Multinomial.cc
// Kevin Karplus
// 23 September 1995

#include "Multinomial.h"
#include "LogGamma.h"

Prob Multinomial(const float *counts, int num_dim)
{
    double SumCounts=0;
    for (int i=num_dim-1; i>=0; i--)
	SumCounts += counts[i];
    return Multinomial(counts, num_dim, SumCounts);
}

Prob Multinomial(const float *counts, int num_dim, double SumCounts)
{
    double log_multi = LogGamma(SumCounts+1);
    for (int i = num_dim-1; i>=0; i--)
        if (counts[i]>0)
	{    log_multi -= LogGamma(counts[i] +1);
	}
    return static_cast<Prob>(static_cast<from_log>(log_multi));
}

// CHANGE LOG:
// 15 Nov 1995 Kevin Karplus
//	eliminated GammaCache, and now use LogGamma, instead of lgamma.
// 15 March 2004 Kevin Karplus
//	Changed old-style cast to static_cast
// 30 March 2004 Kevin Karplus
//	Changed test of counts[i] to >0
