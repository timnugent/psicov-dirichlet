// Probability Package:  Prob.h ProgGuts.h LargeRealGuts.h
//
// Original design by Richard Hughey  rph@ce.ucsc.edu
//
// Latter modifications and improvements by Don Fong, Leslie Grate,
// John Panzer, and Richard Hughey


#ifndef _BIOPROBHEADER_H
#define _BIOPROBHEADER_H

#include <math.h>
#include <assert.h>
#include <iostream>
#include <iomanip>	// for setw

// (Usually) hidden class used to keep compiler from
// helpfully constructing temporary objects.  Under
// CC, you might have to use this to cast doubles
// to probabilities.
class from_double
{
public:
  double val;
  from_double(const double d) {val = d;}
};

// ditto.  if you have a value that is ALREADY a log prob, use
// this class to cast into a prob.
class from_log
{
 public:
  double val;
  from_log(const double d) {val = d;}
};

// Compiler bug workarounds:
#ifdef __BCPLUSPLUS__
#define INLINE 
#else
#define INLINE inline
#endif

class AbstractProb
// Base class for all probability classes.  Right now it just
// defines some useful constants & really basic items.
{
public:
  // Prob and LargeReal constants
  typedef enum {One=1, Unity = One, Zero = 0,Invalid = -1} ConstProb;
};

#ifndef PROB_INVARIANT_CHECK
#define PROB_INVARIANT_CHECK /*assert(OK())*/
#endif

// Logical operators.  Following macro used to keep typing down.
#define PROB_DEF_RELOP(PROBCLASS,OP) \
 \
inline int operator OP(const PROBCLASS & x,const PROBCLASS & y) \
{ return x.ret_log() OP y.ret_log(); } 


// Please see ProgGuts.h for an explaination of LOGTYPEMAXDIF and
// other constants.
#define LOGTYPE double
#define PROBCONSTRUCT Prob
#define PROBCLASS ProbBase
#define OTHERPROB ShortProb
#define LOGTYPEMAXDIFF 37.4
#include <ProbGuts.h>

#undef OTHERPROB
#undef LOGTYPE
#undef PROBCONSTRUCT
#undef PROBCLASS
#undef LOGTYPEMAXDIFF

// NOW do probs using floats
#define LOGTYPE float
#define PROBCONSTRUCT ShortProb
#define PROBCLASS ShortProbBase
#define OTHERPROB Prob
#define LOGTYPEMAXDIFF 16.6
#include <ProbGuts.h>

#undef LOGTYPE
#undef PROBCLASS
#undef PROBCONSTRUCT
#undef OTHERPROB
#undef LOGTYPEMAXDIFF

inline Prob::Prob(const ShortProb p) {set_log (p.ret_log());}
inline ShortProb::ShortProb(const Prob p) {set_log (p.ret_log());}

#undef PROB_DEF_RELOP


#endif




