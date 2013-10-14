/************************************************************************
DO NOT INCLUDE THIS FILE DIRECTLY!!! IT IS PARAMETERIZED SO THAT IT
IMPLEMENTS BOTH FLOAT AND DOUBLE LOG PROBS.
* Probability Class Guts functions.  A class for dealing with
* log-probabilities rather than the probabilities themselves,
* performing multiplication and division as addition and subtraction.
* leslie grate
* 4/20/94
* John Panzer
* 26 October 1993
*
* Richard Hughey
* 4 October 1993
************************************************************************/

/* Needed revisions:
   Variants with different precisions, setting of the 0 cutoff.
   Simple variant that does not allow addition (only comparison and
   multiplication.
   */
class OTHERPROB;
class PROBCONSTRUCT;

// Probstorage currently is just a PROB with no operators except
// casting to the full form.  It is useful if you need a prob that
// does not have a constructor for speed reasons.  

class PROBCLASS : public AbstractProb
{
protected:

  LOGTYPE val;

// Note: inlines defined here, so that they are defined before they
//	are used.
// Zeros and valids are set to maintain the same range in float and
// double probs.  Floats range to about 2E38, while doubles range to
// about 2E307.  We set logzero to be smaller than this: -2E35,
// and loginvalid at -2E37.  Note that when we multiply 2 probs, their
// logarithms are added.  If too many zeros were added together, they
// would become invalid, so we must check for zeros before adding.  If
// NaNs and infinities were checked in hardware (is_NaN() and
// is_inf()) rather than software as on many machines, this code could
// be further sped by making use of IEEE special numbers, eliminating
// the is_zero check from multiplication.

// The smallest IEEE denormal number is about exp(-713)
// or about... 1.0E-323.306215.  This is probably what should be
// used for epsilon, however printf will not print denormals, making them
// difficult to check, so instead the smallest normal is used.  It
// would be nice to use, for example, the Sun mathlib min_normal()
// function, but this is not defined everywhere, so we use the
// actual value instead.
// The smallest IEEE normal number is about 4.45E-308, or exp(-708),
// which produces 1.0E-307.65

  static inline LOGTYPE _logValid(void) {return -1e37;}
  static inline LOGTYPE _logInvalid(void) {return -2e37;} // An invalid value.
  static inline LOGTYPE _logZero(void) {return -2e35;}
  static inline LOGTYPE _logEpsilon(void) {return -708;}

// IEEE 32-bit has 24 bits of precision.  log(2^24) = 16.6, so the sum
// of 2 logs differing in greater than this amount will simply be the
// larger of the two.
  static inline LOGTYPE _maxDiff(void) {return LOGTYPEMAXDIFF;}
  
// approxDiff is used in approx equal.  Two numbers are approximately
// equal if the difference in their logs is smaller than approxDiff.
// It is the same value for both double and float-based logs

  static inline LOGTYPE _approxDiff(void) {return 1e-7;}
  static inline LOGTYPE _calcLog(double v)
  {  // Returns log of given value, cutting off too-small values.
	return (v < ret_epsilon()) ? _logZero() : log(v);
  }
public:
  static double ret_epsilon(void) {return 1.e-307;} 

  // Information functions.
  int valid(void) const {return (val > _logValid());} // representable?
  int non_prob(void) const {return !valid() || val > 0.0;}
  int is_zero(void) const {return (val <= _logZero());} 
  int is_one(void) const {return val>-1.e-100 && val<1.e-100;}             
  double ret_log(void) const {return val;}            // Returns (log(prob)).
  double ret_double(void) const
    {return (val < _logEpsilon()) ? 0.0 : exp(val);} // Returns prob [0,1]

  void set_log (double f) {val = f; PROB_INVARIANT_CHECK;}
  void set_double (double f) {val = _calcLog(f);
				    PROB_INVARIANT_CHECK;}
  void set_const (ConstProb c) {set_log (c == One ? 0 : c == Zero ?
				       _logZero() : _logInvalid());}
  
  // Entropy.  defined as -prob*log(prob).  since internal is log, looks different.
  double entropy(void) const {return (-(ret_double() * val)); };
  // Unary operators.
  PROBCLASS operator~(void) const 
     {PROBCLASS x; x.set_double(1.0-ret_double()); return x;}
     // Returns complement of a probability.

  // Arithmetic operators.
 inline PROBCLASS operator*=(PROBCLASS x); // {val += x.val; return *this;}
 inline PROBCLASS & operator/=(const PROBCLASS& x); // {val -= x.val; return *this;}
 inline PROBCLASS & operator+=(const PROBCLASS& x);
 inline PROBCLASS & operator-=(const PROBCLASS& x);

 // Borland chokes atinline modifier on a friend declaration...
 // DEC Alpha cxx requires it.
 // g++ can take it or leave it.
 friend INLINE PROBCLASS operator*(const PROBCLASS&,const PROBCLASS&);
 friend INLINE PROBCLASS operator/(const PROBCLASS&,const PROBCLASS&);
 friend INLINE PROBCLASS operator+(const PROBCLASS&,const PROBCLASS&);
 friend INLINE PROBCLASS operator-(const PROBCLASS&,const PROBCLASS&);
 friend INLINE int approxEqual(PROBCLASS const &,PROBCLASS const &);
 friend INLINE PROBCLASS pow(PROBCLASS const &, double);
 friend INLINE PROBCLASS pow(PROBCLASS const &a, int i) {return(pow(a, static_cast<double>(i)));};
 friend INLINE double pow_double(PROBCLASS const &, double);
};



inline double
log2(PROBCLASS const & a)
// Returns base-2 logarithm of a Prob.
{ 
  return a.ret_log() / M_LN2;
}

inline double
log(PROBCLASS const & a)
// Returns natural logarithm of a Prob.
{ 
  return a.ret_log();
} 

// min, max, mini, and maxi functions:
inline PROBCLASS
min(PROBCLASS const & a, 
    PROBCLASS const & b)
{
  if (b.ret_log() < a.ret_log())
    return b; // Returns copy of b if b < a.
  else
    return a; // Returns copy of a if a <= b.
} 

inline PROBCLASS
max(PROBCLASS const & a,
    PROBCLASS const & b)
{
  if (b.ret_log() > a.ret_log())
    return b; // Returns copy of b if b > a.
  else
    return a; // Returns copy of a if a >= b.
} 


// X= operators:
inline PROBCLASS 
PROBCLASS::operator*=(PROBCLASS x) 
{ 
  PROB_INVARIANT_CHECK;
  PROBCLASS y;
#ifdef PROB_CHECK_ZERO  
  if (is_zero() || x.is_zero())
    val = _logZero();     
  else
#endif    
  y.val = (val += x.val); 
  PROB_INVARIANT_CHECK;
  return y;
} 

inline PROBCLASS &
PROBCLASS::operator/=(const PROBCLASS& x) 
{
  PROB_INVARIANT_CHECK;
  assert(!x.is_zero()); // Disallow division by zero.
#ifdef PROB_CHECK_ZERO
  if (is_zero())
    val = _logZero();
  else
#endif    
    val -= x.val; 
  PROB_INVARIANT_CHECK;
  return *this;
} 


// X ? Y operators:
inline PROBCLASS 
operator*(const PROBCLASS& xarg,
          const PROBCLASS& y)
{
  PROBCLASS x =xarg;
  return   x *= y;
}  // !!! Make efficient later.

inline PROBCLASS 
operator/(const PROBCLASS& xarg,
          const PROBCLASS& y)
{
  PROBCLASS x=xarg;
  x /= y;
  return x;
}  // !!! Make efficient later.

inline PROBCLASS 
operator+(const PROBCLASS& x,
          const PROBCLASS& y)
{
  if (x.is_zero()) 
    return y;
  else if (y.is_zero()) 
    return x;

  assert(isfinite(x.val));
  assert(isfinite(y.val));
  
  double val, dif;

  if (x.val > y.val) {
    val = x.val;
    dif = y.val - x.val;
  } else {
    val = y.val;
    dif = x.val - y.val;
  }
  assert (dif <= 0);
  // maxDiff is the absolute value of the max log difference that
  // can make a difference in the underlying representation when
  // exp-log-summed. Above, diff is calculated as a negative quantity.
  if (-dif < PROBCLASS::_maxDiff())
    val += log (1.0 + exp(dif));

  //= ((x.val>y.val) ?    x.val + log(1.0+exp(y.val-x.val)) :
  //                                 y.val + log(1.0+exp(x.val-y.val)));
  PROBCLASS P;
  P.val = val;
  return P;
}
inline PROBCLASS 
operator-(const PROBCLASS& x,
          const PROBCLASS& y)
{
  if (y.is_zero())
    return x;
  if (x.val < y.val)  {
    if (approxEqual (x, y)) {
      PROBCLASS tmp;
      tmp.set_const(PROBCLASS::Zero);
      return tmp;
    }
    assert(x.val >= y.val);
  }
  // Make this smarter later:
  double val = exp(x.val) - exp(y.val);
  PROBCLASS P;
  P.val = PROBCLASS::_calcLog(val);
  return P;
} 

inline PROBCLASS &
PROBCLASS::operator+=(const PROBCLASS& x)
{ 
  return (*this = *this + x);
}  // !!! Make efficient later.

inline PROBCLASS &
PROBCLASS::operator-=(const PROBCLASS& x)
{ 
  return (*this = *this - x);
}  // !!! Make efficient later.

// Misc. operators.
inline std::ostream& 
operator<<(std::ostream& o,PROBCLASS const &p)
{ return p.valid() ? (o << p.ret_double()) : (o << "??");} 

PROB_DEF_RELOP(PROBCLASS,==)
PROB_DEF_RELOP(PROBCLASS,<=)
PROB_DEF_RELOP(PROBCLASS,>=)
PROB_DEF_RELOP(PROBCLASS,<)
PROB_DEF_RELOP(PROBCLASS,>)
PROB_DEF_RELOP(PROBCLASS,!=)
// (Translation of code above:  Apply operators to log of prob instead of
//  probability itself.)

inline int
approxEqual(PROBCLASS const & a, PROBCLASS const & b)
// Tells if two logprobs are approximately equal.
{
  if (a.is_zero() && b.is_zero())
    return 1;
  if (a>=b)
    return (a.val-b.val) < PROBCLASS::_approxDiff();
  else
    return (b.val-a.val) < PROBCLASS::_approxDiff();
}
// These are needed to avoid type confusion over
// the various from_double casts.
inline int
approxEqual (PROBCLASS const &a, from_double const &b)
{   PROBCLASS tmp; 
    tmp.set_double(b.val);
    return approxEqual (a, tmp);
}

inline int
approxEqual (from_double const &b, PROBCLASS const &a)
{   PROBCLASS tmp;  
    tmp.set_double(b.val);
    return approxEqual (a,tmp);
} 

inline PROBCLASS
pow(PROBCLASS const &a, double r)
// raise a to the r power
{
  PROBCLASS tmp;
  if (r >= 0.0)  
  {   tmp.set_log(a.val * r);
      return tmp;
  }
// can't raise prob to negative power and return valid prob!
  tmp.set_log(PROBCLASS::_logInvalid());
  return tmp;
}
  
inline double
pow_double(PROBCLASS const &a, double r)
// raise a to the r power
{
  double x;
  if (r >= 0.0) {
    x = (a.val * r); // is negative
    if (x < PROBCLASS::_logEpsilon()) {
// is TOO SMALL!
      return(PROBCLASS::ret_epsilon());
    } 
    return(exp(x));
  }
// negative power.
  x =  (a.val * r);
  if (x > -PROBCLASS::_logEpsilon()) {
// is TOO large !
      return(-PROBCLASS::ret_epsilon());
  } 
  return(exp(x));
}


class PROBCONSTRUCT : public PROBCLASS {
public:
  inline PROBCONSTRUCT(const OTHERPROB p);
  inline PROBCONSTRUCT(PROBCLASS const & OtherProb) 
         {set_log (OtherProb.ret_log());}
  inline PROBCONSTRUCT(void)   {set_log (PROBCLASS::_logZero());}
  inline PROBCONSTRUCT(from_double ProbVal) {set_double (ProbVal.val);}
  inline PROBCONSTRUCT(from_log logVal) {set_log (logVal.val);}
  inline PROBCONSTRUCT(ConstProb c)  {set_const (c);}
private:
  // Make a compilation-time error for attempting to invoke a prob on
  // a double.  rph. 
  inline PROBCONSTRUCT (const double d) {assert (0);}
};

// CHANGE LOG:
// 27 Feb 2004 Kevin Karplus
//	removed old-style casts
// 30 March 2004 Kevin Karplus
//	Changed =0 test for is_one to allow +- tiny error.
// Fri Sep  8 12:43:39 PDT 2006 Kevin Karplus
//	Added isfinite assertions to operator +
