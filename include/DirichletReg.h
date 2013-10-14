// Dirichlet Mixture method as a Regularizer

#ifndef DIRICHLETReg_H
#define DIRICHLETReg_H

#include "Regularizer.h"
#include "MLPReg.h"
#include "MLZReg.h"
#include "Prob.h"
#include "Multinomial.h"
 
#define DirichletReg_No_Extra_Caching

class  DirichletReg : public Regularizer
{
    static const int MaxPlus;	// largest value of "plus" for
    				// LogGammaComp cache
    
    int FreezeComp;	// set if optimization should only adjust mixture
    			// coefficients, not components.
    int FreezeMix;	// adjust only components, not mixture coefficients
    
			
    int NumComp, AllocComp;
    float *MixtureCoeff;
    float *Components;	// pseudocounts for components 
    			// use component(c)[letter] to access components
    // LogGammaComp is a cache of lgamma function of the components
    double *LogGammaComp;	// LogGamma(Components(c)[letter] + plus)
    double *Lgamma1Comp, *Lgamma2Comp;	// first two derivatives of 
    			// lgamma(components)

    double *SumComponents;	// sum of pseudocounts for cth component
    double *LogBetaComp;		// Beta(cth component) = 
    				// prod GammaComp / Gamma(SumComponents)
    double *Lgamma1Sum, *Lgamma2Sum; // derivatives of lgamma(SumComponents)


	float *component(int c)
	    { return Components + c*alphabet_size();}
    	double *lgamma1_comp(int c)
	    {return Lgamma1Comp + c*alphabet_size();}
    	double *lgamma2_comp(int c)
	    {return Lgamma2Comp + c*alphabet_size();}
	double &lgamma_component(int c, int letter, int plus=0)
	    {  
#ifdef  DirichletReg_No_Extra_Caching
	    return LogGammaComp[(c*alphabet_size()+letter)];
#else
	    return LogGammaComp[(c*alphabet_size()+letter)*(MaxPlus+1) + plus];
#endif
	    }

	double lgamma_component(int c, int letter, int plus=0) const
	    {  
#ifdef  DirichletReg_No_Extra_Caching
	    return LogGammaComp[(c*alphabet_size()+letter)];
#else
	    return LogGammaComp[(c*alphabet_size()+letter)*(MaxPlus+1) + plus];
#endif
	    }
	const double *lgamma_component(int c) const
	{    assert(0<=c && c <NumComp);
	    return LogGammaComp + c*alphabet_size();
	}
    
	// cached values set up by get_modified_counts
	//	also used by partials1 and partials2
	double SumCounts;	// cached sum of TrainCounts
	double *X;		// results of last get_modified_counts
	double *CompProbs;	// probability of each component given TrainCounts

	double *Lgamma1_SumCounts;	// lgamma1(SumCounts+SumComponents[c])
	double *Lgamma1_Counts;	// lgamma1(Counts[i] + component(c)[i])
	double &lgamma1_counts(int c, int i)
	{	return Lgamma1_Counts[c*alphabet_size()+i];
	}
	
	int Lgamma1CountsInvalid;
		// set if Lgamma1 not updated for latest Counts yet

	double *Lgamma2_SumCounts;	// lgamma2(SumCounts+SumComponents[c])
	double *Lgamma2_Counts;	// lgamma2(Counts[i] + component(c)[i])
	double &lgamma2_counts(int c, int i)
	{	return Lgamma2_Counts[c*alphabet_size()+i];
	}
	
	int Lgamma2CountsInvalid;
		// set if Lgamma2 not updated for latest Counts yet

	static IdObject ID;
	static NameToPtr* CommandTable;
    private:
	double log_beta_over_beta(const float* TrainCounts,
		double &SumTrainCounts,
		double* Log_Beta_Beta) const;
        // compute Beta(counts+z)/Beta(z) for each component
	// and return the max value (useful for normalizing)

	void compute_lgamma1(const float *Counts);
	// Fill in Lgamma1_Counts so that
	//	lgamma1_counts(c,k)=lgamma1(Counts[k] + component(c)[k])
	//	and Lgamm1_SumCounts[c] = lgamma1(SumCounts + SumComponents[c])

	void compute_lgamma2(const float *Counts);
	// Fill in Lgamma1_Counts and Lgamma2_Counts
	
	void write_knowing_type(ostream &out) const;
	void init_command_table(void);
	NameToPtr *command_table(void) {return CommandTable;}
	
	void delete_all(void);
	void clear(void);	// set all pointer values to 0

	// return 1 if x is a small integer, and set ix to that integer
	// ("small" means in the range usable for the gamma_component cache)
	static int is_small_int(float x, int &ix);    

    protected:
	int num_parameters(void) const 	;
	float min_parameter(int i) const ;
	float max_parameter(int i) const ;
	float parameter(int i) const ;
	void set_parameter(int i, float p);
	
	int zero_second_deriv(void) const {return 0;}
	
	// WARNING: must call get_modified_counts or get_moments first,
	//	to set up cache of component probabilities
	void partials1(float *part1, int i, const float* Counts);
	void partials2(float *part1, float*part2, int i, const float* Counts);

    public:
	DirichletReg(void);
	DirichletReg(const Alphabet *a, istream &in, const char *nm)
		: Regularizer(a,nm)
	{   clear(); 
	    read_knowing_type(in);
	}
	DirichletReg(const Alphabet *a, const char *nm, int size=0);
	DirichletReg(const DirichletReg& dir);	// copy constructor    
	DirichletReg(const MLPReg& dir);	// copy constructor    
	DirichletReg(const MLZReg& dir);	// copy constructor    

	Regularizer* copy(void) const {return new DirichletReg(*this);}

	DirichletReg* posterior_mixture(const float* TrainCounts);
	// creates a new mixture that is the correct posterior distribution
	// using this as a prior and seeing TrainCounts

	~DirichletReg() {delete_all();}
	
	static IdObject* classID(void) {return &ID;}
	virtual IdObject* type(void) const {return &ID;}

	void AddComponent(float MixCoeff,const float*comp);
	
	void print_info(ostream &out) const
	{   Regularizer::print_info(out);
	    out << " (" << num_components() << " components)";
	}

	void print_ordered_component(ostream &out, int c) const;
	// print the Alphabet sorted by how much each letter is favored
	//	by component (relative to background).
	// insert a comma every factor of 2, and >< at background freq.
	
	void get_modified_counts( const float* TrainCounts, 
		 float* probs);        // probs filled in with posterior counts
	// cache the component probabilites and sum of counts

	void get_moments(const float* TrainCounts, 
		 double* ex_prob, double *ex_prob2=0);
	// ex_prob filled in with expected value of probabilites 
	// ex_prob2 filled in with expected value of prob^2
	// cache the component probabilites and sum of counts
	
	int num_components(void) const {return NumComp;}
	void set_component(int c, int lett, float z);
	void scale_component(int c, float multiplier);
	void delete_component(int c);	// moves last component to cth position
	
	void set_mixture(int c, float mix)
	{   assert(0<=c && c <NumComp);
	    MixtureCoeff[c] = mix;
	}
	
	float mixture_coeff(int c) const
	{   assert(0<=c && c <NumComp);
	    return MixtureCoeff[c];
	}
	const float* mixture_coeff(void)const
	{    return MixtureCoeff;
	}
	double sum_component(int c) const
	{   assert(0<=c && c <NumComp);
	    return SumComponents[c];
	}
	const float *component(int c) const
	{    assert(0<=c && c <NumComp);
	    return Components + c*alphabet_size();
	}
	const float component(int c, int lett) const
	{    assert(0<=c && c <NumComp);
	    return Components[c*alphabet_size() + lett];
	}
	
	// make mixture coefficients sum to 1.
	void normalize(void);
	
    	int frozen_components(void) const {return FreezeComp;}
	void freeze_components(void) {FreezeComp=1;}
	void unfreeze_components(void) {FreezeComp=0;}
    	
	int frozen_mixture(void) const {return FreezeMix;}
	void freeze_mixture(void) {FreezeMix=1;}
	void unfreeze_mixture(void) {FreezeMix=0;}

	void component_probs(const float* TrainCounts,
		double &SumTrainCounts,
		double* comp_probs,
		double *log_sum=0);
	// Also return the probability for each component given the sample
	// comp_probs should be pre-allocated with num_components() floats
	// (also return the sum of the TrainCounts, for possible caching)
	// Optionally return
	//	ln sum_c mixture_coeff(c) * Prob(TrainCounts | component(c))
	//	without correction for permutations of the sample.
	// Used for computing the sample_weights array,
	//	also available for external use in optimizing 
	//	the mixture coefficients
	
	const double* component_probs(void) const {return CompProbs;}
	// return the component probabilities cached by the last
	// get_modified_counts, get_moments, or log_probability
	
	double log_probability(const float* TrainCounts,
		float *deriv1=0, float *deriv2=0);
	// returns the natural log of the probability of the count vector
	// (and up to 2 sets of partial derivatives, if the optional
	// arguments are non-zero pointers to arrays of
	// num_parameters() floats)
	// Probability NOT corrected for possible permutations of sample---
	// assumes count vector represents an ordered list of examples.
	
	inline double log_unordered_probability(const float* TrainCounts,
	       float *deriv1=0, float *deriv2=0)
	{   double log_pr=log_probability(TrainCounts, deriv1, deriv2);
	   return log_pr+
	       log(Multinomial(TrainCounts, alphabet_size(), SumCounts));
	}
	// same as log_probability, but with correction for
	// possible permutations of multisets.
	// Note that no correction needed for derivatives, since 
	//  additive correction is independent of paramters, so has deriv=0.


	Prob Probability(const float* TrainCounts)
	{   return static_cast<Prob>(static_cast<from_log>(log_probability(TrainCounts)));
	}
	
	Prob UnorderedProbability(const float* TrainCounts)
	{   double tmp = log_probability(TrainCounts);
	    // tmp used to force execution order, to make sure SumCounts ok
	    return Multinomial(TrainCounts, alphabet_size(), SumCounts)
		* static_cast<Prob>(static_cast<from_log>(tmp));
	}
	// return the probability of the observed count vector, given size
	// (note: assumes that order of observations is irrelevant, so
	// sum over all samples of size k should be 1)


	int verify_log_prob_partials(ostream &logfile,
		const float*TrainCounts, float tolerance);
	// Verify that the partials of log_probability 
	//	are correctly computed for the
	//	given count vector, by doing difference approximation.
	// Report errors to logfile.  Return 1 if no errors, 0 otherwise.

	
	void alloc(int size);	// create arrays with this much room
	// shouldn't really be public, but making it private made
	// the input routines fail.
};

// CHANGE LOG:
// 25 July 1995 Kevin Karplus
//	Added larger table of GammaComp, to allow faster computation
//	for small integer samples (also introduced MaxPlus).
//	Eliminated use_log function (no longer used).
// 2 August 1995 Kevin Karplus
//	Added copy constructor
// 20 Sept 1995 Kevin Karplus
//	Added ComponentProbs(), and changed old "sample_weight"
//		field to "CompProbs"
// 23 Sept 1995 Kevin Karplus
//	Separated out the beta_over_beta routine,
//	added Probability(TrainCounts)
// 31 Oct 1995 Kevin Karplus
//	made is_small_int inline for efficiency
// 1 Nov 1995 Kevin Karplus
//	Added #define DirichletReg_No_Extra_Caching 
//	to turn off caching of small integer gamma-components
// 15 Nov 1995 Kevin Karplus
//	Modified to use NamedObject and NamedClass
// 8 Dec 1995 Kevin Karplus
//	Moved void constructor to DirichletReg.cc and added code to
//	eliminate un-initialized  variables.
// 31 Dec 1995 Kevin Karplus
//	Eliminated Prob from internals, using double of log_probability
//	instead.  
//	Separated print_ordered_component out from write_knowing_type.
//	Restored lgamma1 and lgamma2 caches.
// 4 Jan 1996 Kevin Karplus
//	added caches for lgamma1 of counts+components
//	added flag to indicate caches invalidated by new component_probs
// 12 Jan 1996 Kevin Karplus
//	added log_probability, and defined Probability in terms of it.
// 15 Jan 1996 Kevin Karplus
//	added log_sum return value to component_probs(...)
// 16 Jan 1996 Kevin Karplus
//	made log_probability NOT correct for multinomials.
//	Correction done ONLY in Probability.
//	Eliminated UseFakeDeriv.
//	Added verify_log_prob_partials
// 26 Jan 1996 Kevin Karplus
//	Changed definition of Probability to refer to summary of
//		ordered set, added UnorderedProbability for old meaning
// 30 Jan 1996 Kevin Karplus
//	Added cast from MLZReg (to complement one from MLPReg)
// 7 Feb 1996 Kevin karplus
//	Added posterior_mixture()
// 8 Dec 1996 Kevin Karplus
//	Added frozen_components() and frozen_mixture();
// 15 March 2004 Kevin Karplus
//	Fixed old-style casts
// Thu Jul 21 08:55:11 PDT 2005 Kevin Karplus
//	added mixture_coeff(void) to return array

#endif
