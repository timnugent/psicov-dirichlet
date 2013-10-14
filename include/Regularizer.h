// Regularizer.h

#ifndef Regularizer_H
#define Regularizer_H

#include "NamedClass.h"
#include "AlphabetTuple.h"
#include <iostream>
#include <assert.h>

#define bzero(s, n)	   memset ((s), 0, (n))
#include <string.h>  //for memset

class Regularizer;	// forward declaration

// class for keywords that can be used in Regularizer input
class RegInputCommand: public NamedObject
{
	typedef     int (*fcn)(istream &in, 
		    Regularizer *change,
		    RegInputCommand* self);

	fcn CommandFunction;

	// function to execute when keyword found.
	//	Reading from "in" into "change".
	//	Pass this down to function as 3rd arg, 
	//		so it can report error using self->name().
	//  Return 1 if input should continue, 0 if error or end of input.


    public:
	RegInputCommand(const char *nm, fcn c=0)
	{   set_name(nm);
	    CommandFunction=c;
	}
	
	int execute(istream &in, Regularizer *change)
	{    return (*CommandFunction)(in, change, this);
	}
};

int ReadComment(istream &in, Regularizer *change, RegInputCommand* self);
// function for RegInputCommand that treats keyword as a comment
//	and skips to end of line.

// abstract class for estimating distributions from observed samples
class Regularizer : public NamedClass, public NamedObject
{    
    AlphabetTuple *alph;
    static IdObject ID;
    static NameToPtr *CommandTable;
    int *InputOrder;	// order of symbols for input (from read_order)
    
    friend class OptimizeRegularizer;
    
  public:
    Regularizer(void) {alph=0; InputOrder=0;}
    Regularizer(const Alphabet *a, const char* nm=0)
    {	alph=new AlphabetTuple(a);
	set_name(nm);
	InputOrder=0;
    }
    Regularizer(const AlphabetTuple *a, const char* nm=0)
    {	alph=new AlphabetTuple(*a);
	set_name(nm);
	InputOrder=0;
    }
    virtual ~Regularizer(void) 
    {	delete alph;
	delete [] InputOrder;
    };

    virtual IdObject* type(void) const {return &Regularizer::ID;}
    static  IdObject* classID(void) {return &Regularizer::ID;}

    virtual Regularizer* copy(void) const=0;
    // make a new regularizer of the same type and copy all fields
    
    int alphabet_size() const {return alph->num_normal();}
    const AlphabetTuple* alphabet_tuple() const {return alph;}
    
    // Note: this doesn't copy the tuple, but the Regularizer now owns it.
    // set_alphabet_tuple should only be called with
    // set_alphabet_tuple(new AlphabetTuple ...)
    void set_alphabet_tuple(AlphabetTuple* a)
    {   if (InputOrder)
	{   cerr << "Can't change alphabet after InputOrder set\n"
		<< flush;
	    return;
	}
	delete alph; 
	alph= a;
    }
    void set_alphabet(const Alphabet *a=Alphabet::ret_default())
    {	set_alphabet_tuple(new AlphabetTuple(a));
    }
    
    void print_order(ostream&out) const;
	// print the order of the alphabet
	// (used in several "print" routines)
	// Note: this is NOT (at the moment) the same as the InputOrder.
    
    const int *input_order(void) const {return InputOrder;}
    void read_order(istream& in);
	// read an order for the alphabet and put it in InputOrder
    
    virtual void print_info(ostream& out) const
    {	out << type()->name() << " " << name();
    }
	// print some sort of useful message about your distrubution,
	// like the name and parameters.  This will probably be used
	// to identify what distribution was used in a given experiment
	// so make it useful.

    static Regularizer* read_new(istream& in,
		IdObject* required_type=&Regularizer::ID);
    
    static Regularizer* read_new(const char* filename,
		IdObject* required_type=&Regularizer::ID);

    virtual void get_modified_counts(
	const float* TrainCounts,   // what you use as counts
	float* ModifiedCounts) = 0; // you fill this in with counts
					// to normalize to get probabilities

    void get_probs(const float* TrainCounts,  float* probs);
	// calls get_modified_counts and normalizes the results to sum to 1
	
    float encodingCostForColumnCounts(
	const float* RealProbs,
	const float* TrainCounts,
	float *EstProbs);
	// Same as get_probs, but also
	// return the encoding cost in bits for a column whose
	// probabilities were estimated.
	// (calls get_probs)
    
    virtual void normalize(void) {return;}
	// normalize the parameters
	// (for those Regularizers that have an extra degree of freedom and
	//  may need to be normalized for numeric stability)
    
    int verify_partials1(const float*TrainCounts, float tolerance=0.01);
    int verify_partials2(const float*TrainCounts, float tolerance=0.01);
	// verify that the partial derivatives are within the specified
	// tolerance (for either absolute or relative error, whichever is
	// greater).  The standard of truth is a difference approximation
	// using the parameter incremented by a small amount, so don't try
	// setting the tolerance too small.
	// Return 1 if partials ok, report error and return 0 if not.
    
protected:
	virtual void write_knowing_type(ostream &out) const;
	int read_knowing_type(istream &in);
	
	// This static function is a no-op for the base class,
	// but is useful for setting up the is_a hierarchy for the 
	// derived classes.
	static void init_is_a(IdObject *self)
	{    self->add_is_a(Regularizer::classID());
	}

	virtual NameToPtr * command_table(void) {return CommandTable;}
	virtual void init_command_table(void);	

    int index_2D(int col_len, int row, int col)
    { return col_len*row + col; }

public: //These were originally protected, but a bug in cxx
        //preventing these members from being used in subclasses when
        //unless declared public
    // now many adjustable parameters are there in the model?
    virtual int num_parameters(void) const=0;
    virtual float parameter(int i) const=0;
    virtual void set_parameter(int i, float p)=0;

    // Return the first partial derivatives of posterior_counts[i]
    // with respect to each of the parameters in part[param] 
    // (where posterior_counts is what get_modified_counts would return)
    virtual void partials1(float *part, int i, const float *TrainCounts)=0;
	
    // return 1 if second partial derivative always zero
    virtual int zero_second_deriv(void) const;	
    
    // Return first two partials, default assumes second partial always 0
    virtual void partials2(float *part1, float *part2, int i,
	const float *TrainCounts); 
    
protected:
    // Note: must count DOWN, so that DirichletReg parameters set
    //	in proper order.
    virtual void set_param_vector(const float *v)
    {	for (int i=num_parameters()-1; i>=0; i--)
	    set_parameter(i, v[i]);
    }
    
    virtual void get_param_vector(float *v)
    {	for (int i=num_parameters()-1; i>=0; i--)
	    v[i] = parameter(i);
    }
    
    virtual int use_log(void) const {return 1;}
    // Use log of parameters when optimizing 
    // Regularizers that can have negative parameters should return 0 instead
    
    // limits on the value of each parameter.
    virtual float min_parameter(int i) const {return 1.e-7;}
    virtual float max_parameter(int i) const {return 100.;}

};


//CHANGE LOG:
// 2 May 1995 Kevin Karplus
//	Modified to use AlphabetTuple
// 29 May 1995 Kevin Karplus
//	Moved optimization out to OptimizeRegularizer class.
//	Added set_param_vector and get_param_vector
// 25 July 1995 Kevin Karplus
//	Eliminated use_log() member function (no longer used anywhere).
// 15 Nov 1995 Kevin Karplus
//	Modified to use NamedObject and NamedClass
// 24 Nov 1995 Kevin Karplus
//	Made pure virtual copy function.
// 16 Dec 1995 Kevin Karplus
//	Re-instituted use_log() member function---needed for optimizing
//	GribskovReg, which can have negative parameters.
// 26 Dec 1995 Kevin Karplus
//	Created default print_info, and used it in all the specific
//	print_info for derived classes.
// 4 Jan 1996 Kevin Karplus
//	Added verify_partials1
// 10 Jan 1996 Kevin Karplus
//	Added verify_partials2
// 17 Oct 1996 Spencer Tu
//      Added index_2D.
//      Changed uses of variable length local arrays to arrays
//      allocated by new and delete.  Affected methods were
//      verify_partials1, verify_partials2 and the function
//      ReadAlphabets.
#endif
