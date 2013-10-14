#ifndef MLPREG_H
#define MLPREG_H

#include "Regularizer.h"
#include "Input.h"
#include "EqualStrings.h"
#include <iostream>

// maximum likelihood class (with pseudocounts from file) MLPReg

class  MLPReg : public Regularizer
{
	float *Pseudocounts;
	int FreezeDist;	// if set, only change scaling
  
	static IdObject ID;
	static NameToPtr* CommandTable;    
	
	void init_command_table(void);
	NameToPtr *command_table(void) {return CommandTable;}
	
	void write_knowing_type(ostream &out) const
	{   Regularizer::write_knowing_type(out);
	    out << "Pseudocounts =\n";
	    for (int i=0; i<alphabet_size(); i++)
		out << Pseudocounts[i] << "\n";
	}

  public: 
    MLPReg(void) : Regularizer() {Pseudocounts=0; FreezeDist=0;}
    MLPReg(const Alphabet *a, istream& in, const char* nm)
	: Regularizer(a,nm)
    {	Pseudocounts=0; FreezeDist=0;
	read_knowing_type(in);
    }
    MLPReg(const Alphabet *a, const float *ps, const char*nm=0)
	: Regularizer(a,nm)
    {	Pseudocounts=0; FreezeDist=0;
        set_pseudocounts(ps);
    }
    MLPReg(const AlphabetTuple *a, const float *ps, const char*nm=0)
	: Regularizer(a,nm)
    {	Pseudocounts=0; FreezeDist=0;
        set_pseudocounts(ps);
    }
    
    Regularizer* copy(void) const
    {    return new MLPReg(alphabet_tuple(), Pseudocounts, name());
    }
    ~MLPReg() {delete [] Pseudocounts;}
    
    static IdObject* classID(void) {return &ID;}
    virtual IdObject* type(void) const {return &ID;}
  
    const float* pseudocounts(void) const {return Pseudocounts;}
    void set_pseudocounts(const float *ps)
    {   if (!Pseudocounts) 
	    Pseudocounts = new float[alphabet_size()];
        for (int i=alphabet_size()-1; i>=0; i--)
	    Pseudocounts[i] = ps[i];
    }
    
    void get_modified_counts(
           const float* TrainCounts,   // what you use as counts
           float* probs)        // you fill this in with probs.
      {
    	    int i;
	    for (i=alphabet_size()-1; i >= 0; i--) 
    	    {      probs[i] = TrainCounts[i] + Pseudocounts[i];
    	    }
  	}

    void freeze_dist(void) {FreezeDist=1;}
    void unfreeze_dist(void) {FreezeDist=0;}
    protected:
    
	// if FreezeDist, there is only one parameter---the
	//	sum of the Pseudocounts,
	// otherwise each pseudocount is a separate parameter
	int num_parameters(void) const 	
		{return FreezeDist? 1: alphabet_size();}
	float parameter(int i) const 
	{   if (!FreezeDist) return Pseudocounts[i];
	    float sum=0.;
	    for (int j=alphabet_size()-1; j>=0; j--)
	    	sum+= Pseudocounts[j];
	    return sum;
	}
	void set_parameter(int i, float p) 
	{   if (!FreezeDist)
	    {	Pseudocounts[i] = p;
	    	return;
	    }
	    float scale=0.;
	    int j;
	    for (j=alphabet_size()-1; j>=0; j--)
	    	scale+= Pseudocounts[j];
	    scale = p/scale;
	    for (j=alphabet_size()-1; j>=0; j--)
	    	Pseudocounts[j] *= scale;
	}
	void partials1(float *part1, int i, const float* Counts)
	{   if (FreezeDist) 
	    {   part1[0] = 1.;
	    	return;
	    }
	    for (int j=alphabet_size()-1; j>=0; j--)
	    	part1[j]=0.;
	    part1[i] = 1.;
	}

};

// CHANGE LOG:
//  30 March 2004 Kevin Karplus
//	Fixed parameters that shadow members.
#endif
