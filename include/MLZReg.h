#ifndef MLZReg_H
#define MLZReg_H

#include "Regularizer.h"

// maximum likelihood class (with zero offset) MLZReg

class  MLZReg : public Regularizer
{
	float ZeroOffset;
	static IdObject ID;
      	static NameToPtr* CommandTable;    
	
	void init_command_table(void);
	NameToPtr *command_table(void) {return CommandTable;}

  public: 
    MLZReg(void) : Regularizer() {ZeroOffset=1e-4;}
    MLZReg(const Alphabet *a, istream& in, const char* nm)
    	: Regularizer(a, nm)
    {    read_knowing_type(in);
    }    
    MLZReg(const Alphabet *a, float z=1e-4, const char* nm=0)
    	: Regularizer(a, nm)
    {	set_zero_offset(z);
    }
    MLZReg(const AlphabetTuple *a, float z=1e-4, const char* nm=0)
    	: Regularizer(a, nm)
    {	set_zero_offset(z);
    }
    Regularizer *copy(void) const
    {    return new MLZReg(alphabet_tuple(), zero_offset(), name());
    }
    
    static IdObject *classID(void) {return &ID;}
    virtual IdObject* type(void) const {return &ID;}

    void set_zero_offset(float z) {ZeroOffset=z;}
    float zero_offset(void) const {return ZeroOffset;}
    
  
    void print_info(ostream &out) const
    {	Regularizer::print_info(out);
	out << " zero_offset= " << ZeroOffset;
    }
  
    void get_modified_counts(
           const float* TrainCounts,   // what you use as counts
           float* modified_counts)    
    {  	for (int i=alphabet_size()-1; i >= 0; i--) 
	    modified_counts[i] = TrainCounts[i] + ZeroOffset;
    }

    protected:
	int num_parameters(void) const {return 1;}
	float parameter(int i) const {return ZeroOffset;}
	void set_parameter(int i, float p) {ZeroOffset=p;}
	void partials1(float *part1, int i, const float* Counts)
	{   part1[0] =1.;
	}

	void write_knowing_type(ostream &out) const
	{   Regularizer::write_knowing_type(out);
	    out << "ZeroOffset = " << ZeroOffset << "\n";
	}

};

#endif

// CHANGE LOG:
// 15 Nov 1995 Kevin Karplus
//	Modified to use NamedObject and NamedClass
// CHANGE LOG:
//  30 March 2004 Kevin Karplus
//	Fixed parameters that shadow members.
