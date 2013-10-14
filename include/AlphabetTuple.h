// AlphabetTuple.h
// created: 26 April 1995
// Kevin Karplus

// ABSTRACT
//	An AlphabetTuple is a tuple (pair, triple, ...) of Alphabet
//	It is intended for use with short fixed-length tuples, not sequences.
//	The associated bases are represented with a BaseTuple.
//	
//	Examples:
//	An AlphabetTuple of two RNA Alphabets can be used to represent
//	base pairs for folding RNA.
//
//	An AlphabetTuple containing 3 nucleic acid Alphabets can be
//	used for representing codons.
//
//	An AlphabetTuple containing 2 amino acid alphabets can be used
//	for residue contacts in proteins.

//  Question: should BaseTuple include a  pointer to the
//	corresponding AlphabetTuple? or should it just be raw Bases?
//	Should BasePair, BaseTriple, ... be created as more
//	space-efficient special cases?

#ifndef AlphabetTuple_H
#define AlphabetTuple_H

#include <iostream>
#include "Alphabet.h"
#include "NumericAlphabet.h"
#include "NamedObject.h"

class BaseTuple;

class AlphabetTuple: public NamedObject
{
  private:
  	int NumAlphabets;
	const Alphabet** alphabets;

        NumericAlphabet NumericAlph;

	int NumNormal;	// number of "normal" BaseTuples =
			// product of sizes of alphabets
        int NumNormWc; // number of nomal and wc BaseTuples in this AlphabetTuple 
  public:
				// construct tuple from
	AlphabetTuple(const Alphabet *a0);	// singleton Alphabet
	AlphabetTuple(const Alphabet *a0, 
			const Alphabet *a1);	 // pair
	AlphabetTuple(const Alphabet *a0, 
			const Alphabet *a1,
			const Alphabet *a2);	 // triple
	AlphabetTuple(int i, const Alphabet**a); // array of pointers
	AlphabetTuple(const AlphabetTuple& a);
        AlphabetTuple(const NumericAlphabet &na);
	
	~AlphabetTuple()
	{   delete [] alphabets;
	}

        int num_alphabets() const {return NumAlphabets;}
	int num_normal() const { return NumNormal; }
        int num_normal_wc() const {return NumNormWc;}
	const Alphabet* operator[](int i) const
	{   assert (0<=i && i<NumAlphabets);
	    return alphabets[i];
	}

	int same_as(const AlphabetTuple* other) const
        {   if (isNumeric() != other->isNumeric()) return 0;
            if (isNumeric())
            { return (NumericAlph.length() == other->NumericAlph.length());
            } else {
                if (num_alphabets() != other->num_alphabets()) return 0;
                for (int i=num_alphabets()-1; i>=0; i--)
                    if (alphabets[i] != other->alphabets[i]) return 0;
            }
	    return 1;
	}

        int isNumeric() const { return alphabets == 0; }
        void set_length(int n);
	
	void print_command(ostream & out) const;
	
	// convert BaseTuples to and from integer indices
	//	0<= index <NumNormal
	// Note: all Bases in the BaseTuple must be normal (not wildcard)
	int index(const BaseTuple& bt) const;
	BaseTuple* unindex(int index) const;	// creates a new BaseTuple
        void unindex(int index, BaseTuple &bt) const; // unindex into existing BaseTuple
	void print_unindex(ostream &out, int index) const;	
		// print BaseTuple without creating new one

        int norm_wc_index(const BaseTuple &bt) const;
        void norm_wc_unindex(int ind, BaseTuple &bt) const;
	
};

// Read AlphabetTuple from input (creates a new AlphabetTuple)
// If dimension is specified, then that many Alphabet names are read,
// otherwise the dimension is read first.
extern AlphabetTuple* read_AlphabetTuple(istream &in, int dim=0);

// Reads a command of the form
//	Alphabet= <alphabet_name>
//	AlphabetPair= <alphabet_name> <alphabet_name>
//	AlphabetTriple= <alphabet_name> <alphabet_name> <alphabet_name>
//	AlphabetTuple=  <number> <alphabet_name> ... <alphabet_name>
// as would be output by print_command
// If the firstword is not recognized, it is looked up as an alphabet name,
//	as if preceded by "Alphabet="
extern AlphabetTuple* read_AlphabetTuple_command(istream &in);


extern ostream& operator << (ostream&, const AlphabetTuple&);
// prints just the alphabet tuple---to print a whole command, use
//	member function print_command

class BaseTuple
{    
  private:
    	const AlphabetTuple *al;
    	Base *bases;
  public:
    	BaseTuple(const AlphabetTuple& a) 
	{   al= &a;
	    bases = new Base[al->num_alphabets()];
	}
        BaseTuple(const BaseTuple &bt)
	{
            al = bt.al;
            int num_dims = al->num_alphabets();
            bases = new Base[num_dims];
            for (int i = 0; i < num_dims; i++) bases[i] = bt[i];
        }

	~BaseTuple()
	{   delete [] bases;
	}
  
  	Base& operator [](int i)
	{   assert(0<=i && i<al->num_alphabets());
	    return bases[i];
	}
  	const Base operator [](int i) const
	{   assert(0<=i && i<al->num_alphabets());
	    return bases[i];
	}

        const AlphabetTuple *alphabet_tuple() const { return al; };

        int is_normal() const
	{
            int all_normal = 1;
            int num_dims = al->num_alphabets();
            int i;
            for (i = 0; i < num_dims; i++)
            {
                all_normal = all_normal && bases[i].is_normal();
            }
            return all_normal;
        }

	// just outputs Bases, not AlphabetTuples
	friend ostream& operator << (ostream&, const BaseTuple&);
	
	// must already have AlphabetTuple defined
	friend istream& operator >> (istream&, BaseTuple&);	

};


#endif

// CHANGE LOG:
// 20? Nov 1995  Kevin Karplus
//	Added same_as
// 24 Nov 1995 Kevin Karplus
//	added print_command and copy constructor
// Sat Jun 18 18:18:06 PDT 2005 Kevin Karplus
//	Made AlphabetTuple be a named object.
//	(The names are constructed by concatenating the Alphabet names,
//	separated by commas.)
