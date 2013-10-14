//
// Alph.cc
//

#include "Alph.h"
#include <iostream>
#include <string>

int Alph::NumAlphs = 0;
const NucleicAlphabet  *Alph::NucleicPtr = NULL;
const RNAAlphabet      *Alph::RNAPtr = NULL;
const DNAAlphabet      *Alph::DNAPtr = NULL;
const ExtDNAAlphabet   *Alph::ExtDNAPtr = NULL;
const AminoAlphabet    *Alph::AAPtr = NULL;
const ExtAminoAlphabet *Alph::ExtAAPtr = NULL;
string std_alphabet = "data/std.alphabet";


const Alphabet* Alph::init_and_lookup(const Alphabet* ptr, 
  					const char* name)
{  
    InitAlph();
    ptr = Alphabet::name_to_alphabet(name);
    if (!ptr)
    {    cerr << "Error: Alphabet " << name << " not found in " 
    		<< std_alphabet
      		<< "\n" << flush;
    }
    return ptr;
}
      
void
Alph::InitAlph (void)
{
  //Ensure that at most one meaningful Alph is ever created
  if (NumAlphs == 0)
  {
      int num_loaded = Alphabet::load_alphabet_file(std_alphabet.c_str());
      if (num_loaded==0)
      {    cerr << "Warning: No Alphabets found in " << std_alphabet
      		<< "\n" << flush;
      }

      NucleicPtr = dynamic_cast<const NucleicAlphabet *>(Alphabet::name_to_alphabet("Nucleic"));
      RNAPtr = dynamic_cast<const RNAAlphabet *>( Alphabet::name_to_alphabet("RNA"));
      DNAPtr = dynamic_cast<const DNAAlphabet *>( Alphabet::name_to_alphabet("DNA"));
      ExtDNAPtr = dynamic_cast<const ExtDNAAlphabet *>( Alphabet::name_to_alphabet("ExtDNA"));
      AAPtr = dynamic_cast<const AminoAlphabet *>( Alphabet::name_to_alphabet("AA"));
      ExtAAPtr = dynamic_cast<const ExtAminoAlphabet *>( Alphabet::name_to_alphabet("ExtAA"));
      NumAlphs++;
  }

  /*****Use this code instead of above to print out standard alphabets*******
  NucleicPtr = new const NucleicAlphabet("Nucleic");
  RNAPtr = new const RNAAlphabet("RNA");
  DNAPtr = new const DNAAlphabet("DNA");
  ExtDNAPtr = new const ExtDNAAlphabet("ExtDNA");
  AAPtr = new const AminoAlphabet("AA");
  ExtAAPtr = new const ExtAminoAlphabet("ExtAA");
  *************************************************************/
}



// CHANGE LOG:
// 15 March 2004 Kevin Karplus
//	Changed old-style casts to static_cast, dynamic_cast, or
//	const_cast, as appropriate.
