// Alph.h
//
// Alph is a class providing access to a standard set of basic alphabets
//
//

#ifndef _Alph_h_
#define _Alph_h_

#include "Alphabet.h"
#include <string>

extern string std_alphabet;

/*! \class Alph Alph.h "Alph/Alph.h"
 *  \brief Access to a standard set of basic alphabets
 *
 * Alph is a class providing access to a standard set of basic alphabets
 */

class Alph {
private:
  static int NumAlphs;
  static const NucleicAlphabet  *NucleicPtr;
  static const RNAAlphabet      *RNAPtr;
  static const DNAAlphabet      *DNAPtr;
  static const ExtDNAAlphabet   *ExtDNAPtr;
  static const AminoAlphabet    *AAPtr;
  static const ExtAminoAlphabet *ExtAAPtr;


    static void InitAlph (void);
    static const Alphabet* init_and_lookup(const Alphabet* ptr, 
  					const char* name);
  
    static const Alphabet* check_ptr(const Alphabet *ptr, const char* name)
    {  	if (ptr) return ptr;
       	return init_and_lookup(ptr,name);
    }

public:
  static const NucleicAlphabet  &Nucleic() 
    {   NucleicPtr = dynamic_cast<const NucleicAlphabet*>(check_ptr(NucleicPtr, "Nucleic")); 
	return *NucleicPtr;
    }
	  
  static const RNAAlphabet      &RNA()
    {   RNAPtr =  dynamic_cast<const RNAAlphabet*>(check_ptr(RNAPtr, "RNA")); 
	return *RNAPtr;
    }
	  
  static const DNAAlphabet      &DNA()
    {   DNAPtr = dynamic_cast<const DNAAlphabet*>(check_ptr(DNAPtr, "DNA"));
	return *DNAPtr;
    }
	  
  static const ExtDNAAlphabet   &ExtDNA()
    {   ExtDNAPtr = dynamic_cast<const ExtDNAAlphabet*>(check_ptr(ExtDNAPtr, "ExtDNA")); 
	return *ExtDNAPtr;
    }
	  
  static const AminoAlphabet    &AA()
    {   AAPtr = dynamic_cast<const AminoAlphabet*>(check_ptr(AAPtr, "AA")); 
	return *AAPtr;
    }
	  
  static const ExtAminoAlphabet &ExtAA()
  	{ ExtAAPtr = dynamic_cast<const ExtAminoAlphabet*>(check_ptr(ExtAAPtr, "ExtAA")); 
	return *ExtAAPtr;
	}
	  

  // Set silent conversion
  static inline void silent_convert (int val = 1);
  // Default alphabet
  static inline void set_default (const Alphabet& def);
  static inline void set_default (const char *nm);
  static inline const Alphabet* ret_default (void);

  // Name to alphabet
  static const Alphabet* name_to_alphabet (const char *name);

  // Listing alphabets
  static inline int num_alphabets (void);
  static inline const Alphabet* alphabet (int i);

public:
friend class Alphabet;
};

// Alph Inlines
inline void
Alph::silent_convert (int val) { Alphabet::silent_convert(val); }

inline void
Alph::set_default (const Alphabet& def)
{ Alphabet::set_default(def); }

inline void
Alph::set_default (const char *nm)
{ Alphabet::set_default(nm); }

inline const Alphabet*
Alph::ret_default (void) { return Alphabet::ret_default(); }

inline int
Alph::num_alphabets (void) {return Alphabet::num_alphabets();}

inline const Alphabet*
Alph::alphabet (int i) 
{ return Alphabet::alphabet_list(i); }

inline const Alphabet*
Alph::name_to_alphabet (const char *name)
{  return Alphabet::name_to_alphabet(name, ZeroIfNew);	}


// CHANGE LOG:
// 27 Feb 2004 Kevin Karplus
//	Converted old-style casts.



#endif /* _Alph_h_ */
