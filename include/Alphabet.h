// Alphabet.h
// Revised 21 July 1994 rph@ce.ucsc.edu
// 3 Nov 1993 Richard Hughey rph@ce.ucsc.edu
//
// Defines classes for alphabet-independent biosequence analysis.
// Similar formats are used for nucleotides and proteins so that the
// basic data for either can be represented as a "Base".  Base
// includes provisions for wildcards and character variants as well as
// a non-matching base.  In this first version (Nov 3 1993), character
// variants are generally ignored.  Alphabet and its descendents
// define the actual alphabets (amino acids, nucleotides, and so on).
// These are instantiated (to enable pointers) in 'Alph', which
// has one static class member for each alphabet.  Thus, clients can
// refer to Alph::RNA.complement(), for example.

// Additional comments:
//    1.  == between bases is not defined in favor of forcing the
//    client to consider whether or not wildcards are being used.
//    2.  The Alphabet enum descendents should probably be static
//    const Base members rather than ints which get cast to Base.
//    However, for the time being, this is simpler and may be more
//    compiler-efficient.
//    3.  Consistancy checks have not been implemented.  As only one
//    object ever exists of each class, it is not clear that
//    consistency checks are strictly needed.

#ifndef _ALPHABET_H_
#define _ALPHABET_H_

#include <assert.h>
#include <iostream>
#include <string.h>

#include "NameToPtr.h"
#include "NamedObject.h"
#include "NamedClass.h"

class Base;
class Alphabet;

// Whether or not base variants should be considered (or reduced to
// their canonical forms).  Used by the match routines of Base and
// Alphabet. 
enum VarEnum  {NO_VARS, VARS};	
  // Options are using variant characters and on a non-wildcard match
  // whether or not to check for NULL characters.
enum NullEnum {NULL_NULL_FALSE, NULL_NULL_TRUE};

// Options on wildcard-wildcard matches:  does either wildcard match
// the other (SYMMETRIC) or is only the first wildcard's translation
// table checked.
enum MatchEnum {SYMMETRIC_MATCH, SUBSET_MATCH};


// Alphabet-dependant members of Base should be propagated to
// sequences to hide the Alphabet.
class Base
{
  private:
    // alphabet-less casts.  These are now private, only used by
    // the alphabet class to help prevent programmer errors.
    //      operator int (void) const {return _b;}
    //      Base& operator = (const int b) {_b = b; return *this;}
    // Contructor on int returns base #b in the [unknown] alphabet.
    // Base (const int b) : _b (b) {;};
    // No longer have a cast from int, to avoid accidental type changes.
    // Instead we have the private static function base_int()

    typedef unsigned char BaseType;
    BaseType _b;

    int operator == (const Base rhs) const; //  Use no_wc_match
  
  public:
    enum {
        _lowest_wc = 64,	// encode WC above this number.
	_first_var = 128,	// must be power of 2
	_max_chars = 256,	
	_max_wc    = _first_var - _lowest_wc -1
    };

  protected:
    enum {_null_base = (_first_var-1),	// Non-matching base.  Set here rather
	  // than 0 for indexing reasons.
	  _var_mask  = (_first_var-1)      // 2-bits for variants, 6 for chars
    };
    // Although Alphabet is friend, the above enum (in g++) is not
    // available to the class Alphabet prototype, only to the class
    // Alphabet functions.  Is this a bug, or just a quirk of C++?

  public:
    inline Base(void) {};
    inline Base (const Base &from) { _b = from._b; } // dk 26apr96
    static inline Base
	base_int (const int i)  {Base b;  b.set_int (i);  return b;}


    // Canon returns the canonical form (no variants) of a base
    inline Base canon (void) const
    {return base_int ((_b & static_cast<BaseType>( _var_mask)));}

    // Is_wild and is_normal is orthogonal to is_variant.  Thus, a
    // variant of a normal character will return 1 for is_normal.
    inline int is_wild (void) const 
    {   return  !is_null() && canon()._b >= _lowest_wc;
    }
    inline int is_normal (void) const{return canon()._b  <  _lowest_wc;}
  
    inline int is_variant (void) const{return canon()._b  != _b;}

    inline int raw_int (void) const {return _b;}     // Return the raw integer value
    inline void set_int (const int i){_b = i;} // Set the raw integer value

    static inline Base null (void) {return base_int (_null_base);}
    inline int is_null (void) const {return _b == _null_base;}
    
    static inline char null_char (void){return '-';}	   // ASCII for the null character
    static inline int  is_null_char (char c) {return c == null_char();}

    static inline char bad_char  (void){return '?';}	   // ASCII for a bad character.

    static void limits (ostream& s); // Print constants for limits on
					// numbers of characters on an ostream.
    
    // We do not overload == so that the programmer is forced to
    // consider whether wildcard or non-wildcard matches are used.
    // Match 2 bases with wild cards. 
    inline int wc_match  (const Base b2, const Alphabet *a,
			  const VarEnum v = NO_VARS) const;
    
    int wc_subset  (const Base b2, const Alphabet *a,
			   const VarEnum v = NO_VARS) const;
    
			
    // Match 2 bases without wild cards.
    inline int no_wc_match (const Base b2,
			    const NullEnum null_check = NULL_NULL_FALSE,
			    const VarEnum v = NO_VARS) const
    {   if (null_check == NULL_NULL_FALSE && (is_null() || b2.is_null()))
		return 0;
        if (v == NO_VARS)
		return (b2.canon()._b == canon()._b);
	return (b2._b == _b);
    }
    static inline int no_wc_match (const Base b1, const Base b2,
				   const NullEnum nc = NULL_NULL_FALSE,
				   const VarEnum v = NO_VARS)
    {	return b1.no_wc_match(b2, nc, v);
    }

    friend class Alphabet;

    friend inline ostream& operator<< (ostream& outs, const Base b);
    friend inline istream& operator>> (istream& ins, Base &b);
  
};

/* ********************************************************************** */

class AlphabetInputCommand : public NamedObject
{
        typedef int (*AlphabetIcf) (istream &in, Alphabet *change, AlphabetInputCommand *self);

        AlphabetIcf CommandFunction;
    public:
        AlphabetInputCommand(const char *nm, AlphabetIcf icf)
        {
            set_name(nm);
            CommandFunction = icf;
        }

        inline int execute(istream &in, Alphabet *change)
        {
            return (*CommandFunction) (in, change, this);
        }
};

class Alphabet : public NamedObject, public NamedClass 
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(NamedClass::classID());
    }
    int read_knowing_type(istream &in);
    void write_knowing_type (ostream &os) const;

    static NameToPtr *CommandTable;
    static void init_command_table(void);

private:
    static NamedClass *create_Alphabet(void);

protected:
    static int ReadName(istream &in,
                        Alphabet *change,
                        AlphabetInputCommand *self);
    static int ReadComment(istream &in, 
                Alphabet *change,
                AlphabetInputCommand* self);
    static int VerifyClassName(istream &in, 
                Alphabet *change,
                AlphabetInputCommand* self);
    static int ReadNormalChars(istream &in,
                               Alphabet *change,
                               AlphabetInputCommand *self);
    static int ReadAlias(istream &in, Alphabet *change, AlphabetInputCommand *self);
    static int ReadWildcard(istream &in, Alphabet *change, AlphabetInputCommand *self);
    static int ReadAllMatch(istream &in, Alphabet *change, AlphabetInputCommand *self);
    static int ReadCharName(istream &in, Alphabet *change, AlphabetInputCommand *self);
    static int ReadIsNucleic(istream &in, Alphabet *change, AlphabetInputCommand *self);
    static int ReadIsRNucleic(istream &in, Alphabet *change, AlphabetInputCommand *self);
    static int ReadIsAmino(istream &in, Alphabet *change, AlphabetInputCommand *self);

    void print_normal_chars(ostream &out) const;
    void print_wildcards(ostream &out) const;
    void print_alias_chars(ostream &out) const;
    void print_char_names(ostream &out) const;
    void print_is_values(ostream &out) const;

public:
    static IdObject *classID() { return &::Alphabet::ID; }
    virtual IdObject *type(void) const { return &::Alphabet::ID; }

public:
  // Information functions.
  int is_nucleic  (void) const;
  int is_rnucleic (void) const;
  int is_amino    (void) const;
  
  int norm_length (void) const;	// Number of normal characters
				// Always contiguous from first_char.
  int wc_length   (void) const;	// Number of wc characters.
  int norm_wc_length (void) const; // Number of normal and wc characters.
  int max_num_base  (void) const; // Max number chars not including variants.  
  int max_num_var_base(void) const; // Total number chars including variants.
				    // Note representation can be sparse.
  
  int first_char (void) const;  // Index of first normal character.
  int last_char  (void) const;
  int first_wc   (void) const;  // First wc, can be != first_char + length.
  int last_wc    (void) const;
  int first_var  (void) const;  // Index of first variant character.
  int last_var   (void) const;

    inline Base first_all_match_wildcard(void) const
  	// return (as a Base) the first wildcard that matches everything.
    {   Base wildcard;
	for (int i = norm_length(); i < norm_length() + wc_length(); ++i)
    	{   
		if (is_all_match(unindex(i)))    return unindex(i);
	}
	return Base::null();
    }


  // Return an index, possibly more compact than the Base to int cast.
  // Ignores variants, which are sparse.
  int  index   (const Base b) const; 
  int index_from_char(const char ch) const;
  Base unindex (const int i) const; 
  char unindex_to_char(const int i) const;

  // Matching function (no_wc_match requires no alphabet info).
  int wc_match  (const Base b1, const Base b2,
		 const VarEnum v = NO_VARS) const;
  
  // wc_subset returns 1 if b2 is a subset of b1.
  int wc_subset  (const Base b1, const Base b2,
		 const VarEnum v = NO_VARS) const;
  
  // Match 2 bases without wild cards.
  static inline int no_wc_match (const Base b1, const Base b2,
			  const NullEnum nc=NULL_NULL_FALSE,
			  const VarEnum v=NO_VARS)
    {	return b1.no_wc_match (b2, nc, v);
    }

  // ASCII conversion and alphabet name.
  /*inline*/ char to_char (const Base b1, const VarEnum v = NO_VARS) const;
  char to_ascii(const Base b1, const VarEnum v) const;
  Base to_base  (const char c) const;
  Base to_base  (const char *short_name) const;	
  	//  note: short_name never case-sensitive
  const Alphabet * const id (void) const;

  // Return a null-base-terminated list of non-wildcards the given
  // Base matches. Variant versions not available.
  inline const Base *matches (const Base b1, const VarEnum v) const;
  inline int num_matches (const Base b) const;
  
  // Return a Base or a Null base if the character is invalid for the
  // given alphabet.  
  Base valid_or_null (const Base b) const;
  int is_valid (const Base b) const;
  int is_valid_char(const char ch) const;
  Base null (void) const;

  // Return an abbreviated or full name of a base.
  // Returns NULL (not empty string) if no name has been assigned for this character.
  const char *abbrev (const Base b) const;
  const char *full_name (const Base b) const;

  int is_case_sensitive() const;
  void set_case_sensitive(int cs);
  char if_case_sensitive_char(const char ch) const;

  static inline void silent_convert (int val = 1);

  static const Alphabet *name_to_alphabet(const char *nm,
                                          OptionIfNew ifnew=ErrorIfNew);
  static const Alphabet *alphabet_list(int i);
  static int num_alphabets();

  // Describe the alphabets on an ostream.
  void describe (ostream& s) const;
  static void describe_all (ostream& s);
  static const NameToPtr *name_to_alphabet_table();
  static int load_alphabet_file(const char *filename);
  // Default alphabet
  static inline void set_default (const Alphabet& def);
  static inline void set_default (const char *nm);
  static inline const Alphabet* ret_default (void);

  int is_all_match(Base wcb) const;


protected:
  // Add an alias to a character.
  void add_alias (const char newchar, const char alias);
  // Add normal characters to the alphabet (used by descendent classes)
  void add_normal_char (const char c, const char *s_name = 0,
			const char *l_name = 0);
  // Add a wild card to the alphabet with a list of matches.
  void add_wild_card (const char wc, const char *matches,
		      const char *s_name = 0, const char *l_name = 0);
  // Add a match-everything wild card, including characters not yet added.
  void add_all_match (const char wc,
		  const char *s_name=0, const char *l_name=0);  

  void set_char_name(const char ch,
                     const char *short_name, const char *long_name);
  void set_short_name(const char ch, const char *short_name);
  void set_name(const char *nm);
  void reset_name (const char *nm);
  static char *new_copy_string(const char *nm);

  void set_is_nucleic(int v);
  void set_is_rnucleic(int v);
  void set_is_amino(int v);

  static void add_named_alphabet(const Alphabet *alphab);
  static void add_to_alphabet_list(const Alphabet *alphab);

private:
  char *_name;
  char   _chars[Base::_max_chars]; // Legal characters (BaseType to ASCII)
  Base   _translate[Base::_max_chars];       // ASCII to BaseType
  char   _wc_translate[Base::_max_wc][Base::_max_chars]; // [i][j] = WCi matches char j.
  // Null-terminated lists of normal Bases that match normals(1) & wilds(2)
  // Don't forget to add an extra position for NULL character to wild cards!
  Base  _match1[Base::_lowest_wc][2];
  Base  _match2[Base::_max_wc+1][Base::_lowest_wc+1];
  // How many bases match each wildcard.
  int   _nmatch[Base::_max_wc+1];
  int _length;			// Number of non-wildcard, non-variant chars.
  int _wc_length;		// Number of wildcard chars.
  int _case_sensitive;		// 1 if case-sensitive, 0 if not.

  char *_short_names [Base::_max_chars];
  char *_long_names [Base::_max_chars];
  NameToPtr *ShortNameToBaseTable;

  int _is_nucleic;
  int _is_rnucleic;
  int _is_amino;

  // All access is through the alphabet class.

  static NameToPtr *NameToAlphabetTable;
  static const Alphabet **AlphabetList;
  static int AlphabetListLength;
  static int AlphabetListCapacity;
  static int _noisy_convert;
  static const Alphabet *_default_alphabet;

  void init();

protected:
public:
  Alphabet (const char *name = 0,  const char *chars = 0,
            const int case_sensitive = 0);
  virtual ~Alphabet (void);
};

/* ********************************************************************** */
/* ********************************************************************** */
class NucleicAlphabet : public Alphabet 
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(Alphabet::classID());
    }
    static NamedClass *create_NucleicAlphabet(void);

public:
    static IdObject *classID() { return &::NucleicAlphabet::ID; }
    virtual IdObject *type(void) const { return &::NucleicAlphabet::ID; }

public:
  // Really, this should be a BaseType enum (or base), but integers
  // will do and enum is more efficient.
  enum {A = 0, G = 1, C = 2, TU = 3};
  
  // These depend overly on data
  // representation.  The may have to become virtual as wildcards are
  // added.   They currently deal poorly with wildcards.
  virtual Base complement   (const Base b) const;
  virtual int same_group    (const Base b1, const Base b2) const;
  virtual int is_complement (const Base b1, const Base b2) const;
  virtual int is_pyrimidine (const Base b) const;
  virtual int is_purine     (const Base b) const;
  
protected:
  NucleicAlphabet() {}
public:
  NucleicAlphabet (const char *name);
};

class RNAAlphabet : public NucleicAlphabet
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(NucleicAlphabet::classID());
    }
    static NamedClass *create_RNAAlphabet(void);
public:
    static IdObject *classID() { return &::RNAAlphabet::ID; }
    virtual IdObject *type(void) const { return &::RNAAlphabet::ID; }

public:
  enum {U = TU};
protected:
  RNAAlphabet(){}
public:
  RNAAlphabet (const char *nm);
};

class DNAAlphabet : public NucleicAlphabet
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(NucleicAlphabet::classID());
    }
    static NamedClass *create_DNAAlphabet(void);

public:
    static IdObject *classID() { return &::DNAAlphabet::ID; }
    virtual IdObject *type(void) const { return &::DNAAlphabet::ID; }

public:
  enum {T = TU};
protected:
  DNAAlphabet(){}
public:
  DNAAlphabet (const char *nm);
};

class ExtDNAAlphabet : public DNAAlphabet
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(DNAAlphabet::classID());
    }
    static NamedClass *create_ExtDNAAlphabet(void);

public:
    static IdObject *classID() { return &::ExtDNAAlphabet::ID; }
    virtual IdObject *type(void) const { return &::ExtDNAAlphabet::ID; }

public:
  enum {K=Base::_lowest_wc,W=K+1,Y=K+2,M=K+3,R=K+4,S=K+5,V=K+6,B=K+7,
  			D=K+8,H=K+9,N=K+10,X=K+11};
protected:
  ExtDNAAlphabet() {}
public:
  ExtDNAAlphabet (const char *nm);
  virtual Base complement(const Base b) const;
  virtual int is_pyrimidine (const Base b) const;
  virtual int is_purine     (const Base b) const;
  virtual int same_group (const Base b1, const Base b2) const;
  virtual int is_complement (const Base b1, const Base b2) const;
};

/* ********************************************************************** */
class AminoAlphabet : public Alphabet
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(Alphabet::classID());
    }
    static NamedClass *create_AminoAlphabet(void);

public:
    static IdObject *classID() { return &::AminoAlphabet::ID; }
    virtual IdObject *type(void) const { return &::AminoAlphabet::ID; }

public:
  enum {A=0,C=1,D=2,E=3,F=4,G=5,H=6,I=7,K=8,L=9,M=10,
	N=11,P=12,Q=13,R=14,S=15,T=16,V=17,W=18,Y=19};
protected:
  AminoAlphabet() {}
public:
  AminoAlphabet (const char *nm);
};

class ExtAminoAlphabet : public AminoAlphabet
{
//required NamedClass and NamedObject members
private:
    static IdObject ID;
    static void init_is_a(IdObject *self)
    {
        self->add_is_a(AminoAlphabet::classID());
    }
    static NamedClass *create_ExtAminoAlphabet(void);

public:
    static IdObject *classID() { return &::ExtAminoAlphabet::ID; }
    virtual IdObject *type(void) const { return &::ExtAminoAlphabet::ID; }

public:
  enum {B=Base::_lowest_wc,Z=B+1,X=B+2};
protected:
  ExtAminoAlphabet() {}
public:
  ExtAminoAlphabet(const char *nm);
};

/* ********************************************************************** */
// Base inlines  that can't be declared inside the class Base declaration

// wc_match needs full Alphabet class declaration
inline int Base::wc_match  (const Base b2, const Alphabet *a,
		      const VarEnum v) const
{   if (v != NO_VARS)
    {  cerr <<
	"Base::wc_match:  variant characters ignored (not yet implemented).\n";
    }

    if (b2.canon().raw_int() >= a->first_wc() ||
      canon().raw_int() >= a->first_wc())
	    return a->wc_match (*this, b2, v);
    return no_wc_match (b2, NULL_NULL_FALSE, v); // wc matches check null wc
}


// wc_subset needs full Alphabet class declaration
inline int Base::wc_subset  (const Base b2, const Alphabet *a,
		      const VarEnum v) const
{   if (v != NO_VARS)
    {  cerr <<
	"Base::wc_subset:  variant characters ignored (not yet implemented).\n";
    }

    if (b2.canon().raw_int() >= a->first_wc() ||
      canon().raw_int() >= a->first_wc())
	    return a->wc_subset (*this, b2, v);
    return no_wc_match(b2, NULL_NULL_FALSE, v); // check null wc
}


inline ostream&
	operator<< (ostream& outs, const Base b)
	{return (outs << Alphabet::ret_default()->to_char (b));}
inline istream&
	operator>> (istream& ins, Base &b)
	{char c; ins >> c; b = Alphabet::ret_default()->to_base (c); return ins;}

/************************************************************************/
// Alphabet Inlines.
inline char
Alphabet::to_ascii (const Base b1, const VarEnum v = NO_VARS) const
{return (to_char(b1, v));}
inline int
Alphabet::norm_length      (void) const {return _length;}
inline  int
Alphabet::wc_length   (void) const {return _wc_length;}
inline  int
Alphabet::norm_wc_length (void) const {return _wc_length+_length;}
inline  int
Alphabet::max_num_base (void) const {return Base::_var_mask+1;}
inline  int
Alphabet::max_num_var_base (void) const {return Base::_max_chars;}
inline int
Alphabet::first_char  (void) const {return 0;}
inline int
Alphabet::last_char  (void) const {return norm_length()-1;}
inline int
Alphabet::first_wc    (void) const {return Base::_lowest_wc;}
inline int
Alphabet::last_wc    (void) const {return first_wc() + _wc_length-1;}
inline int
Alphabet::first_var   (void) const {assert (0); return Base::_first_var;}
inline int
Alphabet::last_var   (void) const {assert (0); return max_num_var_base()-1;}

inline int
Alphabet::is_case_sensitive() const { return _case_sensitive; }
inline void
Alphabet::set_case_sensitive(int cs) { _case_sensitive = cs; }
inline char
Alphabet::if_case_sensitive_char(const char ch) const
    { return (_case_sensitive || ! isalpha(static_cast<const int>(ch)))
    	? ch : toupper(static_cast<const int>(ch)); }

inline void
Alphabet::silent_convert (int val) {  _noisy_convert = ! val;}

inline void
Alphabet::set_default (const Alphabet& def)
{  _default_alphabet = def.id(); }
inline void
Alphabet::set_default (const char *nm)
{  _default_alphabet = name_to_alphabet(nm); }
inline const Alphabet*
Alphabet::ret_default (void) { return _default_alphabet;}
inline const NameToPtr *
Alphabet::name_to_alphabet_table() { return NameToAlphabetTable; }

inline int
Alphabet::is_nucleic  (void) const {return _is_nucleic;}
inline int
Alphabet::is_rnucleic (void) const {return _is_rnucleic;}
inline int
Alphabet::is_amino    (void) const {return _is_amino;}
inline void
Alphabet::set_is_nucleic(int v) { _is_nucleic = v; }
inline void
Alphabet::set_is_rnucleic(int v) { _is_rnucleic = v; }
inline void
Alphabet::set_is_amino(int v) { _is_amino = v; }

inline int
Alphabet::is_valid (const Base b) const 
	{return _chars[b.raw_int()] != Base::bad_char();}

inline int
Alphabet::is_valid_char(const char ch) const
    { return !_translate[static_cast<int>(if_case_sensitive_char(ch))].is_null()
             || Base::is_null_char(ch); }

inline const Base *
Alphabet::matches (const Base b1, const VarEnum v = NO_VARS) const
{  assert (v == NO_VARS);
   // QUESTION: what should this do with a null or invalid input???
   // ANSWER: return a null-terminated empty list,
   //	which should be the effect of returning _match2
   return (b1.is_null() || b1.canon().is_wild()) ? 
   	   _match2[b1.canon().raw_int() -first_wc()]
	:  _match1[b1.canon().raw_int()];
}

inline int
Alphabet::num_matches (const Base b) const
{   return is_valid(b) ? 
	(b.is_null()?	0
		: (b.canon().is_wild() ?
			_nmatch[b.canon().raw_int() -first_wc()] 
			: 1
		   )
	)
	:0;
}

inline Base
Alphabet::null        (void) const {return Base::null();}
inline int
Alphabet::index       (const Base b) const
{return  b.canon().raw_int() - 
	(b.is_normal() ? 0 : (first_wc() - norm_length()));}
inline int
Alphabet::index_from_char(const char ch) const
{ return index(to_base(ch)); }

inline Base
Alphabet::unindex     (const int i) const
{return Base::base_int ( (i < norm_length()) ? i : (i + first_wc() - norm_length()));}
inline char
Alphabet::unindex_to_char(const int i) const
{ return to_char(unindex(i)); }

inline const char *
Alphabet::abbrev (const Base b) const 
{   int i = b.raw_int();
    assert(i>=0 && i<Base::_max_chars);
    return _short_names[i];
}

inline const char *
Alphabet::full_name (const Base b) const 
{   int i = b.raw_int();
    assert(i>=0 && i<Base::_max_chars);
    return _long_names[i];
}

inline void
Alphabet::reset_name (const char * nm) { set_name(nm); }

inline char *Alphabet::new_copy_string(const char *str)
{
    if (str == NULL) return NULL;

    int new_len = strlen(str);
    char *new_str = new char[new_len + 1];
    strcpy(new_str, str);
    return new_str;
}

inline Base
NucleicAlphabet::complement (const Base b) const
{   return b.is_null() ? b : Base::base_int ((b.raw_int())^3);}
inline int
NucleicAlphabet::is_pyrimidine (const Base b) const
{return ((b.canon()).raw_int()==A || (b.canon()).raw_int()==G);}
inline int
NucleicAlphabet::is_purine     (const Base b) const
{return ((b.canon()).raw_int()==C || (b.canon()).raw_int()==TU);}
inline int
NucleicAlphabet::same_group (const Base b1, const Base b2) const
{return !b1.is_null() && !b2.is_null() && ( (b1.raw_int() & 0x2) == (b2.raw_int() & 0x2));} // Match in MSB.
inline int
NucleicAlphabet::is_complement (const Base b1, const Base b2) const
{return (b1.canon()).raw_int() == (complement (b2.canon()).raw_int());}


/* ********************************************************************** */
/* ********************************************************************** */
// It may be more efficient at some point to have the SeqList store
// a->first_wc() and pass it as a paramater rather than pointer
// chasing.

// CHANGE LOG:
// 20 August 2003 Kevin Karplus
//	Fixed is_wild() to check for is_null().
// 27 Feb 2004 Kevin Karplus
//	Converted old-style casts.
// 3 March 2004 Kevin Karplus
//	Added assertions to abbrev() and full_name() members.
// 15 March 2004 Kevin Karplus
//	Got rid of old-style casts hidden in isalpha and toupper
// 18 May 2004 Kevin Karplus
//	Fixed num_matches to recognize that is_null no longer
//	implies is_wild.
// 18 May 2004 Kevin Karplus
//	Fixed matches to treat is_null like it does is_wild.
// Fri Mar 25 17:58:17 PST 2005 Kevin Karplus
//	Added static_cast for is_valid_chat call of if_case_sensitive_char
// Sat 12 May 2007 John Archie
//	Made is_all_match() public so an all_match wildcard may be determined
//      from outside of the alphabet class
// Sat 19 May 2007 John Archie
//      Made the read functions protected so that they could be called from
//      a derived class: ReadName, ReadComment, VerifyClassName,
//	ReadNormalChars, ReadAlias, ReadWildcard, ReadAllMatch, ReadCharName,
//      ReadIsNucleic, ReadIsRNucleic, and ReadIsAmino.  Also, made print
//      commands protected for same reason: print_normal_chars,
//	print_wildcards, print_alias_chars, print_char_names, and
//	print_is_values).
//	These changes were necessary to allow significant code reuse in
//      undertaker's AngleVectorAlphabet class.
//
// Fri Apr  3 09:45:26 PDT 2009 Kevin Karplus
//	Added word "void" to void constructor.
// Thu Aug 13 17:07:23 PDT 2009 Kevin Karplus
//	Added first_all_match_wildcard() [from AngleVectorAlphabet]
#endif // ! _ALPHABET_H

