// Alphabet.C
// Revise 21 July 1994 rph@ce.ucsc
// 3 Nov 93 Richard Hughey
//
// Storage and functions for Alphabet.h (classes Base, AlphClass and
// descendents, and Alphabet).

//Change log:
//12 Oct 96, changed uses of String::elem() to String::operator[]
//           so that it would work with both the cxx lib and g++ lib
//           since g++ String has now been fixed for const Strings

// 21 Jan 1999, fixed 3-letter code for isoleucine
// 23 Feb 1999, Changed uses of String class to char*.
//
//              Added unindex_to_char(), which takes the integer
//              index of a base and returns its character represenation.
//
//              Changed Alphabet to a subclass of NamedClass and NamedObject.
//              NamedClass-style I/O is now supported.
//
//              Added automatic registering of Alphabets in a lookup
//              table whenever set_name() is called with a non-NULL string.
//              Alphabets can be looked up using Alphabet::name_to_alphabet()
//
//              Separated Alph class into its own header file and cc file
//              Old Alph::describe() is now Alphabet::describe().
//
//              Added new to_base() which takes a short name as an argument
//              and returns the corresponding base.
//
//  23 July 1999 Kevin Karplus
//	increased size of _match2 so that null_base can be used in "matches"
//
// 16 August 2001 Kevin Karplus
//	Modified load_alphabet_file to use gzifstream and Filenames::
//
// 17 July 2003 George Shackelford
//	Replaced 'ipfx(1)' with 'good()' for gcc 3.x
//
#include "Alphabet.h"
/* #include "zfstream/zfstream.h"  ggs 8/13/03 */
#include "zfstream.h"
#include "Filenames.h"
#include "SkipSeparators.h"

#include <ctype.h>
#include <stdlib.h>
#include <fstream>
/* ********************************************************************** */

//
// class NameBasePair is used for entries into the lookup table
// for translating short names of bases to the bases.
//
class NameBasePair : public NamedObject
{
private:
    Base b;
public:
    NameBasePair(const char *nm, Base b_val)
    {
        set_name(nm);
        this->b = b_val;
    }

    Base get_base() { return b; }
};


/* ********************************************************************** */
int Alphabet::_noisy_convert = 1;
const Alphabet * Alphabet::_default_alphabet = NULL;

/* ********************************************************************** */
/* ********************************************************************** */

static const char *return_empty_if_null(const char str[])
{
    return (str == NULL) ? "" : str;
}

/* ********************************************************************** */
void
Alphabet::describe (ostream& s) const
{
    int oldn = _noisy_convert;
    _noisy_convert = 0;
  
    //const Alphabet *a = this;

    assert (id() == name_to_alphabet(name())->id());
        s <<"The " << name()
          << " alphabet is " <<
          (is_amino() ? "based on amino acids" :
           is_rnucleic() ? "r-nucleic" :
           is_nucleic()  ? "nucleic" :
           "of unknown type") << " and includes the characters:\n    \t"; 

    int charnum = 0;
    int i;
    for (i = first_char() ; i < max_num_base() ; i++)
    {
        if (is_valid (Base::base_int(i)))
        {
	    const char *abr = return_empty_if_null(abbrev(Base::base_int(i)));
	    const char *lnam = return_empty_if_null
                                   (full_name (Base::base_int(i)));
            if (strcmp(abr, lnam) == 0) lnam = "";
            if (4 + strlen(abr) + strlen(lnam) + charnum > 65)
            {
	       s << "\n    \t";
               charnum = 0;
            }
	
            // check conversions to-from int and to-from char
            Base b = Base::base_int(i);
            assert (i == b.raw_int());
            char c = to_char(b);
            assert(to_base(c).raw_int() == i);
            s << c;
	
            charnum+= 2;
            if (strlen(abr) > 0)
            {
                if (strlen(lnam) == 0)
                {
                    s << "(" << abr << ")";
	            charnum += strlen(abr) + 2;
	        }
                else
                {
	            charnum += strlen(abr) + strlen(lnam) + 4;
	            s << "(" <<  abr << "," <<  lnam << ") ";
	        }
            }
            s << " ";
        }
    }

    s << endl;
    s << "    " << name() << " has " << norm_length()
      << " normal characters, " << wc_length()
      << " wildcards, and the null character.\n";
    for (i = 'A' ; i <= 'Z' ; i++)
    {
        char c=0; Base b;
        if ((b = to_base (i)).raw_int() != Base::null().raw_int() 
             && (c = to_char(b)) != i)
        s << "    The letter " << static_cast<char>(i) << " is an alias for "
          << c << ".\n";
    }
    if (wc_length() > 0)
    {
        s << "    Wildcard listing: all characters match and lowercase are not subsets\n";
        for (i = first_wc() ; i <= last_wc() ; i++)
        {
            Base wcb = Base::base_int(i);
	    assert( is_valid(wcb));
	    assert( wcb.canon().is_wild());

	    s << "    Wildcard " << i << " (" 
              << to_char (Base::base_int(i)) <<   ") : ";

            for (int j = 0 ; j < max_num_base() ; j++)
            {
                if (is_valid(Base::base_int(j)))
                {
                    if (wc_match (Base::base_int(i), Base::base_int(j)))
                    {
                        if (wc_subset(Base::base_int(i),Base::base_int(j))) 
                            s << to_char (Base::base_int(j));
                        else
                            s << static_cast<char>(tolower (to_char (Base::base_int(j))));
                    }
                }
            }

            // Print out the match function for the wildcards
            s << "\t matches: ";
            int NumM= 0;
            for (const Base *bp = matches(Base::base_int(i));
                 !bp->is_null() ; bp++)
            {
                s << to_char (*bp);
                NumM++;
            }
	    s << " (" << num_matches(Base::base_int(i))
	    	<< " matches)" << flush;
            assert(NumM==num_matches(Base::base_int(i)));
            if (is_nucleic()) 
            {
                s << "\t complement: ";
                s << to_char(dynamic_cast<const NucleicAlphabet*>(this)
                             ->complement(Base::base_int(i)));
            }
            s << endl;
        }
    }
    // Confirm the match function for the non-wildcards
    for (i = 0 ; i < first_wc() ; i++) 
    {
        Base b = Base::base_int(i);
        const Base *bp = matches(b);
//        cerr << "is_valid(b)=" << is_valid(b)
//		<< " b.canon().raw_int()=" << b.canon().raw_int()
//		<< " _match1[b.canon().raw_int()][0]=" << _match1[b.canon().raw_int()][0]
//		<< " bp->is_null()=" << bp->is_null()
//		<< " bp[0].raw_int()=" << bp[0].raw_int()
//		<< " b.raw_int()=" << b.raw_int()
//		<< " bp[1].is_null()=" << bp[1].is_null()
//		<< "\n" << flush;
	assert ( (! is_valid(b) && bp->is_null())
                 || (bp[0].raw_int() == b.raw_int()) && (bp[1].is_null()));

    }
    s << endl;

    _noisy_convert = oldn;
}

/* ********************************************************************** */
void
Alphabet::describe_all (ostream& s)
{
    int oldn = _noisy_convert;
    _noisy_convert = 0;
  
    s << "There are currently " << num_alphabets()
      << " alphabets available as follows: \n";
    for (int k = 0 ; k < num_alphabets() ; k++)
    {
        const Alphabet *a = alphabet_list(k);

        assert (a->id() == name_to_alphabet(a->name())->id());
        s << k << ".  ";
        a->describe(s);
    }
    _noisy_convert = oldn;
}

/************************************************************************/
void 
Base::limits(ostream& s) {
  s << "The current Base implementation has the following limits:";
  s << "\n   Bytes of storage per base:         " << sizeof(Base);
  s << "\n   Number of non-wildcard characters: " << static_cast<int>(_lowest_wc);
  s << "\n   Number of wildcard characters:     " << _null_base - _lowest_wc;
  s << "\n   The null character is number:      " << null().raw_int() << "   `"
    << Base::null_char() << "'";
  s << "\n   Character variants range from " << static_cast<int>(_first_var) << " to "
    << _max_chars-1;
  s << "\n   Illegal characters are represented as   `" <<
    Base::bad_char() << "'\n";
}

/*
 * Alphabet methods
 *
 */
NamedClass *Alphabet::create_Alphabet(void) { return new Alphabet; }
IdObject Alphabet::ID("Alphabet",
                      Alphabet::create_Alphabet,
                      Alphabet::init_is_a,
                      "the root Alphabet class");

NameToPtr *Alphabet::CommandTable = 0;
NameToPtr *Alphabet::NameToAlphabetTable = NULL;
const Alphabet **Alphabet::AlphabetList = NULL;
int Alphabet::AlphabetListLength = 0;
int Alphabet::AlphabetListCapacity = 0;

NamedClass *NucleicAlphabet::create_NucleicAlphabet(void)
{ return new NucleicAlphabet; }
IdObject NucleicAlphabet::ID("NucleicAlphabet",
                             NucleicAlphabet::create_NucleicAlphabet,
                             NucleicAlphabet::init_is_a,
                             "Nucleic Alphabet class");

NamedClass *RNAAlphabet::create_RNAAlphabet(void) { return new RNAAlphabet; }
IdObject RNAAlphabet::ID("RNAAlphabet",
                         RNAAlphabet::create_RNAAlphabet,
                         RNAAlphabet::init_is_a,
                         "RNA Alphabet class");

NamedClass *DNAAlphabet::create_DNAAlphabet(void) { return new DNAAlphabet; }
IdObject DNAAlphabet::ID("DNAAlphabet",
                         DNAAlphabet::create_DNAAlphabet,
                         DNAAlphabet::init_is_a,
                         "DNA Alphabet class");

NamedClass *ExtDNAAlphabet::create_ExtDNAAlphabet(void)
{ return new ExtDNAAlphabet; }
IdObject ExtDNAAlphabet::ID("ExtDNAAlphabet",
                            ExtDNAAlphabet::create_ExtDNAAlphabet,
                            ExtDNAAlphabet::init_is_a,
                            "Extended DNA Alphabet class");

NamedClass *AminoAlphabet::create_AminoAlphabet(void)
{ return new AminoAlphabet; }
IdObject AminoAlphabet::ID("AminoAlphabet",
                           AminoAlphabet::create_AminoAlphabet,
                           AminoAlphabet::init_is_a,
                           "Amino Acid Alphabet class");

NamedClass *ExtAminoAlphabet::create_ExtAminoAlphabet(void)
{ return new ExtAminoAlphabet; }
IdObject ExtAminoAlphabet::ID("ExtAminoAlphabet",
                              ExtAminoAlphabet::create_ExtAminoAlphabet,
                              ExtAminoAlphabet::init_is_a,
                              "Extended Amino Acid Alphabet class");


void Alphabet::init()
{
    _name = NULL;
    _length = 0;
    _wc_length = 0;
    _case_sensitive = 0;

    _is_nucleic = 0;
    _is_rnucleic = 0;
    _is_amino = 0;

    int i;
    for (i = 0; i < Base::_max_chars; i++)
    {
        _short_names[i] = NULL;
        _long_names[i] = NULL;
    }

    for (i = 0 ; i < Base::_max_chars ; i++) {
        _chars[i] = Base::bad_char();
        _translate[i] = Base::null();
        assert(_translate[i].is_null());
    }
    assert( Base::_lowest_wc >0 );
    for (i = 0 ; i < Base::_lowest_wc ; i++)
    {    _match1[i][0] = Base::null(); 
    	 _match1[i][1] = Base::null();
    }
    assert(_match1[_length][0].is_null());

    for (i = 0 ; i <= Base::_max_wc ; i++) {
        _match2[i][0] = Base::null();
        _nmatch[i] = 0;
    }

    _chars[Base::null().raw_int()] = Base::null_char();
    ShortNameToBaseTable = new NameToPtr(257);
    ShortNameToBaseTable->ignore_case();
   
    assert(_length==0);
    assert(_match1[_length][0].is_null());
}

/* ********************************************************************** */
Alphabet::Alphabet (const char *nme, const char *chars,
		    const int case_sensitive) :
_wc_length (0), _case_sensitive (case_sensitive)
{  // cerr << "Constructing named alphabet " << nme << "\n" << flush;
  init();
  set_name(nme);
  set_case_sensitive(case_sensitive);
  if (chars != NULL)
  {
      int chars_len = strlen(chars);
      for (int j = 0 ; j < chars_len ; j++)
          add_normal_char (chars[j]);
  }
}

const Alphabet *const
Alphabet::id (void) const
{
  return this;
}

Alphabet::~Alphabet()
{
    delete [] _name;

    int i;
    for (i = 0; i < Base::_max_chars; i++)
    {
        delete [] _short_names[i];
        delete [] _long_names[i];
    }

    delete ShortNameToBaseTable;
}


Base
Alphabet::to_base (const char c) const
{						  
  Base rval = _translate[static_cast<int>(if_case_sensitive_char(c))];
  if (Alphabet::_noisy_convert && rval.is_null() && !Base::is_null_char(c))
    cerr << "Bad base " <<  c << " for alphabet " << name() << "\n";
//  assert (! rval.is_null());   // to_base illegal translation.
  return rval;
}

Base
Alphabet::to_base  (const char *short_name) const
{
    if (short_name == NULL)
    {
        return Base::null();
    }

    NameBasePair *nbp = dynamic_cast<NameBasePair *>(
                        ShortNameToBaseTable->FindOldName(short_name,
                                                          ZeroIfNew));

    if (nbp != NULL)
    {
        return nbp->get_base();
    }
    else
    {   if (Alphabet::_noisy_convert)
    	{   cerr << "Bad short name: " << short_name
            	<< " for alphabet: " << name()
             	<< "\n";
        }
	return Base::null();
    }
}

void
Alphabet::add_alias (const char newchar, const char alias)
{
  assert (_translate[static_cast<int>(newchar)].is_null());
  assert (!_translate[static_cast<int>(alias)].is_null());
  _translate[static_cast<int>(newchar)] = _translate[static_cast<int>(alias)];
}
void
Alphabet::add_normal_char (const char ch, const char *s_name,
			   const char *f_name)
{
  assert (_wc_length == 0);// Can't add normal characters after wildcards.

  char c = if_case_sensitive_char(ch);
  if (! _translate[static_cast<int>(c)].is_null() )
  {    cerr << "Error: character " << ch << " converted to " << c 
  	<< " which already had a definition in Alphabet " << name()
	<< "\n" << flush;
  }
  assert (_translate[static_cast<int>(c)].is_null());
  _chars[_length] = c;
  _translate[static_cast<int>(c)] = Base::base_int(_length);
//  cerr << " name=" << name()
//  	<< " _length=" << _length  
//  	<< " _match1[_length][0].raw_int()=" << _match1[_length][0].raw_int()
//	<< "\n" << flush;
  assert(_match1[_length][0].is_null());
  _match1[_length][0] = Base::base_int(_length);  // Normals match only themselves.

  set_short_name(ch, s_name);
  _long_names [_length] = new_copy_string(f_name);

  if (_length > Base::_lowest_wc) {
    Base::limits (cerr);
    assert (_length <= Base::_lowest_wc);
  }

  _length++;
}
void
Alphabet::add_wild_card (const char wch, const char *matches_what,
			 const char *s_name, const char *f_name)
{
  char wc = if_case_sensitive_char(wch);
  assert (_translate[static_cast<int>(wc)].is_null());// Not already in use.

  // Add the character to the character list
  _chars[static_cast<int>(first_wc() + _wc_length)] = wc;
  // Add its translation
  _translate[static_cast<int>(wc)] = Base::base_int(first_wc() + _wc_length);

  set_short_name(wch, s_name);
  _long_names [static_cast<int>(first_wc()  + _wc_length)] = new_copy_string(f_name);

  int i,j;
  for (i = 0 ; i < Base::_max_chars ; i++) {
    _wc_translate[_wc_length][i] = 0;
  }
  int matches_len = strlen(matches_what);
  for (i = 0, j = 0; i < matches_len ; i++) {
    char matches_char = matches_what[i];
    Base b = to_base (matches_char);
    assert (is_valid (b));
    // assert valid character in wc_list.  .
    _wc_translate[_wc_length][b.raw_int()] = 1;
    if (b.is_normal()) {
      _match2[_wc_length][j++] = b;
      _nmatch[_wc_length]++;
    }
  }
  assert(_nmatch[_wc_length] > 0);
  _match2[_wc_length][j] = Base::null();
  
  // Number of wildcards.
  _wc_length++;
  if (_wc_length > Base::_null_base) {
    Base::limits (cerr);
    assert (_wc_length <= Base::_null_base);
  }
}
void
Alphabet::add_all_match (const char wch,
			 const char *s_name, const char *f_name)
{
  char wc = if_case_sensitive_char(wch);
  assert (_translate[static_cast<int>(wc)].is_null());	  // Not already in use.

  _translate[static_cast<int>(wc)] = Base::base_int(first_wc() + _wc_length);
  _chars[first_wc() + _wc_length] = wc;

  set_short_name(wch, s_name);
  _long_names [first_wc()  + _wc_length] = new_copy_string(f_name);
  
  // Match this and all future characters and wildcards.
  // Note: actually will only match all current normal characters
  // because add_normal_char does not allow any more normal chars
  // to be added after a wildcard is added.
  int i;
  for (i = 0 ; i < Base::_max_chars ; i++) {
    _wc_translate[_wc_length][i] = 1;
  }
  int j;
  for (j=0, i = 0 ; Base::base_int (i).is_normal() ; i++)
    if (is_valid(Base::base_int(i))) {
      _match2[_wc_length][j++] = Base::base_int(i);
      _nmatch[_wc_length]++;      
    }
  _match2[_wc_length][j] = Base::null();  
  _wc_translate[_wc_length][Base::null().raw_int()] = 0;
  _wc_length++;
  if (_wc_length > Base::_null_base) {
    Base::limits (cerr);
    assert (_wc_length <= Base::_null_base);
  }
}
Base
Alphabet::valid_or_null (const Base b) const
{
  assert (b.canon()._b == b._b);   // No variants yet!
  if (is_valid(b))
    return b;
  else
    return null();
}
/************************************************************************/
int
Alphabet::wc_match (const Base b1, const Base b2,
		    const VarEnum v) const
{
  if (b1.is_null() || b2.is_null())
    return 0;
  if (! (b1.is_wild() || b2.is_wild()))
    return b1.no_wc_match (b2, NULL_NULL_FALSE, v);

  Base wild = (b1.is_wild() ? b1 : b2).canon();
  Base base = (b1.is_wild() ? b2 : b1);
  if (v == NO_VARS) {
    wild = wild.canon();
    base = base.canon();
  }

  return (_wc_translate [wild.canon().raw_int()-first_wc()][base.raw_int()]
	|| (base.is_wild() &&
	    _wc_translate[base.canon().raw_int()-first_wc()][wild.raw_int()]));

}
/************************************************************************/
int
Alphabet::wc_subset (const Base b1, const Base b2,
		     const VarEnum v) const
{
  if (b1.is_null() || b2.is_null())
    return 0;
  if (! b1.is_wild())
    return b1.no_wc_match (b2, NULL_NULL_FALSE, v);
  
  return static_cast<int> (_wc_translate [b1.canon().raw_int()-first_wc()]
		[v == NO_VARS ? b2.canon().raw_int() : b2.raw_int()]);

}

void
Alphabet::set_char_name(const char ch,
			const char *short_name, const char *long_name)
{
  assert (!_translate[static_cast<int>(ch)].is_null());// Should be already in use.

  int raw_i = _translate[static_cast<int>(ch)].raw_int();
  delete [] _long_names[raw_i];
  _long_names [raw_i] = new_copy_string(long_name);

  set_short_name(ch, short_name);
}

void
Alphabet::set_short_name(const char ch, const char *short_name)
{
  assert (!_translate[static_cast<int>(ch)].is_null());// Should be already in use.

  int raw_i = _translate[static_cast<int>(ch)].raw_int();
  delete [] _short_names[raw_i];

  if (short_name != NULL)
  {
      _short_names [raw_i] = new_copy_string(short_name);
      ShortNameToBaseTable->AddName(new NameBasePair(short_name,
                                                      _translate[static_cast<int>(ch)]),
                                     ErrorIfOld);
  }
}

/*
 * Alphabet NamedClass I/O
 *
 */
enum { TraceRead = 0 };
enum { BUF_LEN = 256 };

void Alphabet::add_named_alphabet(const Alphabet *alphab)
{
    if (NameToAlphabetTable == NULL)
    {
        NameToAlphabetTable = new NameToPtr(113);
    }

    NameToAlphabetTable->AddName( const_cast<Alphabet*>(alphab), ErrorIfOld);
}

void Alphabet::add_to_alphabet_list(const Alphabet *alphab)
{
    if (AlphabetListLength == AlphabetListCapacity)
    {
        int new_capacity = AlphabetListCapacity == 0
                           ? 1 : 2*AlphabetListCapacity;
        const Alphabet **new_list = new const Alphabet*[new_capacity];
        int i;

        for (i = 0; i < AlphabetListLength; i++)
        {
            new_list[i] = AlphabetList[i];
        }

        delete [] AlphabetList;
        AlphabetList = new_list;
        AlphabetListCapacity = new_capacity;
    }

    AlphabetList[AlphabetListLength] = alphab;
    AlphabetListLength++;
}

const Alphabet *Alphabet::name_to_alphabet(const char *nm, OptionIfNew ifnew)
{
    return dynamic_cast<Alphabet *>(NameToAlphabetTable->FindOldName(nm, ifnew));
}

const Alphabet *Alphabet::alphabet_list(int i)
{
    assert(0 <= i && i < num_alphabets());

    return AlphabetList[i];
}

int Alphabet::num_alphabets()
{
    return AlphabetListLength;
}

void Alphabet::set_name(const char *nm)
{
    const char *old_nm = name();

    if (old_nm != NULL)
    {
        NameToAlphabetTable->DeleteName(name(), ErrorIfNew);
    }

    NamedObject::set_name(nm);

    if (name() != NULL)
    {
        add_named_alphabet(this);

        if (old_nm == NULL)
        {
            add_to_alphabet_list(this);
        }
    }
}

int Alphabet::read_knowing_type(istream &in)
{
    if (!CommandTable)
        init_command_table();

    char word[BUF_LEN];
    int non_end_class = 1;

    while (in.good() && non_end_class)
    {
        get_word(in, word, '=');
        AlphabetInputCommand *comm = dynamic_cast<AlphabetInputCommand *>
                                (CommandTable->FindOldName(word, ZeroIfNew));

        if (TraceRead)
	{
            cerr << "Alphabet::read_knowing_type: " << word << "\n" << flush;
        }
        if (comm != 0)
        {
            non_end_class = comm->execute(in, this);
        }
        else
        {
            cerr << "Unknown Alphabet Command: "
                 << word 
                 << "\n" << flush;
        }
    }

    if (non_end_class)
    {
        cerr << "Warning: Alphabet::read_knowing_type(): expected closing"
             << " EndClassName command\n";
    }

    if (name() == NULL)
    {
        cerr << "Warning: Alphabet::read_knowing_type(): alphabet has"
                " no name, not registered in lookup table.\n";
    }

    return !non_end_class;
}

void Alphabet::print_normal_chars(ostream &out) const
{
    int bi;

    out << "NormalChars = ";
    for (bi = 0; bi < norm_length(); bi++)
    {
        char ch = to_char(unindex(bi));

        out << ch;
    }

    out << "\n";
}

void Alphabet::print_alias_chars(ostream &out) const
{
    int ci;

    for (ci = 'A'; ci < 'Z'; ci++)
    {
        Base b = _translate[ci];

        if (!b.is_null())
        {
            char back_ch = to_char(b);

            if ( ci != static_cast<int>(back_ch) )
            {
                out << "Alias = " << static_cast<char>(ci) << " " << back_ch << "\n";
            }
        }
    }
}


int Alphabet::is_all_match(Base wcb) const
{
    int saved_noisy_convert = _noisy_convert;
    _noisy_convert = 0;

    int all_match = 1;
    int raw_i;

    for (raw_i = 0; raw_i < max_num_base(); raw_i++)
    {
        Base b = Base::base_int(raw_i);

        if (is_valid(b) && !b.is_null())
        {
            all_match = all_match && (wc_match(wcb, b) && wc_subset(wcb, b));

            /*******************************************
            if (!all_match)
            {
                cerr << "Alphabet: " << name()
                     << " wcb: " << to_char(wcb)
                     << " fails all_match on char: " << to_char (b)
                     << "\n";
            }
            *****************************************/
        }
    }

    _noisy_convert = saved_noisy_convert;

    return all_match;
}

void Alphabet::print_wildcards(ostream &out) const
{
    int saved_noisy_convert = _noisy_convert;
    _noisy_convert = 0;

    int wci;

    for (wci = norm_length(); wci < norm_wc_length(); wci++)
    {
        Base wcb = unindex(wci);
        char wch = to_char(wcb);

        assert( is_valid(wcb));
	assert( wcb.canon().is_wild());
	assert( index(wcb.canon()) == wci);
	
	if (is_all_match(wcb))
        {
            out << "AllMatch = " << wch;
        }
        else
        {
            out << "Wildcard = " << wch << " ";

            int raw_i;
            for (raw_i = 0; raw_i < max_num_base(); raw_i++)
	    {
                Base b = Base::base_int(raw_i);

                if (is_valid(b))
                {
                    if (wc_match(wcb, b) && wc_subset(wcb, b))
                    {
                        out << to_char(b);
                    }
                }
            }
        }

        out << "\n";
    }

    _noisy_convert = saved_noisy_convert;
}

void Alphabet::print_char_names(ostream &out) const
{
    int bi;

    for (bi = 0; bi < norm_wc_length(); bi++)
    {
        Base b = unindex(bi);
        char ch = to_char(b);
        const char *short_name = abbrev(b);
        const char *long_name = full_name(b);

        if (short_name != NULL
            && strlen(short_name) > 0
            && long_name != NULL
            && strlen(long_name) > 0)
        {
            out << "CharName = " << ch
                << " " << short_name
                << " " << long_name
                << "\n";
        }
    }
}

void Alphabet::print_is_values(ostream &out) const
{
    if (is_nucleic())
    {
        out << "IsNucleic = " << is_nucleic() << "\n";
    }

    if (is_rnucleic())
    {
        out << "IsRNucleic = " << is_rnucleic() << "\n";
    }

    if (is_amino())
    {
        out << "IsAmino = " << is_amino() << "\n";
    }
}

void Alphabet::write_knowing_type (ostream &out) const
{
    out << "\nName = " << name() << "\n";

    print_is_values(out);
    print_normal_chars(out);
    print_wildcards(out);
    print_alias_chars(out);

    print_char_names(out);

    out << "\n";
}

int Alphabet::ReadName(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char word[BUF_LEN];
    get_word(in, word);
    change->set_name(word);
    return 1;
}

int Alphabet::ReadComment(istream &in, 
                Alphabet *change,
                AlphabetInputCommand* self)
{   
    SkipSeparators(in, 1, '\n');
    return 1;
}

int Alphabet::VerifyClassName(istream &in, 
                Alphabet *change,
                AlphabetInputCommand* self)
{
    char word[BUF_LEN];
    get_word(in, word);
    const IdObject *end_id = IdObject::id(word);
    if (end_id != change->type())
    {    cerr << "Warning: " << self->name() << word << " doesn't match "
                << change->type()->name() << "\n" << flush;
    }
    // continue if "ClassName", stop if "EndClassName"
    return EqualStrings(self->name(), "ClassName", 1);
}

int Alphabet::ReadNormalChars(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char norm_chars[BUF_LEN];
    get_word(in, norm_chars);

    int norm_chars_len = strlen(norm_chars);
    for (int j = 0 ; j < norm_chars_len ; j++)
        change->add_normal_char(norm_chars[j]);

    return 1;
}

int Alphabet::ReadAlias(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char new_chars[BUF_LEN];
    get_word(in, new_chars);

    char old_char[BUF_LEN];
    get_word(in, old_char);

    int new_chars_len = strlen(new_chars);
    int i;
    for (i = 0; i < new_chars_len; i++)
    {
        change->add_alias(new_chars[i], old_char[0]);
    }

    return 1;
}

int Alphabet::ReadWildcard(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char wildcard[BUF_LEN];
    get_word(in, wildcard);

    char norm_chars[BUF_LEN];
    get_word(in, norm_chars);

    change->add_wild_card(wildcard[0], norm_chars);

    return 1;
}

int Alphabet::ReadAllMatch(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char wildcard[BUF_LEN];
    get_word(in, wildcard);

    change->add_all_match(wildcard[0]);

    return 1;
}

int Alphabet::ReadCharName(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char ch[BUF_LEN];
    get_word(in, ch);

    char short_name[BUF_LEN];
    get_word(in, short_name);

    char long_name[BUF_LEN];
    get_word(in, long_name);

    // Check first that the character in the CharName command is
    // defined as a base for this Alphabet
    // Second, check that this is not a duplicate CharName for the
    // base
    if (change->to_base(ch[0]).is_null())
    {
        const char *name_of_change = "Unnamed Alphabet";
        if (change->name() != NULL)
        {
            name_of_change = change->name();
        }
        cerr << "Error: input file for "
             << name_of_change
             << " has a CharName command for character " << ch[0]
             << " which in not defined in the Alphabet.\n";

        abort();
    }
    else if (change->abbrev(change->to_base(ch[0])) != NULL)
    {
        const char *name_of_change = "Unnamed Alphabet";
        if (change->name() != NULL)
        {
            name_of_change = change->name();
        }

        cerr << "Warning: input file for "
             << name_of_change
             << " has a duplicate CharName command for character "
             << ch[0]
             << ", command ignored."
             << endl
             << "    old short name: "
             << change->abbrev(change->to_base(ch[0])) 
             << "    old long name: "
             << change->full_name(change->to_base(ch[0]))
             << endl
             << "    duplicate short name: "
             << short_name 
             << "    duplicate long name: "
             << long_name
             << endl;
    }
    else /* No error found to this point */
    {
        change->set_char_name(ch[0], short_name, long_name);
    }

    return 1;
}

int Alphabet::ReadIsNucleic(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char ch[BUF_LEN];
    get_word(in, ch);

    int val = atol(ch);
    change->set_is_nucleic(val);

    return 1;
}

int Alphabet::ReadIsRNucleic(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char ch[BUF_LEN];
    get_word(in, ch);

    int val = atol(ch);
    change->set_is_rnucleic(val);

    return 1;
}

int Alphabet::ReadIsAmino(istream &in, Alphabet *change, AlphabetInputCommand *self)
{
    char ch[BUF_LEN];
    get_word(in, ch);

    int val = atol(ch);
    change->set_is_amino(val);

    return 1;
}

void Alphabet::init_command_table(void)
{
    assert(CommandTable == 0);
    CommandTable = new NameToPtr(23);
    CommandTable->ignore_case();

    CommandTable->AddName(new AlphabetInputCommand("Comment", ReadComment));
    CommandTable->AddName(new AlphabetInputCommand("ClassName", VerifyClassName));
    CommandTable->AddName(new AlphabetInputCommand("EndClassName", VerifyClassName));
    CommandTable->AddName(new AlphabetInputCommand("Name", ReadName));
    CommandTable->AddName(new AlphabetInputCommand("NormalChars",
                                                    ReadNormalChars));
    CommandTable->AddName(new AlphabetInputCommand("Alias",
                                                   ReadAlias));
    CommandTable->AddName(new AlphabetInputCommand("Wildcard",
                                                   ReadWildcard));
    CommandTable->AddName(new AlphabetInputCommand("AllMatch",
                                                   ReadAllMatch));
    CommandTable->AddName(new AlphabetInputCommand("CharName",
                                                   ReadCharName));
    CommandTable->AddName(new AlphabetInputCommand("IsNucleic",
                                                   ReadIsNucleic));
    CommandTable->AddName(new AlphabetInputCommand("IsRNucleic",
                                                   ReadIsRNucleic));
    CommandTable->AddName(new AlphabetInputCommand("IsAmino",
                                                   ReadIsAmino));
}

int Alphabet::load_alphabet_file(const char *filename)
{
    vector<NamedClass*> alphabets;
    int num_alphabets_read = NamedClass::read_many_new(alphabets,filename,Alphabet::classID());

//    cerr << "Alphabet::load_alphabet_file: num_alphabets_read: "
//         << num_alphabets_read
//         << "\n";

    // The vector is automatically discarded, but it does not own the
    // alphabets read in, so they remain, accessible through the name_to_alphabet_table.

    return num_alphabets_read;
}


/* ********************************************************************** */
NucleicAlphabet::NucleicAlphabet (const char *nme) :
Alphabet (nme)
{
  set_is_nucleic(1);
  add_normal_char ('A', "Adenine", "Adenine");
  add_normal_char ('G', "Guanine", "Guanine");
  add_normal_char ('C', "Cytosine", "Cytosine");

  if (nme == "Nucleic") {
    add_normal_char ('T', "Thymine", "Thymine");
    add_alias ('U', 'T');
  }
}
RNAAlphabet::RNAAlphabet (const char *nm) :
NucleicAlphabet (nm)
{
  set_is_rnucleic(1);
  add_normal_char ('U', "Uracil", "Uracil");
  add_alias ('T', 'U');
}
DNAAlphabet::DNAAlphabet (const char *nm) :
NucleicAlphabet (nm)
{
  add_normal_char ('T', "Thymine", "Thymine");
  add_alias ('U', 'T');
}
AminoAlphabet::AminoAlphabet (const char *nm)  :
Alphabet (nm)
{
  set_is_amino(1);
  add_normal_char ('A', "Ala","Alanine");
  add_normal_char ('C', "Cys","Cysteine");
  add_normal_char ('D', "Asp","Aspartic_acid");
  add_normal_char ('E', "Glu","Glutamic_acid");
  add_normal_char ('F', "Phe","Phenylalanine");
  add_normal_char ('G', "Gly","Glycine");
  add_normal_char ('H', "His","Histidine");
  add_normal_char ('I', "Ile","Isoleucine");
  add_normal_char ('K', "Lys","Lysine");
  add_normal_char ('L', "Leu","Leucine");
  add_normal_char ('M', "Met","Methionine");
  add_normal_char ('N', "Asn","Asparagine");
  add_normal_char ('P', "Pro","Proline");
  add_normal_char ('Q', "Gln","Glutamine");
  add_normal_char ('R', "Arg","Arginine");
  add_normal_char ('S', "Ser","Serine");
  add_normal_char ('T', "Thr","Threonine");
  add_normal_char ('V', "Val","Valine");
  add_normal_char ('W', "Trp","Tryptophan");
  add_normal_char ('Y', "Tyr","Tyrosine");
}
ExtAminoAlphabet::ExtAminoAlphabet (const char *nm) :
AminoAlphabet(NULL)
{
  reset_name (nm);
  add_wild_card ('B', "ND", "Asn/Asp", "Asn/Asp");
  add_wild_card ('Z', "QE", "Gln/Glu", "Gln/Glu");
  add_all_match ('X');
}
ExtDNAAlphabet::ExtDNAAlphabet(const char *nm):
DNAAlphabet(NULL)
{
    reset_name(nm);
    add_wild_card('K', "GT");
    add_wild_card('W', "AT");
    add_wild_card('Y', "CT", "Pyrimidine", "Pyrimidine");
    add_wild_card('M', "AC");
    add_wild_card('R', "AG", "Purine", "Purine");
    add_wild_card('S', "GC");
    add_wild_card('V', "AGCRMS");
    add_wild_card('B', "GCTSKY");
    add_wild_card('D', "AGTRWK");
    add_wild_card('H', "ACTMWY");
    add_all_match('N');
    add_alias('X','N');
}
Base
ExtDNAAlphabet::complement(const Base b) const
{
  if (b.is_null())
    return b;
  if (b.is_wild()) {
    switch (b.raw_int()) {
    case K: return Base::base_int(M); // GT <-> CA
    case M: return Base::base_int(K);
    case W: return Base::base_int(W); // AT <-> TA
    case Y: return Base::base_int(R); // AG <-> TC
    case R: return Base::base_int(Y);
    case S: return Base::base_int(S); // GC <-> CG
    case V: return Base::base_int(B); // AGC <-> TCG
    case B: return Base::base_int(V); // AGC <-> TCG
    case D: return Base::base_int(H); // AGT <-> TCA
    case H: return Base::base_int(D); // AGT <-> TCA
    case N: return Base::base_int(N);
    case X: return Base::base_int(X);
    default:
      cerr << "ExtDNAAlphabet::complement:  unknown wildcard\n";
      assert (0);
    }
  }
  return Base::base_int (b.raw_int()^3);
}
int
ExtDNAAlphabet::is_pyrimidine (const Base b) const
{return (b.canon()).raw_int()==A || b.canon().raw_int()==G || b.canon().raw_int() == R;}
int
ExtDNAAlphabet::is_purine     (const Base b) const
{return b.canon().raw_int()==C || b.canon().raw_int()==T
 || b.canon().raw_int() == Y;}
int
ExtDNAAlphabet::same_group (const Base b1, const Base b2) const
{return ( (is_pyrimidine(b1) && is_pyrimidine (b2)) ||
	  is_purine(b1) && is_purine (b2));}
int
ExtDNAAlphabet::is_complement (const Base b1, const Base b2) const
{return b1.canon().raw_int() == complement (b2.canon()).raw_int();}


/*inline*/ char
Alphabet::to_char (const Base b1, const VarEnum v) const
{
  assert (v == NO_VARS);   // to_char cannot deal with variants.

  if (b1.is_null())
    return Base::null_char();
  else
    return _chars[b1.raw_int()];	// String[] is not const.
}

// CHANGE LOG:
// 30 March 2004 Kevin Karplus
//	Fixed some variable names that shadowed member names.
// Fri Mar 25 17:58:17 PST 2005 Kevin Karplus
//	Added static_cast for array subscripts that were char
// Sat Jun 18 10:38:51 PDT 2005 Kevin Karplus
//	Modified load_alphabet_file to use NamedObject:read_many_new
