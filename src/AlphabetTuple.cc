//AlphabetTuple.cc
// Kevin Karplus
// created 26 April 1995

// ABSTRACT
//	Contains the implementation for AlphabetTuple and BaseTuple classes

#include "EqualStrings.h"
#include "AlphabetTuple.h"
#include "Input.h"
#include <string.h>
#include <stdio.h>

// construct from singleton Alphabet
AlphabetTuple::AlphabetTuple(const Alphabet *a0)
{
    assert (a0 != NULL);
    NumAlphabets=1;
    alphabets = new const Alphabet*[NumAlphabets];
    alphabets[0] = a0;
    NumNormal = a0->norm_length();
    NumNormWc = a0->norm_wc_length();
    set_name(a0->name());
}

// construct from pair of Alphabets
AlphabetTuple::AlphabetTuple(const Alphabet *a0, const Alphabet *a1)
{
    assert (a0 != NULL && a1 != NULL);
    NumAlphabets=2;
    alphabets = new const Alphabet*[NumAlphabets];
    alphabets[0] = a0;
    alphabets[1] = a1;
    NumNormal = a0->norm_length() * a1->norm_length();
    NumNormWc = a0->norm_wc_length() * a1->norm_wc_length();
    
    // concatenate the alphabet names
    char *nm = new char[strlen(a0->name()) + strlen(a1->name()) + 2];
    int i=0;
    for(const char *nm0=a0->name(); *nm0; nm0++)
    {    nm[i++] = *nm0;
    }
    nm[i++] = ':';
    for(const char *nm1=a1->name(); *nm1; nm1++)
    {    nm[i++] = *nm1;
    }
    nm[i++] = 0;
    set_name(nm);
    delete [] nm;
}

// construct from triple of Alphabets
AlphabetTuple::AlphabetTuple(const Alphabet *a0, const Alphabet *a1, 
		const Alphabet *a2)
{
    assert (a0 != NULL  && a1 != NULL  && a2 != NULL);
    NumAlphabets=3;
    alphabets = new const Alphabet*[NumAlphabets];
    alphabets[0] = a0;
    alphabets[1] = a1;
    alphabets[2] = a2;
    NumNormal = a0->norm_length() * a1->norm_length() * a2->norm_length();
    NumNormWc = a0->norm_wc_length() * a1->norm_wc_length() * a2->norm_wc_length();

    // concatenate the alphabet names
    char *nm = new char[strlen(a0->name()) + strlen(a1->name()) + strlen(a2->name()) + 3];
    int i=0;
    for(const char *nm0=a0->name(); *nm0; nm0++)
    {    nm[i++] = *nm0;
    }
    nm[i++] = ':';
    for(const char *nm1=a1->name(); *nm1; nm1++)
    {    nm[i++] = *nm1;
    }
    nm[i++] = ':';
    for(const char *nm2=a2->name(); *nm2; nm2++)
    {    nm[i++] = *nm2;
    }
    nm[i++] = 0;
    set_name(nm);
    delete [] nm;
}

// construct from array of Alphabets
AlphabetTuple::AlphabetTuple(int n, const Alphabet **a)
{
    assert (n>=0);
    NumAlphabets=n;
    NumNormal=1;
    NumNormWc = 1;
    alphabets = new const Alphabet*[NumAlphabets];
    int sum_name_lengths=0;
    for (int i=n-1; i>=0; i--)
    {	assert(a[i] != NULL);
    	alphabets[i] = a[i];
	NumNormal *= a[i]->norm_length();
	NumNormWc *= a[i]->norm_wc_length();
	sum_name_lengths+=strlen(a[i]->name());
    }
    
    // concatenate the alphabet names
    char *nm = new char[sum_name_lengths+n];
    int i=0;
    for (int al=0; al<n; al++)
    {   for(const char *nm_index=a[al]->name(); *nm_index; nm_index++)
	{    nm[i++] = *nm_index;
	}
	nm[i++] = ':';
    }
    nm[i-1] = 0;
    set_name(nm);
    delete [] nm;

}


// copy constructor
AlphabetTuple::AlphabetTuple(const AlphabetTuple &a)
	: NamedObject(a.name(), a.help())
{
    if (a.isNumeric())
    {
        NumAlphabets = 1;
        alphabets = 0;
        NumericAlph = a.NumericAlph;
        NumNormal = a.NumNormal;
        NumNormWc = a.NumNormWc;
        return;
    }
        
    int n=a.num_alphabets();
    NumAlphabets=n;
    NumNormal=1;
    NumNormWc = 1;
    alphabets = new const Alphabet*[NumAlphabets];
    for (int i=n-1; i>=0; i--)
    {	alphabets[i] = a[i];
	NumNormal *= a[i]->norm_length();
	NumNormWc *= a[i]->norm_wc_length();
    }
    
}

AlphabetTuple::AlphabetTuple(const NumericAlphabet &na)
{
    alphabets = 0;
    NumericAlph = na;
    NumAlphabets = 1;
    NumNormal = NumericAlph.length();
    NumNormWc = NumNormal;
    char nm[100];
    sprintf(nm,"%d",NumNormal);
    set_name(nm);
}

void AlphabetTuple::set_length(int n)
{
    assert(n >= 0  &&  isNumeric());
    NumericAlph.set_length(n);
    NumNormal = n;
    NumNormWc = NumNormal;
    char nm[100];
    sprintf(nm,"%d",NumNormal);
    set_name(nm);
}

// indexes are computed like array subscripts, with the last part of the
//	tuple being the least significant part
// Note: all Bases in the BaseTuple must be normal (not wildcard)
int AlphabetTuple::index(const BaseTuple &bt) const
{
    assert(!isNumeric());

    int ofs = 0; //1-dimensional offset accumulator of base tuple
    int i;

    for (i = 0; i < NumAlphabets; i++)
    {
        const Alphabet *a = alphabets[i];

        assert(bt[i].is_normal());
        ofs *= a->norm_length();
        ofs += a->index(bt[i]);
    }

    return ofs;
}

BaseTuple* AlphabetTuple::unindex(int index_val) const
{
    BaseTuple *bt=new BaseTuple(*this);
    unindex(index_val, *bt);
    return bt;
}


void AlphabetTuple::unindex(int ind, BaseTuple &bt) const
{
    assert(!isNumeric());
    assert( 0 <= ind && ind <= NumNormal );

    int i;

    for (i = NumAlphabets - 1; i >= 0; i--)
    {
        const Alphabet *a = alphabets[i];
        int len = a->norm_length();
        bt[i] = a->unindex(ind % len);
        ind /= len;
    }
    assert(ind == 0);
}



void AlphabetTuple::print_unindex(ostream& out, int index_val) const
{
    assert(0<=index_val && index_val<NumNormal);
    if (isNumeric())
    {   out << index_val;
    }
    else
    {   int shift = NumNormal;
        for (int i=0; i<NumAlphabets; i++)
        {   const Alphabet* a= alphabets[i];
            int mask = a->norm_length();
            shift /= mask;
            out << a->to_char(a->unindex((index_val/shift) % mask));
        }
        assert(shift==1);
    }
}


int AlphabetTuple::norm_wc_index(const BaseTuple &bt) const
{
    assert(!isNumeric());

    int ofs = 0; //1-dimensional offset accumulator of base tuple
    int i;

    for (i = 0; i < NumAlphabets; i++)
    {
        const Alphabet *a = alphabets[i];

        ofs *= a->norm_wc_length();
        ofs += a->index(bt[i]);
    }

    return ofs;
}

void AlphabetTuple::norm_wc_unindex(int ind, BaseTuple &bt) const
{
    assert(!isNumeric());

    int i;

    assert( 0 <= ind && ind <= NumNormWc);
    for (i = NumAlphabets - 1; i >= 0; i--)
    {
        const Alphabet *a = alphabets[i];
        int len = a->norm_wc_length();
        bt[i] = a->unindex(ind % len);
        ind /= len;
    }
}



AlphabetTuple* read_AlphabetTuple(istream &in, int dimension)
{
    if (dimension<=0) in >> dimension;
    assert(dimension >0);
    
    const Alphabet** alphabets = new const Alphabet* [dimension];
    char word[100];
    for (int i=0; i< dimension; i++)
    {   get_word(in, word);
	alphabets[i] = Alphabet::name_to_alphabet(word);
	if (!alphabets[i])
	{   cerr << "Alphabet " << word << " not known.\n";
	    return 0;
	}
    }
    AlphabetTuple *ret_tuple = new AlphabetTuple(dimension,alphabets);
    delete [] alphabets;
    return ret_tuple;
}

AlphabetTuple *read_NumericAlphabetTuple(istream &in)
{
    int len;

    in >> len;

    assert(len >= 0);
    NumericAlphabet na(len);

    return new AlphabetTuple(na);
}

// Reads a command of the form
//	Alphabet= <alphabet_name>
//	AlphabetPair= <alphabet_name> <alphabet_name>
//	AlphabetTriple= <alphabet_name> <alphabet_name> <alphabet_name>
//	AlphabetTuple=  <number> <alphabet_name> ... <alphabet_name>
// If the firstword is not recognized, it is looked up as an alphabet name,
//	as if preceded by "Alphabet="
AlphabetTuple* read_AlphabetTuple_command(istream &in)
{
    char word[100];
    get_word(in,word);
    int len = strlen(word);
    if (len && word[len-1] == '=') 
    {    word[--len] = '\0';
    }
    else // swallow the =-sign
    {    char c;
   	in >> c;
	if (c!='=') in.putback(c);
    }
    if (EqualStrings(word, "Alphabet",1))
	return read_AlphabetTuple(in,1);
    else if (EqualStrings(word, "AlphabetPair",1))
	return read_AlphabetTuple(in,2);
    else if (EqualStrings(word, "AlphabetTriple",1))
	return read_AlphabetTuple(in,3);
    else if (EqualStrings(word, "AlphabetTuple",1))
	return read_AlphabetTuple(in);
    else if (EqualStrings(word, "NumericAlphabet"))
        return read_NumericAlphabetTuple(in);
    else
    {   const Alphabet* alphabets;
	alphabets = Alphabet::name_to_alphabet(word, ZeroIfNew);
	if (!alphabets)
	{   cerr << "Alphabet " << word << " not known.\n";
	    return 0;
	}
	return new AlphabetTuple(alphabets);
    }    	
}


ostream& operator << (ostream& out, const AlphabetTuple& at)
{
    if (at.isNumeric())
    {
        out << at.num_normal();
        return out;
    }

    out << at.num_alphabets();
    for (int i=0; i<at.num_alphabets(); i++)
    	out << " "  << at[i]->name();
    out << " ";
    return out;
}

ostream& operator << (ostream& out, const BaseTuple& bt)
{
    assert(!bt.al->isNumeric());

    for (int i=0; i<bt.al->num_alphabets(); i++)
    {   const Alphabet* a=   (*(bt.al))[i];
    	out << a->to_char(bt[i]);
    }
    return out;
}

istream& operator >> (istream& in, BaseTuple& bt)
{
    assert(!bt.al->isNumeric());

    for (int i=0; i<bt.al->num_alphabets(); i++)
    {   char c;
   	in >> c;
	const Alphabet* a=   (*(bt.al))[i];
    	bt[i] = a->to_base(c);
    }
    return in;
}

void AlphabetTuple::print_command(ostream & out) const
{
    if (isNumeric())
    {
        out << "NumericAlphabet = " << NumericAlph.length() << "\n";
        return;
    }

    switch(num_alphabets())
    {   case 1:  out << "Alphabet = " << alphabets[0]->name() << "\n";
		break;
	case 2:  out << "AlphabetPair = " << alphabets[0]->name() 
		<< " " << alphabets[1]->name() 
		<< "\n";
		break;
	case 3:  out << "AlphabetTriple = " << alphabets[0]->name() 
		<< " " << alphabets[1]->name() 
		<< " " << alphabets[2]->name() 
		<< "\n";
		break;
	default:  out << "AlphabetTuple = " << (*this) << "\n";
		break;
    }

}

//CHANGE LOG
// 24 Nov 1995 Kevin Karplus
//	added print_command and copy constructor
// 	eliminated strcasecmp, using EqualStrings instead
// 26 Feb 1996 Spencer Tu
//      Added NumericAlphabet, and wild card indexing functions.
//      NumericAlphabet here is just a hack, should be replaced
//      later.  Fixed problem with index().
// 19 Oct 1996 Spencer Tu
//      Replaced print_unindex() with the fixed version in ultimate src.
// 30 March 2004 Kevin Karplus
//	Fixed some variable names that shadowed member names.
// Sat Jun 18 18:29:11 PDT 2005 Kevin Karplus
//	Added names to constructors for AlphabetTuple.
//	Removed some commented-out old code.
