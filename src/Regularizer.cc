// Regularizer.cc

#include <stdlib.h>     // for atoi
#include <memory.h>	// for bzero
#include <fstream>

#include "Regularizer.h"
#include "log2.h"
#include "minmax.h"
#include "abs.h"	// abs

const float REALBIGCOST = 1.0e33;

// To do:
//	Modify verify_partials1 and verify_partials2 to
//	take an ostream for output of error messages


IdObject Regularizer::ID("Regularizer",0,0,
"a pure virtual type, used as the base for all regularizer classes.\nA regularizer is a method for converting an observed sample to an\nestimate of a distribution.");

NameToPtr* Regularizer::CommandTable = 0;

int Regularizer::zero_second_deriv(void) const
{    return 1;
}

// default implementation of 2nd partials, assumes they are 0.
// Should be overloaded with correct 2nd partial computation
void Regularizer::partials2(float *part1, float *part2, int i,
    	const float *TrainCounts)
{
    bzero(part2, num_parameters()*sizeof(float));
    partials1(part1, i, TrainCounts);
}

int Regularizer::verify_partials1(const float*TrainCounts, float tolerance)
{
    int ret=1;	// return value, set to 0 if bad partial found
    float *x = new float[alphabet_size()];
    get_modified_counts(TrainCounts, x);
    int part1_col_len = num_parameters();
    float *part1 = new float[alphabet_size()*part1_col_len];
    int i;
    for (i=alphabet_size()-1; i>=0; i--)
	partials1(&part1[index_2D(part1_col_len, i, 0)], i, TrainCounts);

    float *x_plus = new float[alphabet_size()];
    float *x_minus = new float[alphabet_size()];
    for (int p=num_parameters()-1; p>=0; p--)
    {   float save=parameter(p);
	double delta = kk::abs(save*0.01);
	set_parameter(p, save+delta);
	get_modified_counts(TrainCounts, x_plus);
	set_parameter(p, save-delta);
	get_modified_counts(TrainCounts, x_minus);
	for (i=alphabet_size()-1; i>=0; i--)
	{   double diff= 0.5*(x_plus[i] - x_minus[i]);
	    double loosened_tolerance = 
		max(static_cast<double>(tolerance), 1.e-6 *x[i] / kk::abs(diff));
	    double deriv = diff/delta;
	    double abserr = kk::abs(deriv - part1[index_2D(part1_col_len, i, p)]);
	    if (abserr > loosened_tolerance 
		&& kk::abs(abserr/ deriv) > loosened_tolerance)
	    {    cerr << "First partial outside tolerance for "
			<< type()->name() << " " << name()
			<< "\n parameter " << p << " = " << save
			<< " letter " << i << "(";
		alphabet_tuple()->print_unindex(cerr, i);
		cerr << ")\n";
		cerr << "counts=";
		for (int j=0; j<alphabet_size(); j++)
		    cerr << " " << TrainCounts[j];
		cerr << "\n";
		cerr << "computed partial = "
                        << part1[index_2D(part1_col_len, i, p)]
			<< " but X_plus[" << i << "]= " << x_plus[i]
			<< " X[" << i << "]= " << x[i]
			<< " deriv_from_diff = " << deriv
			<< "\n";
	        ret = 0;
	    }
	}
	set_parameter(p, save);
    }

    delete [] x_plus;
    delete [] x_minus;
    delete [] part1;
    delete [] x;
    return ret;
}

int Regularizer::verify_partials2(const float*TrainCounts, float tolerance)
{
    int ret=1;	// return value, set to 0 if bad partial found
    float *x = new float[alphabet_size()];
    get_modified_counts(TrainCounts, x);
    const int part_col_len = num_parameters();
    float *part1 = new float[alphabet_size()*part_col_len];
    float *part2 = new float[alphabet_size()*part_col_len];
    
    float *part1_from_partials1 = new float[num_parameters()];
    int i;
    for (i=alphabet_size()-1; i>=0; i--)
    {   partials2(&part1[index_2D(part_col_len, i, 0)],
                  &part2[index_2D(part_col_len, i, 0)], i, TrainCounts);
	partials1(part1_from_partials1, i, TrainCounts);
	for (int p=num_parameters()-1; p>=0; p--)
	{   if ( fabs(part1[index_2D(part_col_len, i, p)] - part1_from_partials1[p])
		> 1.e-8)
	    {	cerr << "partials1 and partials2 disagree on first partial "
			<< type()->name() << " " << name()
			<< "\n parameter " << p << " = " << parameter(p)
			<< " letter " << i << "(";
		alphabet_tuple()->print_unindex(cerr, i);
		cerr << ")\n";
		cerr << "counts=";
		for (int j=0; j<alphabet_size(); j++)
		    cerr << " " << TrainCounts[j];
		cerr << "\n";
		cerr << " partials1 says " << part1_from_partials1[p]
			<< " but partials2 says "
                        << part1[index_2D(part_col_len, i, p)] << "\n";
		ret=0;
	    }
	}
    }

    float *x_plus = new float[alphabet_size()];
    float *x_minus = new float[alphabet_size()];
    for (int p=num_parameters()-1; p>=0; p--)
    {   float save=parameter(p);
	double delta = kk::abs(save*0.01);
	set_parameter(p, save+delta);
	get_modified_counts(TrainCounts, x_plus);
	set_parameter(p, save-delta);
	get_modified_counts(TrainCounts, x_minus);
	
	for (i=alphabet_size()-1; i>=0; i--)
	{   double diff= 0.5*(x_plus[i] - x_minus[i]);
	    double loosened_tolerance = 
		max(static_cast<double>(tolerance), 1.e-6 *x[i] / kk::abs(diff));
	    double deriv = diff/delta;
	    double abserr = kk::abs(deriv - part1[index_2D(part_col_len, i, p)]);
	    if (abserr > loosened_tolerance 
		&& kk::abs(abserr/ deriv) > loosened_tolerance)
	    {    cerr << "First partial outside tolerance for "
			<< type()->name() << " " << name()
			<< "\n parameter " << p << " = " << save
			<< " letter " << i << "(";
		alphabet_tuple()->print_unindex(cerr, i);
		cerr << ")\n";
		cerr << "counts=";
		for (int j=0; j<alphabet_size(); j++)
		    cerr << " " << TrainCounts[j];
		cerr << "\n";
		cerr << "computed partial = "
                        << part1[index_2D(part_col_len, i ,p)]
			<< " but deriv_from_diff = " << deriv << "\n"
			<< " X_plus[" << i << "]= " << x_plus[i]
			<< " X[" << i << "]= " << x[i]
			<< " X_minus[" << i << "]= " << x_minus[i]
			<< "\n";
	        ret = 0;
	    }
	    
	    double diff2 = x_plus[i] + x_minus[i] - 2*x[i];
	    loosened_tolerance = max(static_cast<double>(tolerance),  1.e-6*x[i]/kk::abs(diff2));
	    double deriv2 = diff2/(delta*delta);
	    double abserr2 = kk::abs(deriv2 - part2[index_2D(part_col_len, i, p)]);
	    if (abserr2 > loosened_tolerance
		&& kk::abs(abserr2/ deriv2) > loosened_tolerance)
	    {    cerr << "Second partial outside tolerance for "
			<< type()->name() << " " << name()
			<< "\n parameter " << p << " = " << save
			<< " letter " << i << "(";
		alphabet_tuple()->print_unindex(cerr, i);
		cerr << ")\n";
		cerr << "counts=";
		for (int j=0; j<alphabet_size(); j++)
		    cerr << " " << TrainCounts[j];
		cerr << "\n";
		cerr << "computed partial2 = "
                        << part2[index_2D(part_col_len, i, p)]
			<< " but deriv2_from_diff = " << deriv2 << "\n"
			<< " X_plus[" << i << "]= " << x_plus[i]
			<< " X[" << i << "]= " << x[i]
			<< " X_minus[" << i << "]= " << x_minus[i]
			<< "\n";
	        ret = 0;
	    }
	    
	}
	set_parameter(p, save);
    }

    delete [] x_minus;
    delete [] x_plus;
    delete [] part1_from_partials1;
    delete [] part2;
    delete [] part1;
    delete [] x;
    return ret;
}

void Regularizer::get_probs(const float* TrainCounts,  float* EstProbs)
{
    get_modified_counts(TrainCounts, EstProbs);
    double sum=0.0;
    int i;
    for(i = alphabet_size()-1; i >= 0; i--) 
        sum += EstProbs[i];

    if (sum<=0)
    {    // predict flat for all zeros
        float flat = 1./alphabet_size();
	for(i = alphabet_size()-1; i >= 0; i--) 
	    EstProbs[i] = flat;
        return;
    }
    
    double scale = 1./sum;
    for(i = alphabet_size()-1; i >= 0; i--) 
    	EstProbs[i] *= scale;
    
}

// return cost of column in bits,
//	after estimating the posterior counts (probs).
float Regularizer::encodingCostForColumnCounts(
  const float *realProbs,
  const float *counts,
  float *EstProbs) 
{
    get_probs(counts, EstProbs);
    register float ent=0.;
    for(int i = alphabet_size()-1; i >= 0; i--) 
    {	if (realProbs[i]<=0.) continue;
        if (EstProbs[i] <=0.) return(REALBIGCOST);
      	ent -= realProbs[i] * log2(EstProbs[i]);
    }
    return(ent);
}

void Regularizer::print_order(ostream & out) const
{
    out << "Order = ";
    for (int i=0; i<alphabet_size(); i++)
    {    out << " ";
        alphabet_tuple()->print_unindex(out,i);
    }
    out << "\n";

}



void Regularizer::read_order(istream &in)
{
    if (InputOrder) delete [] InputOrder;
    // read the alphabet order
    InputOrder = new int[alphabet_size()];
    BaseTuple bt(*alphabet_tuple());
    for (int i=0; i<alphabet_size(); i++)
    {   in >> bt;
	InputOrder[i] = alphabet_tuple()->index(bt);
	assert( 0<= InputOrder[i] && InputOrder[i] < alphabet_size());
    }
}    

// the following static functions are commands that are understood by all
// Regularizer classes

static int ReadAlphabets(int num, istream &in, 	Regularizer *change)
{
    const Alphabet** alphs = new const Alphabet* [num];
    
    char word[100];
    for (int i=1; i<=num; i++)
    {   get_word(in, word);
	const Alphabet* al= Alphabet::name_to_alphabet(word, ZeroIfNew);
	if (!al)
        {
            cerr << "Unrecognized " << i << "th alphabet " << word 
                 << " for regularizer " 
                 << (change->name()? change->name(): " that has no name")
                 << "\n" << flush;
            delete [] alphs;
	    return 0;
        }
        alphs[i-1] = al;
    }
    change->set_alphabet_tuple(new AlphabetTuple(num, alphs));
    delete [] alphs;
    return 1;
}


static int ReadAlphabet(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   return ReadAlphabets(1, in, change);
}

static int ReadAlphabetPair(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   return ReadAlphabets(2, in, change);
}

static int ReadAlphabetTriple(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   return ReadAlphabets(3, in, change);
}

static int ReadAlphabetTuple(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   int num;
    in >> num;
    if (num<=0)
    {    cerr << "Error: " << self->name()
		<< " needs to start with the number of alphabets, not "
		<< num << "\n" << flush;
	return 0;
    }
    return ReadAlphabets(num, in, change);
}

static int ReadNumericAlphabet(istream &in,
                Regularizer *change,
                RegInputCommand *self)
{
    char word[128];
    get_word(in, word);
    
    int len = atoi(word);
    NumericAlphabet nalph(len);

    change->set_alphabet_tuple(new AlphabetTuple(nalph));
    return 1;
}

static int ReadName(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   char word[500];
    get_word(in, word);
    change->set_name(word);
    return 1;
}

static int ReadOrder(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   
    change->read_order(in);
    return 1;
}

// not static, so that derived Regularizers can treat  certain
//	keywords as comments
int ReadComment(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   
    SkipSeparators(in, 1, '\n');
    return 1;
}

static int VerifyClassName(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{   char word[100];
    get_word(in, word);
    const IdObject *end_id = IdObject::id(word);
    if (end_id != change->type())
    {    cerr << "Warning: " << self->name() << word << " doesn't match "
		<< change->type()->name() << "\n" << flush;
    }
    // continue if "ClassName", stop if "EndClassName"
    return EqualStrings(self->name(), "ClassName", 1);
}

void Regularizer::init_command_table()
{   assert(!CommandTable);
    CommandTable = new NameToPtr(13);
    CommandTable->ignore_case();
    CommandTable->AddName(new RegInputCommand("Alphabet", ReadAlphabet));
    CommandTable->AddName(new RegInputCommand("AlphabetPair",
				ReadAlphabetPair)); 
    CommandTable->AddName(new RegInputCommand("AlphabetTriple",
				ReadAlphabetTriple)); 
    CommandTable->AddName(new RegInputCommand("AlphabetTuple",
				ReadAlphabetTuple)); 
    CommandTable->AddName(new RegInputCommand("NumericAlphabet", ReadNumericAlphabet));
    CommandTable->AddName(new RegInputCommand("Name", ReadName));
    CommandTable->AddName(new RegInputCommand("Order", ReadOrder));
    CommandTable->AddName(new RegInputCommand("Comment", ReadComment));
    CommandTable->AddName(new RegInputCommand("ClassName", VerifyClassName));
    CommandTable->AddName(new RegInputCommand("EndClassName",
				VerifyClassName)); 
}

int Regularizer::read_knowing_type(istream &in)
{   
// THE FOLLOWING TEST which forces alphabet to be set before reading
// in the regularizer seems to be buggy.  Let's try eliminating it.
//    if (! alphabet_tuple()) set_alphabet();	// set default alphabet
//    if (! alphabet_tuple())
//    {    cerr << "Error: no alphabet specified, and no default set\n"
//		<< flush;
//	return 0;
//    }
    
    if (! command_table()) init_command_table();
    if (! Regularizer::CommandTable) Regularizer::init_command_table();

    char word[100];
    while (in.good())
    {   get_word(in, word, '=');
	RegInputCommand *comm = dynamic_cast<RegInputCommand *>(
		command_table()->FindOldName(word, ZeroIfNew));
        if (!comm) comm = dynamic_cast<RegInputCommand *>(
		Regularizer::CommandTable->FindOldName(word, ZeroIfNew));
	if (comm)
	{   if (!comm->execute(in, this)) return 1;
	}
	else
	{    cerr << "Unrecognized keyword: " << word 
		<< " for type " << type()->name()
		<< " " << name()
		<< "\n" << flush;
	}
    }
    return 0;
}

Regularizer* Regularizer::read_new(istream& in,
	    IdObject* required_type)
{   NamedClass *n=NamedClass::read_new(in);
    if (n && !n->is_a(required_type))
    {    cerr << "Error: was expecting a " << required_type->name()
	    << " not a " << n->type()->name() << "\n" << flush;
	 delete n;
	 return 0;
    }
    return dynamic_cast<Regularizer*>(n);
}


Regularizer* Regularizer::read_new(const char* filename,
	    IdObject* required_type)
{   ifstream in(filename);
    if (!in.good()) 
    {    cerr << "Error: couldn't open " << filename
	    << " to read Regularizer.\n" << flush; 
	return 0;
    }
    Regularizer *r=read_new(in, required_type);
    if (r && !r->name())
	r->set_name(filename);
    return r;
}


void Regularizer::write_knowing_type(ostream &out) const
{   if (name()) out << "Name = " << name() << "\n";
    if (alphabet_tuple()) 
    { alphabet_tuple()->print_command(out);
        if (!alphabet_tuple()->isNumeric())
            print_order(out);
    }
    
    // Note: print_order must be last, for this to work with SubstPseudoReg
}


//ChangLog
//
// 19 Oct 1996 Spencer Tu
//      Changed variable length arrays in verify_partials1
//      and verify_partials1 and ReadAlphabets to
//      dynamically allocated arrays.
//      Added member index_2D to support the above change.
// 30 Oct 2001 Kevin Karplus
//	Removed forcing AlphabetTuple to be set before reading.
// 17 July 2003 George Shackelford
//	Changed abs to kk::abs to use namespace of Abs.h
// 15 March 2004 Kevin Karplus
//	Changed old-style cast to static_cast
// 30 March 2004 Kevin Karplus
//	Fixed != test to abs(x-y)>small number.
