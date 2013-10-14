#include <math.h>

#include "DirichletReg.h"

#include "OrderSort.h"// for comment on symbol frequency in printout
#include "Input.h"
#include "LogGamma.h"
#include "EqualStrings.h"

#include "log2.h"
#include "minmax.h"
#include "abs.h" // for abs

// To Do:
// 2) Turn off extra caching, now that the LogGamma function is
//	working better?  Maybe not--the extra caching only involves
//	log, not lgamma calls, and is pretty cheap.


static NamedClass* create_DirichletReg(void) {return new DirichletReg;}
IdObject DirichletReg::ID("DirichletReg", create_DirichletReg,
	Regularizer::init_is_a,
"a regularizer based on a Dirichlet mixture. \nEach component of the mixture is a set of pseudocounts defining a \nDirichlet prior on the possible distributions of letters. \nThe mixture coefficients specify the probability associated with each \ndifferent set of psuedocounts. \nNote that a single-component Dirichlet mixture is equivalent to a \nsimple pseudocount regularizer."
);
NameToPtr* DirichletReg::CommandTable=0;

// Pre-compute gamma(z_c(i) + s) for s=0, 1, 2, ... , MaxPlus
#ifdef  DirichletReg_No_Extra_Caching
const int DirichletReg::MaxPlus = 0;
#else
const int DirichletReg::MaxPlus = 30;
#endif


void DirichletReg::clear(void)
{
    FreezeComp = 0;
    FreezeMix = 0;
    
    NumComp=AllocComp=0;
    MixtureCoeff=0;
    Components=0;
}

DirichletReg::DirichletReg(const MLPReg &dir)
	: Regularizer(dir.alphabet_tuple(),dir.name())
{   clear();
    alloc(1);
    AddComponent(1, dir.pseudocounts());
}

DirichletReg::DirichletReg(const MLZReg &dir)
	: Regularizer(dir.alphabet_tuple(),dir.name())
{   clear();
    alloc(1);
    float *counts = new float[alphabet_size()];
    for (int i=alphabet_size()-1; i>=0; i--)
	counts[i] = dir.zero_offset();
    AddComponent(1, counts);
    delete [] counts;
}
// copy constructor

DirichletReg::DirichletReg(const DirichletReg &dir)
	: Regularizer(dir.alphabet_tuple(),dir.name())
{
    FreezeComp = dir.FreezeComp;
    FreezeMix = dir.FreezeMix;
        
    NumComp=AllocComp=0;
    MixtureCoeff=0;
    Components=0;
    alloc(dir.num_components());

    for (int c=0; c<dir.num_components(); c++)
    {    AddComponent(dir.mixture_coeff(c), dir.component(c));
    }
}

DirichletReg::DirichletReg(void) : Regularizer() 
{   clear();
    alloc(0);
}
	
	
DirichletReg::DirichletReg(const Alphabet *a, const char *nm, int size)
    	: Regularizer(a,nm)
{
    clear();
    alloc(size);
}

void DirichletReg::alloc(int num)
   {
    NumComp=0;
    if (num==0)
    {	MixtureCoeff = 0;
	Components = 0;
	LogGammaComp = 0;
        Lgamma1Comp = 0;
        Lgamma2Comp = 0;

	SumComponents = 0;
        Lgamma1Sum = 0;
        Lgamma2Sum = 0;
	LogBetaComp = 0;
	
	CompProbs = 0;
	X=0;
	Lgamma1_SumCounts=0;
	Lgamma1_Counts=0;
	Lgamma1CountsInvalid=1;
	Lgamma2_SumCounts=0;
	Lgamma2_Counts=0;
	Lgamma2CountsInvalid=1;
	
	return;
    }
    assert(alphabet_tuple() != 0);
    AllocComp = num;
    int AllocAlpha= num*alphabet_size();
    MixtureCoeff = new float[AllocComp];
    Components = new float[AllocAlpha];
    LogGammaComp = new double[AllocAlpha*(MaxPlus+1)];
    Lgamma1Comp = new double[AllocAlpha];
    Lgamma2Comp = new double[AllocAlpha];
    
    SumComponents = new double[AllocComp];
    Lgamma1Sum = new double[AllocComp];
    Lgamma2Sum = new double[AllocComp];
    LogBetaComp = new double[AllocComp];
    
    // the following are caches for a specific sample.
    CompProbs = new double[AllocComp];
    X = new double[alphabet_size()];
    Lgamma1_SumCounts = new double[AllocComp];
    Lgamma1_Counts= new double[AllocAlpha];
    Lgamma1CountsInvalid=1;
    Lgamma2_SumCounts = new double[AllocComp];
    Lgamma2_Counts= new double[AllocAlpha];
    Lgamma2CountsInvalid=1;
}


void DirichletReg:: delete_all(void)
{   NumComp=AllocComp=0;
    delete [] MixtureCoeff;
    delete [] Components;
    delete [] LogGammaComp;
    delete [] Lgamma1Comp;
    delete [] Lgamma2Comp;
    delete [] SumComponents;
    delete [] Lgamma1Sum;
    delete [] Lgamma2Sum;
    delete [] LogBetaComp;
    delete [] X;
    delete [] CompProbs;
    delete [] Lgamma1_SumCounts;
    delete [] Lgamma1_Counts;
    delete [] Lgamma2_SumCounts;
    delete [] Lgamma2_Counts;
}

DirichletReg* DirichletReg::posterior_mixture(const float* TrainCounts)
{
    DirichletReg* dir = new DirichletReg;
    dir->set_alphabet_tuple(new AlphabetTuple(* alphabet_tuple()));
    dir->set_name(name());
    dir->clear();
    dir->alloc(num_components());
    
    component_probs(TrainCounts, SumCounts, CompProbs);
    
    float *new_comp = new float[alphabet_size()];
    for (int c=0; c<num_components(); c++)
    {   float* comp = component(c);
	for (int j= alphabet_size()-1; j>=0; j--)
	    new_comp[j] = comp[j] + TrainCounts[j];
	dir->AddComponent(CompProbs[c], new_comp);
    }
    delete [] new_comp;
    return dir;
}


// return 1 if x is a small integer, and set ix to that integer
// (ix undefined if x not a small integer)
inline int DirichletReg::is_small_int(float x, int &ix)
{   
    if (x==0)
    {   ix=0;		// short-cut for most common case
	return 1;
    }
#ifndef  DirichletReg_No_Extra_Caching
    const float MAX_ERROR=0.005;	// allow for slight error
    float fl = floor(x);
    if (x-fl < MAX_ERROR)
    {    ix = (int) fl;
	return (ix>= 0 && ix <=MaxPlus);
    }
    else if (x-fl > 1.-MAX_ERROR)
    {    ix = (int) fl +1;
	return (ix>= 0 && ix <=MaxPlus);
    }
#endif
    return 0;
}

// Note: if c==NumComp, then adding new component, not updating old one
void DirichletReg::set_component(int c, int lett, float z)
{   assert(0<=lett && lett<alphabet_size());
    if (z<=0.) z=1.e-8;		// replace zeros with small number
    assert(0<z && z <1.e25);
    
    // remove old gamma component from BetaComp
    if (c<NumComp)
	LogBetaComp[c] -= lgamma_component(c,lett);
    
    // fill in the new gamma_component cache table
    double lgamma_tmp;
    LogGamma_derivs(z, lgamma_tmp,
		lgamma1_comp(c)[lett],lgamma2_comp(c)[lett]); 
    for (int plus=0; plus<=MaxPlus; plus++)
    {    lgamma_component(c,lett,plus) = lgamma_tmp;
    	lgamma_tmp += log(z+plus);
    }
    
    // Restore BetaComp
    LogBetaComp[c] += lgamma_component(c,lett);

    double diff = z - component(c)[lett];
    component(c)[lett] = z;

    if (SumComponents[c]>0) 
    {	LogBetaComp[c] += LogGamma(SumComponents[c]);
    }
    SumComponents[c] += diff;
    if (SumComponents[c] < 1e-4 || (-1e-5 *diff < z && z < 1e-5*diff))
    {    // possible accumulated round-off error, re-sum
    	double sum=0., log_beta=0.;
//	cerr << "Setting component " << c << " letter " << lett
//		<< " triggers re-computing sum\n" 
//		<< "Old SumComponents= " << SumComponents[c]
//		<< " LogBetaComp= " << LogBetaComp[c]
//		<< "\n"
//		<< flush;
	for (int i=alphabet_size()-1; i>=0; i--)
	{   double comp = component(c)[i];
	    sum += comp;
//	    cerr << i << " sum= " << sum << "\n" << flush;
	    // don't add in the uninitialized lgamma_components
	    if (comp > 0)
		log_beta += lgamma_component(c,i);
//	    cerr << i << " log_beta= " << log_beta << "\n" << flush;
	}
//	cerr << "New SumComponents is " << sum << "\n" << flush;
	SumComponents[c] = sum;
	LogBetaComp[c] = log_beta;
    }
    
    LogGamma_derivs(SumComponents[c], lgamma_tmp, Lgamma1Sum[c],
		Lgamma2Sum[c]); 
    LogBetaComp[c] -= lgamma_tmp;
}

void DirichletReg::scale_component(int i, float multiplier)
{
    for(int letter= alphabet_size()-1; letter>=0; letter--)
    	set_component(i, letter, multiplier*component(i)[letter]);
}






void DirichletReg::delete_component(int i)
{
    assert(0<=i && i < NumComp);
//    cerr << " Deleting component " << i << " from " << name() << " .\n";
    if (i < NumComp)
    {   const float *comp = component(NumComp -1);
    	for (int letter= alphabet_size()-1; letter>=0; letter--)
    	    set_component(i, letter, comp[letter]);
	set_mixture(i, MixtureCoeff[NumComp-1]);
    }
    NumComp --;
}

int DirichletReg::num_parameters(void) const 	
{   return FreezeComp? NumComp:
	FreezeMix? (alphabet_size()*NumComp):
		(alphabet_size()+1)*NumComp;
}
float DirichletReg::min_parameter(int i) const 
{    return 0;	// let things disappear!
}
float DirichletReg::max_parameter(int i) const
{   if (FreezeComp || 
	    (!FreezeMix && i % (alphabet_size()+1) == alphabet_size()))
	return 1.e4;	// max mixture coeff
		// allow some room--normalize() will rescale
    else return 1.e5;	// max pseudocount
}

float DirichletReg::parameter(int i) const
{   if (FreezeComp)
	return MixtureCoeff[i];
    int row_size = alphabet_size()+ (FreezeMix?0:1);
    int which_dist=i/row_size;
    int subscript=i-which_dist*row_size;
    return subscript==alphabet_size()?
		MixtureCoeff[which_dist]:
		component(which_dist)[subscript];
}


void DirichletReg::set_parameter(int i, float p)
{   
    assert(p >= min_parameter(i));
    assert(p <= max_parameter(i));
    if (FreezeComp)
    {	MixtureCoeff[i] = p;
	return;
    }
    int row_size = alphabet_size()+ (FreezeMix?0:1);
    int which_dist=i/row_size;
    int subscript=i-which_dist*row_size;
    if (subscript==alphabet_size())
	 MixtureCoeff[which_dist]= p;
    else set_component(which_dist,subscript,p);
}

// Make mixture coefficients sum to 1.
// Make components have no element less than 1.e-6.
// Recompute SumComponents and LogBeta to avoid roundoff errors.
void DirichletReg::normalize(void)
{   double sum= 0.;
    int c;
    for (c=NumComp-1; c>=0; c--)
	sum += MixtureCoeff[c];
    for (c=NumComp-1; c>=0; c--)
    {	MixtureCoeff[c] /= sum;
	float min_of_component=1.e30;
	int i;
	for (i= alphabet_size()-1; i>=0; i--)
	{   min_of_component= min(component(c)[i], min_of_component);
	}
	if (min_of_component < 1.e-6)
	{   float offset = 1e-5 - min_of_component;
	    for (i= alphabet_size()-1; i>=0; i--)
		set_component(c,i, component(c)[i]+offset);
    	}
    	double log_beta=0.;
	double sum_c=0;
	for (i=alphabet_size()-1; i>=0; i--)
	{   sum_c += component(c)[i];
	    log_beta += lgamma_component(c,i);
	}
	SumComponents[c] = sum_c;
	double lgamma_tmp;
	LogGamma_derivs(sum_c, lgamma_tmp, Lgamma1Sum[c],
		    Lgamma2Sum[c]); 
	log_beta -= lgamma_tmp;
	LogBetaComp[c] = log_beta;


    }
}


static int ReadAlphaChar(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{  
    int NumChar;
    in >> NumChar;
    if (NumChar != change->alphabet_size())
    {    cerr << "Number of characters " << NumChar 
	    << " from " << self->name() << " command, doesn't match "
	    << change->alphabet_size() << " for " 
	    <<  *(change->alphabet_tuple()) << "\n" << flush;
    }	
    return 1;
}

static int ReadNumDistr(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{  
    int NumDistr;
    in >> NumDistr;
    dynamic_cast<DirichletReg*>(change)->alloc(NumDistr+5);	// leave some room to grow
    return 1;
}

static int ReadNumber(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{  
    DirichletReg *g = dynamic_cast<DirichletReg *>(change);
    int Number;
    in >> Number;
    if (Number != g->num_components())
    {	cerr << "Warning: reading " << change->type()->name()
	    << " " << g->name() << " " << self->name()
	    << " = " << Number << ", but "
	    << g->num_components() << " already defined.\n"
	    << flush;
    }
    return 1;
}

static float MixtureCoeffSaved= -1.0;	
	// location for saving mixture coeff between Mixture and Alpha
	//	commands
	// -1 is sentinel to indicate no mixture coeff read
static int ReadMixture(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{
    in >> MixtureCoeffSaved;
    return 1;
}


static int ReadAlpha(istream &in, 
		Regularizer *change,
		RegInputCommand* self)
{
    DirichletReg *g = dynamic_cast<DirichletReg *>(change);
    if (MixtureCoeffSaved <=0)
    {	cerr << "Warning: reading " << change->type()->name()
	    << " " << g->name() << " for component " 
	    << g->num_components() 
	    << 	" need Mixture command before " << self->name()
	    <<  "\n";
	cerr << "Assuming Mixture=1\n";
	MixtureCoeffSaved = 1;
    }
    
    double sum_alpha;	
    float *comp = new float[g->alphabet_size()];
    const int* order=g->input_order();
    
    in >> sum_alpha;
    for (int i=0; i<g->alphabet_size(); i ++)
    {   in >> comp[order? order[i]: i];
    }
    g->AddComponent(MixtureCoeffSaved, comp);
    
    double sum_c = g->sum_component(g->num_components()-1);
    if (kk::abs(sum_c - sum_alpha) > 1.e-5* sum_c)
    {	cerr << "Warning: reading " << change->type()->name()
	    << " " << g->name() << " for component " 
	    << g->num_components()-1 << "\n"
	    << "   first argument of " << self->name()
	    << " was " << sum_alpha << " instead of " << sum_c
	    << "\n";
    }
    MixtureCoeffSaved = -1;
    delete [] comp;
    return 1;
}


void DirichletReg::init_command_table()
{   assert(!CommandTable);
    CommandTable = new NameToPtr(13);
    CommandTable->ignore_case();
    CommandTable->AddName(new
		RegInputCommand("AlphaChar", ReadAlphaChar));
    CommandTable->AddName(new
		RegInputCommand("NumDistr", ReadNumDistr));
    CommandTable->AddName(new
		RegInputCommand("Number", ReadNumber));
    CommandTable->AddName(new
		RegInputCommand("Mixture", ReadMixture));
    CommandTable->AddName(new
		RegInputCommand("Alpha", ReadAlpha));
    CommandTable->AddName(new RegInputCommand("FullUpdate", ReadComment));
    CommandTable->AddName(new RegInputCommand("QUpdate", ReadComment));
    CommandTable->AddName(new RegInputCommand("StructID", ReadComment));
}

void DirichletReg::AddComponent(float MixCoeff,const float*comp)
{
//    cerr << "Adding component to " << name() 
//	<< " with MixCoeff= " << MixCoeff << "\n" << flush;
    if (NumComp>=AllocComp)
    {	// make room for 10 more
    	float* old_comp = Components;
	Components=0;
	float* old_mix = MixtureCoeff;
	MixtureCoeff=0;
	int old_num=NumComp;
	delete_all();
	alloc(old_num+10);
	for (int i=0; i<old_num; i++)
	    AddComponent(old_mix[i], old_comp + alphabet_size()*i);
	delete [] old_mix;
	delete [] old_comp;
    }
    int c=NumComp;
    MixtureCoeff[c] = MixCoeff;
    SumComponents[c] =0.;
    LogBetaComp[c] = 0;
    int lett;
    for (lett=alphabet_size()-1; lett>=0; lett--)
    {
	component(c)[lett] = 0.;
    }
    for (lett=alphabet_size()-1; lett>=0; lett--)
    	set_component(c,lett,comp[lett]);
    NumComp++;
//    cerr << "Sum of component=" << SumComponents[c] << "\n" << flush;
}

void DirichletReg::print_ordered_component(ostream &out, int c) const
{
    // Get estimated background probabilities for
    //	adding comments
    float *probs = new float[alphabet_size()];
    int letter;
    for (letter=0; letter< alphabet_size(); letter++)
    	probs[letter] = NumComp==1? 1.: 0.;
    if (NumComp >1)
    {    for (int d=0; d<NumComp; d++)
	{   const float *comp = component(d);
	    for (letter=0; letter< alphabet_size(); letter++)
	    {    probs[letter] += MixtureCoeff[d] * 
				comp[letter]/ SumComponents[d];
	    }
	}
    }
    double sum=0.;
    for (letter=0; letter< alphabet_size(); letter++)
    	sum += probs[letter];
    for (letter=0; letter< alphabet_size(); letter++)
	probs[letter] /= sum;    
    
    float *cost = new float[alphabet_size()];	// -log2 component/probs
    const float *comp = component(c);
    for (letter=0; letter< alphabet_size(); letter++)
 	    cost[letter] = -log2(comp[letter] / 
		    	(SumComponents[c] * probs[letter])); 
    
    int * low_cost_first = OrderSort(cost, alphabet_size());
    for (int ord=0; ord< alphabet_size(); ord++)
    {   alphabet_tuple()->print_unindex(out,low_cost_first[ord]);
	if (ord+1 < alphabet_size())
	{   // put out separator symbol(s)
	    out << " "; 
	    float f = cost[low_cost_first[ord]];
	    float fnext = cost[low_cost_first[ord+1]];
	    int trunc_f = static_cast<int>(f);
	    int trunc_fnext = static_cast<int>(fnext);
	    for (int i=trunc_f ;
			i <trunc_fnext || (i==trunc_fnext && i==0) ;
			 i++)
	    {   if (i==0 && f*fnext<0)
		{   out << "><";
		    trunc_fnext++;	// zero crossing produces extra symbol
		}
		else if (i<trunc_fnext)
		    out << ",";
	    }
	    if (trunc_fnext > trunc_f)
		out << " ";
	}
    }
    delete [] low_cost_first;
    delete [] cost;
    delete [] probs;
}

void DirichletReg::write_knowing_type(ostream &out) const
{    
    Regularizer::write_knowing_type(out);
    out << "AlphaChar= " << alphabet_size() << "\n"
    	<< "NumDistr= " << NumComp <<"\n";
    
    
    for (int c=0; c<NumComp; c++)
    {	
	out << "\n"
	    << "Number= " << c << "\n"
	    << "Mixture= " << MixtureCoeff[c] << "\n";
	
	   
	out << "Alpha= " << SumComponents[c] ;
	const float *comp = component(c);
	for (int letter=0; letter< alphabet_size(); letter++)
	    out << " " << comp[letter];
    	out << "\n";

	out << "Comment= ";
	print_ordered_component(out, c);
	out << "\n" << flush;
    }
}

void DirichletReg::get_moments( const float* TrainCounts,
	     double* ex_prob, double* ex_prob2)
{
    if (NumComp==0)
    {   double share = 1./alphabet_size();
	if (ex_prob)
	{   for (int letter = alphabet_size()-1; letter >= 0; letter--) 
	        ex_prob[letter] = share ;
	}
	if (ex_prob2)
	{   double share2= share*share;
	    for (int letter = alphabet_size()-1; letter >= 0; letter--) 
	        ex_prob2[letter] = share2 ;
	}
    	return;
    }

    component_probs(TrainCounts, SumCounts, CompProbs);
    
    int letter;
    for (letter = alphabet_size()-1; letter >= 0; letter--) 
	X[letter] =  0.;
    if (ex_prob2)
    {   for (letter = alphabet_size()-1; letter >= 0; letter--) 
	    ex_prob2[letter] =  0.;
    }
    
    for (int c=NumComp-1; c>=0; c--)
    {   float *comp=component(c);
	double weight= CompProbs[c]
			/ (SumCounts+SumComponents[c]);
	for (letter = alphabet_size()-1; letter >= 0; letter--) 
		X[letter] +=  weight * (TrainCounts[letter]+comp[letter]);
	if (ex_prob2)
	{   weight /= (SumCounts+SumComponents[c]+1);
	    for (letter = alphabet_size()-1; letter >= 0; letter--) 
	    {	double num = (TrainCounts[letter]+comp[letter]);
		ex_prob2[letter] +=  weight * num * (num+1);
	    }
	}
    }
    if (ex_prob)
    {	for (letter = alphabet_size()-1; letter >= 0; letter--) 
	    ex_prob[letter] = X[letter];
    }
}

void DirichletReg::get_modified_counts( const float* TrainCounts,
	     float* probs)
{
    if (NumComp==0)
    {   double share = 1./alphabet_size();
	for (int letter = alphabet_size()-1; letter >= 0; letter--) 
	    probs[letter] = X[letter] =  share;
    	return;
    }

    component_probs(TrainCounts, SumCounts, CompProbs);
    
    int alph_size_minus_1=alphabet_size()-1;
    
    int letter;
    for (letter = alphabet_size()-1; letter >= 0; letter--) 
	X[letter] =  0.;

    for (int c=NumComp-1; c>=0; c--)
    {  	register double weight= CompProbs[c]
			/ (SumCounts+SumComponents[c]);
        register const float*comp = component(c);

	for (letter = alph_size_minus_1; letter >= 0; letter--) 
	{   X[letter] += weight*(TrainCounts[letter] + comp[letter]);
	}
    }
    for (letter = alphabet_size()-1; letter >= 0; letter--) 
	probs[letter] = X[letter];
}



// Fill in Lgamma1_Counts so that
//		lgamma1_counts(c,k)=lgamma1(Counts[k] + component(c)[k])
//	and Lgamma1_SumCounts[c] = lgamma1(SumCounts + SumComponents[c])
void DirichletReg::compute_lgamma1(const float *Counts)
{
    for (int c=num_components()-1; c>=0; c--)
    {   Lgamma1_SumCounts[c] = LogGamma_1(SumCounts + SumComponents[c]);
	for (int k=alphabet_size()-1; k>=0; k--)
	{   if (Counts[k]>0)
		lgamma1_counts(c,k) = LogGamma_1(Counts[k] +component(c)[k]);
	    else lgamma1_counts(c,k) = lgamma1_comp(c)[k];
	}
    }
    Lgamma1CountsInvalid=0;
}

// Fill in Lgamma1_Counts so that
//		lgamma1_counts(c,k)=lgamma1(Counts[k] + component(c)[k])
//	and Lgamma1_SumCounts[c] = lgamma1(SumCounts + SumComponents[c])
// also 
//		lgamma2_counts(c,k)=lgamma2(Counts[k] + component(c)[k])
//	and Lgamma2_SumCounts[c] = lgamma2(SumCounts + SumComponents[c])
void DirichletReg::compute_lgamma2(const float *Counts)
{
    for (int c=num_components()-1; c>=0; c--)
    {   double lg;
        LogGamma_derivs(SumCounts + SumComponents[c], lg,
		Lgamma1_SumCounts[c], 	Lgamma2_SumCounts[c]);
	for (int k=alphabet_size()-1; k>=0; k--)
	{   if (Counts[k]>0)
		LogGamma_derivs(Counts[k] +component(c)[k], lg, 
			lgamma1_counts(c,k), lgamma2_counts(c,k));
	    else
	    {	lgamma1_counts(c,k) = lgamma1_comp(c)[k];
		lgamma2_counts(c,k) = lgamma2_comp(c)[k];
	    }
	}
    }
    Lgamma1CountsInvalid=0;
    Lgamma2CountsInvalid=0;
}


void DirichletReg::partials1(float *part1, int i, const float* Counts)
{   
    // X_i = sum_d Prob(d) (Counts[i] + component(d)[i])/
    //				(SumCounts+SumComponents[d]))
    //
    // X'_i = sum_d Prob(d) (component'(d)[i] / (SumCounts+SumComponents[d])
    //			- SumComponents'[d](Counts[i] + component(d)[i])/
    //					    (SumCounts+SumComponents[d])^2)
    //		   + Prob'(d) (Counts[i] + component(d)[i])/
    //				(SumCounts+SumComponents[d]))
    //	
    // For derivative with respect to MixtureCoeff[c],
    // X'_i = sum_d Prob'(d) (Counts[i] + component(d)[i])
    //			/ (SumCounts+SumComponents[d])
    //
    // For derivative w.r.t component(c)[j],
    // X'_i = Prob(c) (delta(i==j) / (SumCounts+SumComponents[d]) 
    //			- (Counts[i] + component(c)[i])/
    //					    (SumCounts+SumComponents[c])^2)
    //	    + sum_d Prob'(d) (Counts[i] + component(d)[i])/
    //				(SumCounts+SumComponents[d]))

    // Note: derivative of Prob(d) w.r.t. component(c,k) or MixtureCoeff[c]
    // can be derived as follows:
    //	(letting q_d = MixtureCoeff[d] and alpha_d=component(d))
    // Prob(d) = q_d B(alpha_d+s)/B(alpha_d) / 
    //			sum_e q_e (B(alpha_e+s)/B(alpha_e))
    //	= exp(ln q_d + log_beta_over_beta(d)) /
    //			sum_e exp(ln q_e + log_beta_over_beta(e))
    // Letting u_e = ln q_e + log_beta_over_beta(e), 
    // we have Prob(d)=exp(u_d)/ sum_e exp(u_e), then
    //	Prob'(d) = exp(u_d) u'_d / sum_e exp(u_e)
    //		- exp(u_d) (sum_e exp(u_e) u'_e) / (sum_e exp(u_e))^2 
    //	   = Prob(d) (u'_d - (sum_e exp(u_e) u'_e) / sum_e exp(u_e) )
    //	   = Prob(d) (u'_d - sum_e Prob(e) u'_e)
    //
    //	Notice that u'_e = 0 for e!= c (where the derivative is with respect
    //	to Mixture coefficient q_c or any part of component c),
    // so Prob'(d) = Prob(d) (delta(d==c) -Prob(c)) u'_c 
    
    //	sum_d Prob'(d) (Counts[i] + component(d)[i]) 
    //				/ (SumCounts+SumComponents[d])
    //  = u'(c) Prob(c) (Counts[i] + component(c)[i]) 
    //				/ (SumCounts+SumComponents[c])
    //		- u'(c) Prob(c) sum_d Prob(d) (Counts[i] + component(d)[i]) 
    //				/ (SumCounts+SumComponents[d])
    // = u'(c) Prob(c) [(Counts[i] + component(c)[i]) 
    //				/ (SumCounts+SumComponents[c])
    //			- X_i]


    // u'_c = 1/q_c when taking partial with respect to mixture coeff q_c
    // u'_c = lgamma1(Counts[j] + component(c)[j]) - lgamma1_comp(c)[j]
    //		-lgamma1(SumCounts+SumComponents[c]) + Lgamma1Sum
    //		for partial with respect to component(c)[j]

    // So we can simplify the derivatives of X'_i:
    //
    // For derivative with respect to MixtureCoeff[c],
    // X'_i = Prob(c) /MixtureCoeff[c]
    //		 [(Counts[i] + component(c)[i]) / (SumCounts+SumComponents[c])
    //		- X_i]
    //
    // For derivative w.r.t component(c)[j],
    // X'_i = Prob(c) (delta(i==j) / (SumCounts+SumComponents[c])
    //			- (Counts[i] + component(c)[i])/
    //					    (SumCounts+SumComponents[c])^2)
    //		+Prob(c) 
    //		* (lgamma1(Counts[j] + component(c)[j]) - lgamma1_comp(c)[j]
    //			-lgamma1(SumCounts+SumComponents[c]) + Lgamma1Sum)
    //		* [(Counts[i] + component(c)[i])/ (SumCounts+SumComponents[c])
    //		-X_i]
    
    
//    if (Lgamma1CountsInvalid && num_components()>1)
    if (Lgamma1CountsInvalid)
	compute_lgamma1(Counts);
    int as1=alphabet_size() + (FreezeMix? 0:1);
    for (int c=NumComp-1; c>=0; c--)
    {	    
	double counts_plus=(Counts[i] + component(c)[i]);
	double sum_counts_plus = (SumCounts+SumComponents[c]);
	double normalized_post = counts_plus / sum_counts_plus;
	double this_minus_x = normalized_post - X[i];
	// force this_minus_x to be zero if it is not just because of
	// bad round-off error
	if (kk::abs(this_minus_x) < 1.0e-8* normalized_post)
		this_minus_x=0;
		
	double prob_c_this_minus_x = CompProbs[c]*this_minus_x;
	
	// first handle the derivatives w.r.t. the mixture coefficients.
	if (FreezeComp)
	{   part1[c] = prob_c_this_minus_x / MixtureCoeff[c];
	    continue;
	}
	if (!FreezeMix)
	    part1[c*as1+alphabet_size()] =prob_c_this_minus_x/MixtureCoeff[c];

	double prob_norm_over_sum =
		CompProbs[c]*normalized_post/sum_counts_plus; 
	for (int k=alphabet_size()-1; k>=0; k--)
	{   double u_prime_c= lgamma1_counts(c,k) - lgamma1_comp(c)[k]
				-Lgamma1_SumCounts[c] + Lgamma1Sum[c];
	    part1[c*as1+k] = prob_c_this_minus_x *u_prime_c 
		- prob_norm_over_sum;
	}
	part1[c*as1+i] += CompProbs[c]/sum_counts_plus;
    }
    
}

void DirichletReg::partials2(float *part1, float*part2, int i, 
	const float* Counts)
{   
    // Letting u_c = ln q_c + log_beta_over_beta(c), 
    // counts_plus= (Counts[i] + component(c)[i])
    // sum_counts_plus = (SumCounts+SumComponents[c])
    // normalized_post = counts_plus/sum_counts_plus
    
    
    // For derivative w.r.t. MixtureCoeff[c] or component(c)[j]
    // Prob'(d) = Prob(d) u'(c) (delta(d==c) - Prob(c))
    // In particular, Prob'(c)=Prob(c)(1-Prob(c))u'(c)
    
    // for derivative w.r.t. MixtureCoeff[c] 
    // u'(c) = 1/MixtureCoeff[c], and u''(c) = -1/MixtureCoeff[c]^2
    // normalized_post' = 0
    //
    // X'_i = Prob(c) u'(c) (normalized_post-X_i)
    // X''_i = (1-2 Prob(c)) u'(c) X'_i + Prob(c) u''(c) (normalized_post-X_i)
    //		= -2 Prob(c) X'_i / MixtureCoeff[c]
    
    // for derivative w.r.t. component[c](j)
    // normalized_post' = (delta(i==j) - normalized_post)/sum_counts_plus
    // normalized_post'' = -2 (delta(i==j) - normalized_post)/sum_counts_plus^2
    //
    // X'_i = Prob(c) 
    //		[ u'(c)  (normalized_post-X_i) +
    //			(delta(i==j) - normalized_post)/sum_counts_plus]
    // X''_i = Prob'(c)/Prob(c) X'_i
    //		+ Prob(c) u''(c)  (normalized_post-X_i) 
    //		+ Prob(c) u'(c)(delta(i==j) - normalized_post)/sum_counts_plus)
    //		- Prob(c) u'(c) X'_i
    //		-2 Prob(c) (delta(i==j)-normalized_post)/ sum_counts_plus ^2
    //
    // X''_i = (1-2 Prob(c))u'(c) X'_i
    //		+ Prob(c) u''(c)  (normalized_post-X_i) 
    //		+ Prob(c) (u'(c)-2/sum_counts_plus) 
    //			(delta(i==j) - normalized_post)/sum_counts_plus)
    //
    // u'(c) = lgamma1(counts_plus) - lgamma1_comp(c)[j]
    //			-lgamma1(sum_counts_plus) + Lgamma1Sum
    // u''(c) = lgamma2(counts_plus) - lgamma2_comp(c)[j]
    //			-lgamma2(sum_counts_plus) + Lgamma2Sum

//    if (Lgamma2CountsInvalid && num_components()>1)
    if (Lgamma2CountsInvalid)
	compute_lgamma2(Counts);
    int as1=alphabet_size() + (FreezeMix? 0:1);
    double *u_prime_c = new double[alphabet_size()];
    for (int c=NumComp-1; c>=0; c--)
    {	    
	double counts_plus=(Counts[i] + component(c)[i]);
	double sum_counts_plus = (SumCounts+SumComponents[c]);
	double normalized_post = counts_plus / sum_counts_plus;
	double this_minus_x = normalized_post - X[i];
	// force this_minus_x to be zero if it is not just because of
	// bad round-off error
	if (kk::abs(this_minus_x) < 1.0e-8* normalized_post)
		this_minus_x=0;
	double prob_c_this_minus_x = CompProbs[c]*this_minus_x;
	
	// first handle the derivatives w.r.t. the mixture coefficients.
	if (FreezeComp)
	{   part1[c] = prob_c_this_minus_x / MixtureCoeff[c];
	    part2[c]= -2 *part1[c] * ( CompProbs[c]/MixtureCoeff[c]);
	    continue;
	}
	if (!FreezeMix)
	{   part1[c*as1+alphabet_size()] =prob_c_this_minus_x/MixtureCoeff[c];
	    part2[c*as1+alphabet_size()] = part1[c*as1+alphabet_size()] *
		-2 * ( CompProbs[c]/MixtureCoeff[c]);
	}
	double prob_norm_over_sum =
		CompProbs[c]*normalized_post/sum_counts_plus; 
	double Lgamma1Sums = Lgamma1Sum[c]-Lgamma1_SumCounts[c] ;
	int k;
	for (k=alphabet_size()-1; k>=0; k--)
	{   u_prime_c[k] = lgamma1_counts(c,k)-lgamma1_comp(c)[k] +Lgamma1Sums;
	    part1[c*as1+k] = prob_c_this_minus_x *u_prime_c[k]
		- prob_norm_over_sum;
	}
	part1[c*as1+i] += CompProbs[c]/sum_counts_plus;

	double Lgamma2Sums = Lgamma2Sum[c]-Lgamma2_SumCounts[c] ;
	for (k=alphabet_size()-1; k>=0; k--)
	{   double u_2_c = lgamma2_counts(c,k)-lgamma2_comp(c)[k] +Lgamma2Sums;
	    part2[c*as1+k] = (1-2*CompProbs[c])*part1[c*as1+k] *u_prime_c[k]
		+ u_2_c*prob_c_this_minus_x
		- prob_norm_over_sum*(u_prime_c[k] - 2/sum_counts_plus);
	}
	part2[c*as1+i] += CompProbs[c]*(u_prime_c[i] - 2/sum_counts_plus)
		/sum_counts_plus;

    }
    delete [] u_prime_c;
}

// compute Beta(counts+z)/Beta(z) for each component
// and return the max value for normalizing
//
double DirichletReg::log_beta_over_beta(const float* TrainCounts,
	double &SumTrainCounts,
	double* log_Beta_Beta) const
{
    double  max_log_Beta_Beta= -1.e99;
    int alph_size = alphabet_size();
    
    int gamma_component_row_size = alph_size*(MaxPlus+1);
    double *lgamma_component_base = LogGammaComp +
		NumComp*gamma_component_row_size; 
    
    for (short int c=NumComp-1; c>=0; c--)
    {   lgamma_component_base -= gamma_component_row_size;
	register double log_beta_beta = 
		-LogGamma(SumTrainCounts+SumComponents[c]); 
	
	register const float * comp = component(c);
	for(int letter=alph_size-1; letter>=0; letter--)
	{   int icount;
	    double tc = TrainCounts[letter];
	    if (is_small_int(tc, icount))
	    {    log_beta_beta +=
		     lgamma_component_base[(MaxPlus+1)*letter+icount];
	    }
	    else
	    {   log_beta_beta += LogGamma(tc + comp[letter]);
	    }
	}
	log_beta_beta -=  LogBetaComp[c];
	if (log_beta_beta > max_log_Beta_Beta)
		max_log_Beta_Beta = log_beta_beta;
	log_Beta_Beta[c] = log_beta_beta;
    }
    return max_log_Beta_Beta;
}

// Return the probability for each component given the sample.
// comp_probs should be pre-allocated with num_components() floats.
//
//	probability of component c computed by normalizing
//		mixture_coeff(c)  * Prob(TrainCounts | component(c))
//	(TrainCounts treated as the summary of a list of characters---no
//	correction for permutations)
//	
//	optionally (in *log_sum)
//	return	ln sum_c mixture_coeff(c) * Prob(TrainCounts | component(c))
//
void DirichletReg::component_probs(const float* TrainCounts,
		double& SumTrainCounts,
		double* probs,
		double* log_sum) 
{
    Lgamma1CountsInvalid=1;	// the derivatives are not computed yet
    Lgamma2CountsInvalid=1;
    SumTrainCounts=0.;
    for (int letter=alphabet_size()-1; letter>=0; letter--)
	SumTrainCounts += TrainCounts[letter];

    if (NumComp==1 && !log_sum) 
    {   probs[0] =1;	// speed up for single component
	return ;
    }
    
    double *sw = new double[NumComp];   
    // fill sw[c] with log Beta(comp+counts)/Beta(comp)
    double max_sw = log_beta_over_beta(TrainCounts, SumTrainCounts, sw);
    
    // fill comp_prob[c] 
    //		with a scaled version of Beta(comp+counts)/Beta(comp)
    //		multiplied by the existing mixture coefficient
    double sum=0;
    int c;
    for (c=NumComp-1; c>=0; c--)
    {  	probs[c] = exp(sw[c]-max_sw) * MixtureCoeff[c];
	sum += probs[c];
    }
    
    // normalize to sum to one
    double scale = 1./sum;
    for (c=NumComp-1; c>=0; c--)
	probs[c] *= scale;   
    if (log_sum) (*log_sum) =  max_sw + log(sum);

    delete [] sw;
}

// return the natural log of the probability of the count vector
// along with the partial derivatives of the natural log (if arrays provided)
// Note: second derivative only taken with respect to same parameter twice,
// not off-diagonal.
//
double DirichletReg::log_probability(const float* TrainCounts,
	float *deriv1, float*deriv2) 
{
    double log_prob;
    component_probs(TrainCounts,SumCounts, CompProbs, &log_prob);
    
    if (SumCounts<=0)	
    {	if (deriv1) 
	{    for (int i=num_parameters()-1; i>=0; i--)
		deriv1[i] = 0;
	}
	if (deriv2) 
	{    for (int i=num_parameters()-2; i>=0; i--)
		deriv2[i] = 0;
	}
	return 0;	// Probability of the 0 vector is always 1
    }

    int c;
    double sum_mix = 0.;	// sum of mixture coefficients
    for (c=num_components()-1; c>=0; c--)
	sum_mix += MixtureCoeff[c];    
    log_prob -= log(sum_mix);
    
    if (!deriv1 && !deriv2)     return log_prob;
    
    if (deriv2) compute_lgamma2(TrainCounts);
    else if (deriv1) compute_lgamma1(TrainCounts);
    
    int as1=alphabet_size() + (FreezeMix? 0:1);
    
    for (c=num_components()-1; c>=0; c--)
    {   
	// which parameter is the mixture coefficient
	int mix_subscript = FreezeComp? c : (as1*c + alphabet_size());
	if (!FreezeMix)
	{   double mult =  CompProbs[c]/MixtureCoeff[c];
	    if (deriv1) deriv1[mix_subscript] = mult - 1/sum_mix;
	    if (deriv2) 
		deriv2[mix_subscript] = 1/(sum_mix*sum_mix) -mult*mult;
	}
	if (FreezeComp) continue;
	
	float *deriv1c = deriv1+c*as1;	// beginning of row of param
	float *deriv2c = deriv2+c*as1;
	double Lgamma1Sums = Lgamma1Sum[c]-Lgamma1_SumCounts[c] ;
	double Lgamma2Sums = deriv2? (Lgamma2Sum[c]-Lgamma2_SumCounts[c]) : 0;

	for (int j=alphabet_size()-1; j>=0; j--)
	{   double lgm1 = lgamma1_counts(c,j)-lgamma1_comp(c)[j] +Lgamma1Sums;
	    if (deriv1) deriv1c[j] = CompProbs[c] * lgm1;
	    if (deriv2)
	    {   double lgm2 = lgamma2_counts(c,j)-lgamma2_comp(c)[j] 
			+Lgamma2Sums;
		deriv2c[j] = CompProbs[c] * (lgm1*lgm1*(1-CompProbs[c])+lgm2);
	    }
	}
    }    
    
    return log_prob;
}

int DirichletReg::verify_log_prob_partials(ostream &logfile,
	const float*TrainCounts, float tolerance)
{
    int ret=1;	// return value, set to 0 if bad partial found
    float *part1 = new float[num_parameters()];
    float *part2 = new float[num_parameters()];
    
    double log_prob = log_probability(TrainCounts, part1, part2);
    for (int p=num_parameters()-1; p>=0; p--)
    {   float save=parameter(p);
	double delta = save< 0.1? kk::abs(save*0.1) : kk::abs(save*0.01);
	set_parameter(p, save+delta);
	double log_prob_plus = log_probability(TrainCounts);
	set_parameter(p, save-delta);
	double log_prob_minus = log_probability(TrainCounts);
	set_parameter(p, save);

        double diff= 0.5*(log_prob_plus - log_prob_minus);
	double loosened_tolerance = 
	    max(static_cast<double>(tolerance), 1.e-6 *kk::abs(log_prob / diff));
	double deriv = diff/delta;
	double abserr = kk::abs(deriv - part1[p]);
	if (abserr > loosened_tolerance
	    && kk::abs(abserr/ deriv) > loosened_tolerance)
	{    logfile << "First log_prob partial outside tolerance "
		    << loosened_tolerance
		    << "\n for "
		    << type()->name() << " " << name()
		    << "\n parameter " << p << " = " << save
		    << "\n";
	    logfile << "counts=";
	    for (int j=0; j<alphabet_size(); j++)
		logfile << " " << TrainCounts[j];
	    logfile << "\n";
	    logfile << "computed partial = " << part1[p]
		    << " but deriv_from_diff = " << deriv << "\n"
		    << " log_prob_plus= " << log_prob_plus
		    << " log_prob= " << log_prob
		    << " log_prob_minus= " << log_prob_minus
		    << "\n\n";
	    ret = 0;
	}

	double diff2 = log_prob_plus + log_prob_minus - 2*log_prob;
	loosened_tolerance = max(static_cast<double>(tolerance),
                                 1.e-6*kk::abs(log_prob/diff2));
	double deriv2 = diff2/(delta*delta);
	double abserr2 = kk::abs(deriv2 - part2[p]);
	if (abserr2 > loosened_tolerance
	    && kk::abs(abserr2/ deriv2) > loosened_tolerance)
	{    logfile << "Second log_prob partial outside tolerance "
		    << loosened_tolerance
		    << "\n for "
		    << type()->name() << " " << name()
		    << "\n parameter " << p << " = " << save
		    << "\n";
	    logfile << "counts=";
	    for (int j=0; j<alphabet_size(); j++)
		logfile << " " << TrainCounts[j];
	    logfile << "\n";
	    logfile << "computed partial2 = " << part2[p]
		    << " but deriv2_from_diff = " << deriv2 << "\n"
		    << " log_prob_plus= " << log_prob_plus
		    << " log_prob= " << log_prob
		    << " log_prob_minus= " << log_prob_minus
		    << "\n\n";
	    ret = 0;
	}

    }

    delete [] part2;
    delete [] part1;
    return ret;
}

// CHANGE LOG
// 7 March 1995 Kevin Karplus
//	Fixed Prob casts to work with new Prob.h file
//	Note: *= doesn't work with cast right-hand side any more.
// 1 June 1995 Kevin Karplus
//	Changed normalization to keep components away from 0
//	Changed parameters to be  q_c and z_c(i) * q_c 
//	Changed normalize() to eliminate very small components.
// 2 June 1995 Kevin Karplus
//	Removed the component deletion from normalize()---it was too
//		aggressive. 
//	Added comment to output that reports high and low probability
//		letters in each component.
// 9 June 1995 Kevin Karplus
//	Added Order= to input and output
//	Modified comment on components to be more descriptive
//	Modified normalize() to rescale extremely small components.
// 25 July 1995 Kevin Karplus
//	Corrected formula for Dirichlet mixture posterior
//	Added larger table of GammaComp, to allow faster computation
//	for small integer samples (also introduced MaxPlus).
//	Changed parameters back to being q_c and z_c(i).
//	Eliminated "real" derivatives (which are now more complicated),
//	and force use of "fake" derivatives that assume sample_weight
//	is independent of z_c.
// 2 August 1995 Kevin Karplus
//	Added copy constructor
// 20 Sept 1995 Kevin Karplus
//	Added component_probs(), and changed old "sample_weight"
//		field to "CompProbs"
// 23 Sept 1995 Kevin Karplus
//	Separated out the beta_over_beta routine,
//	added Probability(TrainCounts)
// 28 Sept 1995 Kevin Karplus
//	Added get_moments
// 6 Dec 1995 Kevin Karplus
//	Fixed initialization bug for max_log_Beta_Beta in beta_over_beta.  
// 8 Dec 1995 Kevin Karplus
//	Fixed uninitialized variable problems for void constructor,
//	by making alloc(0) clear the pointers.
// 31 Dec 1995 Kevin Karplus
//	Removed all use of Prob except return value from Probability,
//		using double of log_probability instead.  
//	Separated print_ordered_component out from write_knowing_type.
//	Restored lgamma1 and lgamma2 caches.
// 2 Jan 1996 Kevin Karplus
//	Fixed derivative w.r.t. mixture coeff, though derivative w.r.t.
//	components is still assuming that beta_over_beta is constant.
// 4 Jan 1996 Kevin Karplus
//	added caches for lgamma1 of counts+components
//	added flag to indicate caches invalidated by new component_probs
//	Now have true first partials
// 11 Jan 1996 Kevin Karplus
//	Now have true second partials
// 12 Jan 1996 Kevin Karplus
//	added log_probability, and defined Probability in terms of it.
// 16 Jan 1996 Kevin Karplus
//	made log_probability NOT correct for multinomials.
//	Added verify_log_prob_partials
// 30 Jan 1996 Kevin Karplus
//	Modified normalize to recompute SumComponents and LogBeta,
//	to control roundoff error.  (Also recompute LogBeta whenever
//	SumComponents recomputed in set_component).
// 7 Feb 1996 Kevin karplus
//	Added posterior_mixture()
// 7 March 1996 Kevin Karplus
//	Fixed constructors to copy whole AlphabetTuple
//	Modified print_ordered component to assume flat background when
//		only 1 component.
// 19 Oct 1996 Spencer Tu
//      Changed variable length arrays to dynamically allocated arrays in:
//          DirichletReg(const MLZReg &)
//          posterior_mixture
//          ReadAlpha
//          print_ordered_component
//          partials2
//          component_probs
//          verify_log_prob_partials
// 4 Dec 1996 Kevin Karplus
//	Fixed test for validity of Lgamma2Counts in partials2,
//		needs to be computed even if only one component,
//		since subsequent computation uses lgamma2_counts cache
// 13 Aug 2003 George Shackelford
//	Changed 'abs'  'kk::abs'
// 15 March 2004 Kevin Karplus
//	Fixed old-style casts
// 30 March 2004 Kevin Karplus
//	Fixed some ==0 tests to >=0.
// 9 April 2004 Kevin Karplus
//	Changed include because minmax.h and abs.h moved to Utilities/.
// Thu Jan 26 15:26:07 PST 2006 Kevin Karplus
//	Switch abs back to kk::abs, since g++ 3.4.4 gets confused by abs.
