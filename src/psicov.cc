/* PSICOV - Protein Sparse Inverse COVariance analysis program */

/* by David T. Jones August 2011 - Copyright (C) 2011 University College London */

/* This code is licensed under the terms of GNU General Public License v2.0 */

/* Version 1.10 - Last Edit 19/12/12 */

/* Modified by Tim Nugent to add dirichlet priors 01/06/13 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <fstream>
#include <ctype.h>
#include "DirichletReg.h"
#include "Alph.h"

#define FALSE 0
#define TRUE 1

#define SQR(x) ((x)*(x))
#define MAX(x,y) ((x)>(y)?(x):(y))
#define MIN(x,y) ((x)<(y)?(x):(y))

#define MAXSEQLEN 5000
#define MINEFSEQS (seqlen)

extern "C" void glasso_(int *, double *, double *, int *, int *, int *, int *, double *, int *, double *, double *, int *, double *, int *);

/* Dump a rude message to standard error and exit */
void fail(char *fmt, ...)
{
    va_list ap;
    
    va_start(ap, fmt) ;
    fprintf(stderr, "*** ");
    vfprintf(stderr, fmt, ap);
    fputc('\n', stderr);
    
    exit(-1);
}

/* Convert AA letter to numeric code (0-21) */
/*
int  aanum(int ch)
{
    const static int aacvs[] =
    {
	999, 0, 3, 4, 3, 6, 13, 7, 8, 9, 21, 11, 10, 12, 2,
	21, 14, 5, 1, 15, 16, 21, 19, 17, 21, 18, 6
    };

    return (isalpha(ch) ? aacvs[ch & 31] : 20);
}
*/

/* Convert AA letter to numeric code (0-21) - matches dirichlet encoding */
int aanum(char ch){

	switch(ch){
		case 'A': return 0;
		case 'C': return 1;
		case 'D': return 2;
		case 'E': return 3;
		case 'F': return 4;
		case 'G': return 5;
		case 'H': return 6;
		case 'I': return 7;
		case 'K': return 8;
		case 'L': return 9;
		case 'M': return 10;
		case 'N': return 11;
		case 'P': return 12;
		case 'Q': return 13;
		case 'R': return 14;
		case 'S': return 15;
		case 'T': return 16;
		case 'V': return 17;
		case 'W': return 18;
		case 'Y': return 19;
		default: return 20;
	}
}

/* Allocate matrix */
/*
void *allocmat(int rows, int columns, int size)
{
    int i;
    void  **p, *rp;

    rp = malloc(rows * sizeof(void *) + sizeof(int));

    if (rp == NULL)
	fail("allocmat: malloc [] failed!");

    *((int *)rp) = rows;

    p = (void**)rp + sizeof(int);

    for (i = 0; i < rows; i++)
	if ((p[i] = calloc(columns, size)) == NULL)
	    fail("allocmat: malloc [][] failed!");

    return p;
}

*/

/* Free matrix */
/*
void freemat(void *rp)
{
    int             rows;
    void **p = (void**)rp;

    rows = *((int *)(rp - sizeof(int)));

    while (rows--)
	free(p[rows]);

    free(rp - sizeof(int));
}
*/

/* Allocate matrix */
double** allocmat(int rows, int columns){
	double** p = new double*[rows];
	if(p == NULL){
		fail("allocmat failed!");
	}
	for(int i = 0; i < rows; i++){
	   p[i] = new double[columns];
	   if(p[i] == NULL){
	   		fail("allocmat failed!");
	   }
	}
	return p;
}

/* Free matrix */
void freemat(double** p, int rows, int columns){
	for(int i = 0; i < rows; i++){
	   delete [] p[i];
	}   
	delete [] p;
}

/* Allocate Fortran compatible square matrix */
void *f_matrix_calloc(int ndim, int size)
{
    void *p;
    
    p = calloc(ndim * ndim, size);
    
    if (p == NULL)
	fail("f_matrix_calloc: calloc failed!");

    return p;
}

/* Allocate vector */
void *allocvec(int columns, int size)
{
    void          *p;

    p = calloc(columns, size);

    if (p == NULL)
	fail("allocvec: calloc failed!");

    return p;
}

/* Perform Cholesky decomposition on matrix */
int test_cholesky(double *a, const int n) 
{
    int i, j, k;
    double sum;
    static double *diag;

    if (diag == NULL)
	diag = (double*)allocvec(n, sizeof(double));

    for (i=0; i<n; i++)
    {
	for (j=i; j<n; j++)
	{
	    sum = a[i*n+j];

	    for (k=i-1; k >= 0; k--)
		sum -= a[i*n+k]*a[j*n+k];

	    if (i == j)
	    {
		if (sum <= 0.0)
		    return TRUE;

		diag[i] = sqrt(sum);
	    }
	    else
		a[j*n+i] = sum / diag[i];
	}
    }
    
    return FALSE;
}

struct sc_entry
{
    double sc;
    int i, j;
} *sclist;

/* Sort descending */
int cmpfn(const void *a, const void *b)
{
    if (((struct sc_entry *)a)->sc == ((struct sc_entry *)b)->sc)
	return 0;

    if (((struct sc_entry *)a)->sc < ((struct sc_entry *)b)->sc)
	return 1;

    return -1;
}
    
/* Read Dirichlet mixture file */    
Regularizer* read_reg(const char* filename, const char* name, IdObject* reg_type){

    ifstream in(filename);
    if (!in.good())
    {	cerr << "Error: can't open " << filename
		<< " to read a regularizer\n"
		<< flush;
		return 0;
    }

    NamedClass* x= NamedClass::read_new(in);
    assert(x->is_a(Regularizer::classID()));
    assert(x->is_a(reg_type));
    Regularizer* r = dynamic_cast<Regularizer*>(x);
    r->set_name(name);
    return r;
}

/* Get 1D array index for 4D array, i.e. array[ww][ww][w][w] */
int index4D(int i, int j, int a, int b, int w, int ww){
	int row = (i*w)+a;
	int column = (j*w)+b;
	return((row*(w*ww))+column);
}

/* Get 1D array index for 2D array, i.e. array[ww][w] */
int index2D(int i, int j, int w){
	return(i*w + j);
}

int main(int argc, char **argv)
{
    int a, b, i, j, k, seqlen, nids, s, nseqs, ncon, opt, ndim, filtflg=0, approxflg=0, initflg=0, debugflg=0, diagpenflg=1, apcflg=1, maxit=10000, npair, nnzero, niter, jerr, shrinkflg=1, rawscflg = 1, pseudoc = 1, minseqsep = 5;
    unsigned int *wtcount, ccount[MAXSEQLEN];
    double thresh=1e-4, del, *pcsum, pcmean, pc, trialrho, rhodefault = -1.0, dscale = 5.0;
    double sum, score, wtsum, lambda, smean, fnzero, lastfnzero, rfact, r2, targfnzero = 0.0, scsum, scsumsq, mean, sd, zscore, ppv;    
    double *weight, idthresh = -1.0, maxgapf = 0.9;
    char buf[4096], seq[MAXSEQLEN], *blockfn = NULL, **aln;
    string mixture = "data/byst-4.5-0-3.9comp";
    FILE *ifp;

    while ((opt = getopt(argc, argv, "aflnpr:b:x:s:z:i:t:c:g:d:j:")) >= 0)
	switch (opt)
	{
	case 'a':
	    approxflg = 1;
	    break;
	case 'n':
	    shrinkflg = 0;
	    break;
	case 'p':
	    rawscflg = 0;
	    break;
	case 'f':
	    filtflg = 1;
	    break;
	case 'l':
	    apcflg = 0;
	    break;
	case 'r':
	    rhodefault = atof(optarg);
	    break;
	case 'd':
	    targfnzero = atof(optarg);
	    break;
	case 't':
	    thresh = atof(optarg);
	    break;
	case 'i':
	    idthresh = 1.0 - atof(optarg)/100.0;
	    break;
	case 'c':
	    pseudoc = atoi(optarg);
	    break;
	case 'j':
	    minseqsep = atoi(optarg);
	    break;
	case 'b':
	    blockfn = strdup(optarg);
	    break;
	case 'x':
	    std_alphabet = optarg;
	    break;
	case 'z':
	    mixture = optarg;
	    break;
	case 'g':
	    maxgapf = atof(optarg);
	    break;
	case 's':
	    dscale = atof(optarg);
	    break;
	case '?':
	    exit(-1);
	}

    if (optind >= argc)
		fail("Usage: psicov [options] alnfile\n\nOptions:\n-a\t: use approximate Lasso algorithm\n-n\t: don't pre-shrink the sample covariance matrix\n-f\t: filter low-scoring contacts\n-p\t: output PPV estimates rather than raw scores\n-l\t: don't apply APC to Lasso output\n-r nnn\t: set initial rho parameter\n-d nnn\t: set target precision matrix sparsity (default 0 = not specified)\n-t nnn\t: set Lasso convergence threshold (default 1e-4)\n-i nnn\t: select BLOSUM weighting with given identity threshold (default selects threshold automatically)\n-c nnn\t: set pseudocount value (default 1)\n-j nnn\t: set minimum sequence separation (default 5)\n-g nnn\t: maximum fraction of gaps (default 0.9)\n-b file\t: read rho parameter file\n-x file\t: location of alphabet file (default data/std.alphabet)\n-z file\t: location of mixture file (default data/byst-4.5-0-3.9comp)\n-s nnn\t: dirichlet scaling factor (default 5)\n");

    ifp = fopen(argv[optind], "r");
    if (!ifp)
	fail("Unable to open alignment file!");

    for (nseqs=0;; nseqs++)
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;

    aln = (char**)allocvec(nseqs, sizeof(char *));
    
    weight = (double*)allocvec(nseqs, sizeof(double));

    wtcount = (unsigned int*)allocvec(nseqs, sizeof(unsigned int));
    
    rewind(ifp);
    
    if (!fgets(seq, MAXSEQLEN, ifp))
	fail("Bad alignment file!");
    
    seqlen = strlen(seq)-1;

    if (!(aln[0] = (char*)malloc(seqlen)))
	fail("Out of memory!");

    for (j=0; j<seqlen; j++)
		aln[0][j] = aanum(seq[j]);

    for (i=1; i<nseqs; i++)
    {
	if (!fgets(seq, MAXSEQLEN, ifp))
	    break;
	
	if (seqlen != strlen(seq)-1)
	    fail("Length mismatch in alignment file!");
	
	if (!(aln[i] = (char*)malloc(seqlen)))
	    fail("Out of memory!");
	
	for (j=0; j<seqlen; j++)
	    aln[i][j] = aanum(seq[j]);

	}

	double **pa = allocmat(seqlen,21);
	double *pab = new double[seqlen*seqlen*21*21];

    /* Calculate pseudocounts from dirichlet mixtures*/
	
    Alphabet::set_default(Alph::ExtAA());    
    DirichletReg* Comp = dynamic_cast<DirichletReg*> (read_reg(mixture.c_str(), "dist.20", DirichletReg::classID()));
	double **dirichlet_probs = allocmat(seqlen,20);
	double dirichlet_scale = MIN(dscale,nseqs);

	for (j=0; j<seqlen; j++){

	    float probs[20] = {0};
	    float NumCounts[20] = {0};
	    for (i=0; i<nseqs; i++){
	    	NumCounts[aln[i][j]]++;
	    }
	    Comp->get_probs(NumCounts, probs);
	    //cout << "seq position " << j << ":\nint\tAA\tObs\tProb\n";
	    for(unsigned int p=0; p<20; p++){	
	    //	char lett = Alph::ExtAA().to_char(Alph::ExtAA().unindex(p));		
	    //	printf("%i\t%c\t%i\t%1.8f\n",p,lett,(int)NumCounts[p],probs[p]);
	    	dirichlet_probs[j][p] = probs[p];
		}
		//cout << "\n";
	}   	

    /* Calculate sequence weights */

    if (idthresh < 0.0)
    {
	double meanfracid = 0.0;
	
	for (i=0; i<nseqs; i++)
	    for (j=i+1; j<nseqs; j++)
	    {
		int nids;
		double fracid;
		
		for (nids=k=0; k<seqlen; k++)
		    if (aln[i][k] == aln[j][k])
			nids++;
		
		fracid = (double)nids / seqlen;
		
		meanfracid += fracid;
	    }
	
	meanfracid /= 0.5 * nseqs * (nseqs - 1.0);
	
	idthresh = 0.38 * 0.32 / meanfracid;
    }

    for (i=0; i<nseqs; i++)
	for (j=i+1; j<nseqs; j++)
	{
	    int nthresh = seqlen * (int)idthresh;

	    for (k=0; nthresh > 0 && k<seqlen; k++){
			if (aln[i][k] != aln[j][k])
		    	nthresh--;
	    }

	    if (nthresh > 0)
	    {
		wtcount[i]++;
		wtcount[j]++;
	    }
	}

    for (wtsum=i=0; i<nseqs; i++)
		wtsum += (weight[i] = 1.0 / (1 + wtcount[i]));
    
    //if (wtsum < MINEFSEQS)
	//	fail("Sorry - not enough sequences or sequence diversity to proceed!\nNeff (%f) < MINEFSEQS (%d)\nIf you want to force a calculation, try reducing MINEFSEQS at your own risk.\n", wtsum, MINEFSEQS);
  
    /* Calculate singlet frequencies with pseudocount */
    for (i=0; i<seqlen; i++)
    {
	for (a=0; a<21; a++){
	    if(a != 20){
	   		pa[i][a] = MAX(pseudoc,(dirichlet_probs[i][a]*dirichlet_scale));
		}else{
			pa[i][a] = pseudoc;
		}
	}

	for (k=0; k<nseqs; k++)
	{
	    a = aln[k][i];
	    if (a < 21)
		pa[i][a] += weight[k];
	}
	
	for (a=0; a<21; a++){
		if(a != 20){
			pa[i][a] /= MAX(pseudoc,(dirichlet_probs[i][a]*dirichlet_scale)) * 21.0 + wtsum;
		}else{	
	    	pa[i][a] /= pseudoc * 21.0 + wtsum;
		}
	}

    }

    /* Calculate pair frequencies with pseudocount */
    for (i=0; i<seqlen; i++)
    {
	for (j=i+1; j<seqlen; j++)
	{
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++){
			if(a < 20 && b < 20){
				pab[index4D(i,j,a,b,21,seqlen)] = MAX(pseudoc,(dirichlet_probs[i][a]*dirichlet_probs[j][b]*dirichlet_scale)) / 21.0;
			}else{	
		    	pab[index4D(i,j,a,b,21,seqlen)] = pseudoc / 21.0;
			}			
		}

	    for (k=0; k<nseqs; k++)
	    {
		a = aln[k][i];
		b = aln[k][j];
		if (a < 21 && b < 21)
		    pab[index4D(i,j,a,b,21,seqlen)] += weight[k];
	    }
	    
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		{
			if(a < 20 && b < 20){
				pab[index4D(i,j,a,b,21,seqlen)] /= MAX(pseudoc,(dirichlet_probs[i][a]*dirichlet_probs[j][b]*dirichlet_scale)) * 21.0 + wtsum;
			}else{	
			    pab[index4D(i,j,a,b,21,seqlen)] /= pseudoc * 21.0 + wtsum;
			}
		    pab[index4D(j,i,b,a,21,seqlen)] = pab[index4D(i,j,a,b,21,seqlen)];		
		}
	}
    }

    freemat(dirichlet_probs,seqlen,20);

    for (i=0; i<seqlen; i++)
	for (a=0; a<21; a++)
	    for (b=0; b<21; b++)
		pab[index4D(i,i,a,b,21,seqlen)] = (a == b) ? pa[i][a] : 0.0;

    double *cmat, *rho, *ww, *wwi, *tempmat;

    ndim = seqlen * 21;

    cmat = (double*)f_matrix_calloc(ndim, sizeof(double));
    tempmat = (double*)f_matrix_calloc(ndim, sizeof(double));

    /* Form the covariance matrix */
    for (i=0; i<seqlen; i++)
	for (j=0; j<seqlen; j++)
	    for (a=0; a<21; a++)
		for (b=0; b<21; b++)
		    if (i != j || a == b)
			cmat[(i*21+a) * ndim + j*21+b] = pab[index4D(i,j,a,b,21,seqlen)] - pa[i][a] * pa[j][b];

	delete [] pab;	

    /* Shrink sample covariance matrix towards shrinkage target F = Diag(1,1,1,...,1) * smean */

    if (shrinkflg)
    {
	for (smean=i=0; i<ndim; i++)
	    smean += cmat[i*ndim+i];
	
	smean /= (double)ndim;
	lambda = 0.1;

	for (;;)
	{
	    memcpy(tempmat, cmat, ndim*ndim*sizeof(double));
	    
	    /* Test if positive definite using Cholesky decomposition */
	    if (!test_cholesky(tempmat, ndim))
		break;
	    
	    for (i=0; i<seqlen; i++)
		for (j=0; j<seqlen; j++)
		    for (a=0; a<21; a++)
			for (b=0; b<21; b++)
			    if (i != j)
				cmat[(i*21+a)*ndim + j*21+b] *= 1.0 - lambda;
			    else if (a == b)
				cmat[(i*21+a)*ndim + j*21+b] = smean * lambda + (1.0 - lambda) * cmat[(i*21+a)*ndim + j*21+b];
	}
    }

    rho = (double*)f_matrix_calloc(ndim, sizeof(double));
    ww = (double*)f_matrix_calloc(ndim, sizeof(double));
    wwi = (double*)f_matrix_calloc(ndim, sizeof(double));

    lastfnzero=0.0;

    /* Guess at a reasonable starting rho value if undefined */
    if (rhodefault < 0.0)
	trialrho = MAX(0.001, 1.0 / wtsum);
    else
	trialrho = rhodefault;

    rfact = 0.0;

    for (;;)
    {
	if (trialrho <= 0.0 || trialrho >= 1.0)
	    fail("Sorry - failed to find suitable value for rho (0 < rho < 1)!");

	for (i=0; i<ndim; i++)
	    for (j=0; j<ndim; j++)
		rho[i*ndim + j] = trialrho;
	
	for (i=0; i<seqlen; i++)
	    for (j=0; j<seqlen; j++)
		for (a=0; a<21; a++)
		    for (b=0; b<21; b++)
			if ((a != b && i == j) || pa[i][20] > maxgapf || pa[j][20] > maxgapf)
			    rho[(i*21+a)*ndim + j*21+b] = 1e9;

	/* Mask out regions if block-out list provided */
	if (blockfn != NULL)
	{
	    ifp = fopen(blockfn, "r");
	    
	    for (;;)
	    {
		if (fscanf(ifp, "%d %d %lf", &i, &j, &score) != 3)
		    break;
		
		for (a=0; a<21; a++)
		    for (b=0; b<21; b++)
		    {
			rho[((i-1)*21+a)*ndim + (j-1)*21+b] = score;
			rho[((j-1)*21+b)*ndim + (i-1)*21+a] = score;
		    }
	    }
	    
	    fclose(ifp);
	}
    
	/* All matrices are symmetric so no need to transpose before/after calling Fortran code */

	glasso_(&ndim, cmat, rho, &approxflg, &initflg, &debugflg, &diagpenflg, &thresh, &maxit, ww, wwi, &niter, &del, &jerr);

	if (targfnzero <= 0.0)
	    break;

	for (npair=nnzero=i=0; i<ndim; i++)
	    for (j=i+1; j<ndim; j++,npair++)
		if (wwi[i*ndim+j] != 0.0)
		    nnzero++;

	fnzero = (double) nnzero / npair;

//  printf("rho=%f fnzero = %f\n", trialrho, fnzero);

	/* Stop iterating if we have achieved the target sparsity level */
	if (fabs(fnzero - targfnzero)/targfnzero < 0.01)
	    break;
	
	if (fnzero == 0.0)
	{
	    /* As we have guessed far too high, halve rho and try again */
	    trialrho *= 0.5;
	    continue;
	}
	
	if (lastfnzero > 0.0 && fnzero != lastfnzero)
	{
//	    printf("fnzero=%f lastfnzero=%f trialrho=%f oldtrialrho=%f\n", fnzero, lastfnzero, trialrho, trialrho/rfact);
	    
	    rfact = pow(rfact, log(targfnzero / fnzero) / log(fnzero / lastfnzero));

//	    printf("New rfact = %f\n", rfact);
	}

	lastfnzero = fnzero;

	/* Make a small trial step in the appropriate direction */

	if (rfact == 0.0)
	    rfact = (fnzero < targfnzero) ? 0.9 : 1.1;
	
	trialrho *= rfact;
    }

    free(rho);
    free(ww);
    freemat(pa,seqlen,21);		


    /* Calculate background corrected scores using average product correction */

    double **pcmat = allocmat(seqlen,seqlen);
    //pcmat = (double**)allocmat(seqlen, seqlen, sizeof(double));
    pcsum = (double*)allocvec(seqlen, sizeof(double));
    
    pcmean = 0.0;
    
    for (i=0; i<seqlen; i++)
	for (j=i+1; j<seqlen; j++)
	{	
	    for (pc=a=0; a<20; a++)
		for (b=0; b<20; b++)
		    pc += fabs(wwi[(i*21+a)*ndim + j*21+b]);

	    pcmat[i][j] = pcmat[j][i] = pc;
	    pcsum[i] += pc;
	    pcsum[j] += pc;

	    pcmean += pc;
	}

    pcmean /= seqlen * (seqlen - 1) * 0.5;

    /* Build final list of predicted contacts */

    sclist = (sc_entry*)allocvec(seqlen * (seqlen - 1) / 2, sizeof(struct sc_entry));

    for (scsum=scsumsq=ncon=i=0; i<seqlen; i++)
	for (j=i+minseqsep; j<seqlen; j++)
	    if (pcmat[i][j] > 0.0)
	    {
		/* Calculate APC score */
		if (apcflg)
		    sclist[ncon].sc = pcmat[i][j] - pcsum[i] * pcsum[j] / SQR(seqlen - 1.0) / pcmean;
		else
		    sclist[ncon].sc = pcmat[i][j];
		scsum += sclist[ncon].sc;
		scsumsq += SQR(sclist[ncon].sc);
		sclist[ncon].i = i;
		sclist[ncon++].j = j;
	    }

    qsort(sclist, ncon, sizeof(struct sc_entry), cmpfn);

    mean = scsum / ncon;
    sd = 1.25 * sqrt(scsumsq / ncon - SQR(mean)); /* Corrected for extreme-value bias */

    for (i=0; i<seqlen; i++)
	ccount[i] = 0;

    /* Print output in CASP RR format with optional PPV estimated from final Z-score */
    if (rawscflg)
	for (i=0; i<ncon; i++)
	    printf("%d %d 0 8 %f\n", sclist[i].i+1, sclist[i].j+1, sclist[i].sc);
    else
	for (i=0; i<ncon; i++)
	{
	    zscore = (sclist[i].sc - mean) / sd;
	    ppv = 0.904 / (1.0 + 16.61 * exp(-0.8105 * zscore));
	    if (ppv >= 0.5 || (!ccount[sclist[i].i] || !ccount[sclist[i].j]) || !filtflg)
	    {
		printf("%d %d 0 8 %f\n", sclist[i].i+1, sclist[i].j+1, ppv);
		ccount[sclist[i].i]++;
		ccount[sclist[i].j]++;
	    }
	}

    freemat(pcmat,seqlen,seqlen);

    return 0;
}
