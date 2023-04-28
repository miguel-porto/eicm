#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <stdbool.h>
#include "todos.h"

#define LOW_MEMORY
#define MAXRND 32768

#define NELEMS(x)  (sizeof(x) / sizeof((x)[0]))

#ifdef LOW_MEMORY
	#define DIVISOR 100000.0
	#define LOOKUPTABLESIZE 4000000		// must be even. This value must be large, otherwise the gradient-based minimizer will get stuck
	#define MIDDLELOOKUPTABLE 2000000.0	// must be LOOKUPTABLESIZE / 2
	#define LOWERBOUND -20	// this must be -MIDDLELOOKUPTABLE / DIVISOR
	#define UPPERBOUND 20	// this must be MIDDLELOOKUPTABLE / DIVISOR
#endif
#ifndef LOW_MEMORY
	#define DIVISOR 1000000.0
	#define LOOKUPTABLESIZE 40000000		// must be even. This value must be large, otherwise the gradient-based minimizer will get stuck
	#define MIDDLELOOKUPTABLE 20000000.0	// must be LOOKUPTABLESIZE / 2
	#define LOWERBOUND -20	// this must be -MIDDLELOOKUPTABLE / DIVISOR. These are the LP bounds, outside of which, prob is set to 0 or 1.
	#define UPPERBOUND 20	// this must be MIDDLELOOKUPTABLE / DIVISOR
#endif

extern double *logInvLogitTab;
extern double *logSymmInvLogitTab;

SEXP _simulate_community_probability(SEXP _nrepetitions, SEXP _env, SEXP _sp, SEXP _envcoefs, SEXP _spcoefs
	, SEXP _cycles, SEXP _randomseed);
SEXP _likelihood(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed);
SEXP _likelihood_NAallowed(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed);
SEXP _mathematical_likelihood_fast(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed);
SEXP _likelihood_superfast(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed);
SEXP _likelihood_superfast_NAallowed(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed);
SEXP _makeModelMatricesFromPars(SEXP _pars, SEXP _nspecies, SEXP _nenv, SEXP _mask, SEXP _offset);
SEXP _makeParsFromModelMatrices(SEXP _matrices, SEXP _mask);
SEXP _getNumberOfParameters(SEXP _nspecies, SEXP _nenv, SEXP _mask);
SEXP _getMostSimilarModel(SEXP _popToEval, SEXP _cachedModelList);
SEXP _isCyclic(SEXP _spcoefs);

void createInverseLinkFunctionTable(unsigned int *table);
void createInverseLinkFunctionTableProb(float *table);
//void matProd(double *m1, double *m2, double *out, int m1rows, int m2rows, int m1cols);
void matProd(const double *m1, const double *m2, double *out, const int m1rows, const int m2rows, const int m1cols);
void matProdExclSpecies(const double *m1, const double *m2, double *out, const int m1rows, const int m2rows, const int m1cols, const bool *exclSpecies);
void matProdInt(int *m1, double *m2, double *out, int m1rows, int m2rows, int m1cols
	, int *whichm1cols, int lenwhichm1cols, int *whichm2rows, int lenwhichm2rows);
void matProdShort(const short *m1, const double *m2, double *out, const int m1rows, const int m2rows, const int m1cols
	, const int *whichm1cols, const int lenwhichm1cols, const int *whichm2rows, const int lenwhichm2rows);
void matProdShortExclSamp(short *m1, double *m2, double *out, int m1rows, int m2rows, int m1cols
	, int *whichm1cols, int lenwhichm1cols, int *whichm2rows, int lenwhichm2rows, short *exclude, short maxErrors);
void matProdShortExclSamp2(short *prM, double *coM, double *out, int prMrows, int coMrows, int prMcols
	, bool *filledSpecies, int *whichcoMrows, int lenwhichcoMrows, bool *excludeSamples);
short *computeDependencyMatrix(SEXP _spcoefs);

static unsigned int g_seed;

//Used to seed the generator.
static inline void fast_srand( int seed ) {
	g_seed = seed;
}

//fastrand routine returns one integer, similar output value range as C lib.
static inline int fastrand(void) {
	g_seed = (214013*g_seed+2531011);
	return (g_seed>>16)&0x7FFF;
}

inline double invLogit(const double n) {
	return(1 / (1 + exp(-n)));
}

