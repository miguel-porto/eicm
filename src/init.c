#include "eicm.h"
#include <R_ext/Rdynload.h>

static const R_CallMethodDef callMethods[]  = {
	{"_simulate_community_probability", (DL_FUNC) &_simulate_community_probability, 7},
	{"_likelihood", (DL_FUNC) &_likelihood, 4},
	{"_likelihood_NAallowed", (DL_FUNC) &_likelihood_NAallowed, 4},
	{"_mathematical_likelihood_fast", (DL_FUNC) &_mathematical_likelihood_fast, 4},
	{"_likelihood_superfast", (DL_FUNC) &_likelihood_superfast, 4},
	{"_likelihood_superfast_NAallowed", (DL_FUNC) &_likelihood_superfast_NAallowed, 4},
	{"_makeModelMatricesFromPars", (DL_FUNC) &_makeModelMatricesFromPars, 5},
	{"_makeParsFromModelMatrices", (DL_FUNC) &_makeParsFromModelMatrices, 2},
	{"_getNumberOfParameters", (DL_FUNC) &_getNumberOfParameters, 3},
	{"_getMostSimilarModel", (DL_FUNC) &_getMostSimilarModel, 2},
	{"_isCyclic", (DL_FUNC) &_isCyclic, 1}, 
	{NULL, NULL, 0}
};

void R_init_eicm(DllInfo *info) {

	R_registerRoutines(info, NULL, callMethods, NULL, NULL);
	R_useDynamicSymbols(info, FALSE);
	R_forceSymbols(info, TRUE);
	
	if(logInvLogitTab == NULL) {
//		Rprintf("Creating logit lookup tables for fast estimation");
//		R_FlushConsole();
		logInvLogitTab = malloc(LOOKUPTABLESIZE * sizeof(double));
		if(logInvLogitTab == NULL)
			Rf_error("Could not allocate memory. Please contact the maintainer.");
			
		logSymmInvLogitTab = malloc(LOOKUPTABLESIZE * sizeof(double));
		if(logSymmInvLogitTab == NULL)
			Rf_error("Could not allocate memory. Please contact the maintainer.");
			
//		Rprintf("...");
//		R_FlushConsole();
		
		for(int i=0; i<LOOKUPTABLESIZE; i++) {
			logInvLogitTab[i] = log(1 / (1 + exp(-(double) (i - MIDDLELOOKUPTABLE) / DIVISOR)));
			logSymmInvLogitTab[i] = log(1- 1 / (1 + exp(-(double) (i - MIDDLELOOKUPTABLE) / DIVISOR)));
		}
/*		Rprintf("done!\nTail:\n");
		for(int i=0; i<10; i++)
			Rprintf("%f ", logInvLogitTab[i]);
		Rprintf("\nMiddle:\n");
		for(int i=MIDDLELOOKUPTABLE - 5; i<MIDDLELOOKUPTABLE + 5; i++)
			Rprintf("%f ", logInvLogitTab[i]);
		Rprintf("\n");
*/
	}
}

void R_unload_eicm(DllInfo *info) {
	Rprintf("Releasing logit lookup tables\n");
	if(logInvLogitTab != NULL) {
		free(logInvLogitTab);
		logInvLogitTab = NULL;
	}
	if(logSymmInvLogitTab != NULL) {
		free(logSymmInvLogitTab);
		logSymmInvLogitTab = NULL;
	}
}
