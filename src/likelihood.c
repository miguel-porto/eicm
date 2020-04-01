#include "eicm.h"
#include<stdio.h>
double *logInvLogitTab;
double *logSymmInvLogitTab;

/**
 * Computes the true likelihood of the model, given a presence/absence matrix.
 * The true likelihood is the probability of observing, in each sample, the exact community provided,
 * i.e. the exact presence-absence pattern across all species in each sample.
 * This function is optimized for speed, so it is less readable than the prediction function.
*/
SEXP _likelihood(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed) {
	int nspecies = 				Rf_nrows(_envcoefs);
	int nsamples = 				Rf_nrows(_env);
	int nenv = 					Rf_ncols(_envcoefs);
	int sizeX = 				nspecies * nsamples;

	double *env = 				NUMERIC_POINTER(_env);
	double *envcoefs = 			NUMERIC_POINTER(_envcoefs);
	double *spcoefs = 			NUMERIC_POINTER(_spcoefs);
	int *observed = 			INTEGER_POINTER(_observed);	

	double *envLP = 			malloc(sizeof(double) * sizeX);
	if(envLP == NULL) return R_NilValue;
	short *dependency = 		computeDependencyMatrix(_spcoefs);
	if(dependency == NULL) return R_NilValue;
	
	SEXP out;
	
	out = PROTECT(				NEW_NUMERIC(nsamples));
	double *pOut = 				NUMERIC_POINTER(out);
	

	// Compute the environmental Linear Predictor in each sample, for all species
	matProd(env, envcoefs, envLP, nsamples, nspecies, nenv);
//	memset(envLP, 0, sizeof(double) * sizeX);

	// for each sample	
	for(int sa=0; sa<nsamples; sa++) {
		double prob = 0;
		// for each species
		for(int spi=0, ptr2=sa; spi<nspecies; spi++, ptr2 += nsamples) {
//			ptr = spi * nsamples + sa;
			
			double LP = envLP[ptr2];	// environmental contribution
			
			// For this species, iterate through all dependencies (species upon which this one depends)
			for(short depi=0; depi<nspecies; depi++) {
				short depsp = dependency[spi + depi * nspecies];	// this species on which it depends
				if(depsp == -1) break;		// no more dependencies
				if(observed[depsp * nsamples + sa] == 1) {	// is this other species present? if yes, sum coefficient to LP
					LP += spcoefs[spi + depsp * nspecies];		// the coefficient of this species depending on each other
//					hash += 1 << depi;
				}
			}
//			Rprintf("LP: %f + %f = %f\n", LP, condP, LP + condP);
			
			double tmp;// = 1 / (1 + exp(-LP));
			if(LP >= UPPERBOUND)
				tmp = 0.999999999;
			else if(LP <= LOWERBOUND)
				tmp = 0.000000001;
			else
				tmp = 1 / (1 + exp(-LP));

			prob += observed[ptr2] == 1 ?
				log(tmp)
				: log(1 - tmp);

//			log(invlink(p1[3])) + log(1-invlink(p2_p1[3])) + log(1-invlink(p3_p1a2[3]))

		}
		pOut[sa] = prob;
	}

	free(dependency);
	free(envLP);
	UNPROTECT(1);
	return out;
}

// This function allows for sparse occurrence data, i.e. there can be scatterd NAs
// anywhere (no need for entire species be NA)
SEXP _likelihood_NAallowed(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed) {
	int nspecies = 				Rf_nrows(_envcoefs);
	int nsamples = 				Rf_nrows(_env);
	int nenv = 					Rf_ncols(_envcoefs);
	int sizeX = 				nspecies * nsamples;

	double *env = 				NUMERIC_POINTER(_env);
	double *envcoefs = 			NUMERIC_POINTER(_envcoefs);
	double *spcoefs = 			NUMERIC_POINTER(_spcoefs);
	int *observed = 			INTEGER_POINTER(_observed);	

	double *envLP = 			malloc(sizeof(double) * sizeX);
	if(envLP == NULL) return R_NilValue;
	short *dependency = 		computeDependencyMatrix(_spcoefs);
	if(dependency == NULL) return R_NilValue;
	
	SEXP out;
	
	out = PROTECT(				NEW_NUMERIC(nsamples));
	double *pOut = 				NUMERIC_POINTER(out);
	

	// Compute the environmental Linear Predictor in each sample, for all species
	matProd(env, envcoefs, envLP, nsamples, nspecies, nenv);
//	memset(envLP, 0, sizeof(double) * sizeX);

	// for each sample	
	for(int sa=0; sa<nsamples; sa++) {
		double prob = 0;
		// for each species
		for(int spi=0, ptr2=sa; spi<nspecies; spi++, ptr2 += nsamples) {
			if(observed[ptr2] == NA_INTEGER) {	// no data of this species in this sample
				continue;	// ignore, probability won't take it into account
				// note that NA is not the same as 0. If the species is absent, the probability is changed
				// (see end of loop)
			}
			double LP = envLP[ptr2];	// environmental contribution
			
			// For this species, iterate through all dependencies (species upon which this one depends)
			for(short depi=0; depi<nspecies; depi++) {
				short depsp = dependency[spi + depi * nspecies];	// this species on which it depends
				if(depsp == -1) break;		// no more dependencies
				if(observed[depsp * nsamples + sa] == 1) {	// is this other species present? if yes, sum coefficient to LP
					LP += spcoefs[spi + depsp * nspecies];		// the coefficient of this species depending on each other
//					hash += 1 << depi;
				}
			}
//			Rprintf("LP: %f + %f = %f\n", LP, condP, LP + condP);
			
			double tmp;// = 1 / (1 + exp(-LP));
			if(LP >= UPPERBOUND)
				tmp = 0.999999999;
			else if(LP <= LOWERBOUND)
				tmp = 0.000000001;
			else
				tmp = 1 / (1 + exp(-LP));

			prob += observed[ptr2] == 1 ?
				log(tmp)
				: log(1 - tmp);

//			log(invlink(p1[3])) + log(1-invlink(p2_p1[3])) + log(1-invlink(p3_p1a2[3]))

		}
		pOut[sa] = prob;
	}

	free(dependency);
	free(envLP);
	UNPROTECT(1);
	return out;
}

// TODO further optimize for speed. Consider using only fixed-point math?
SEXP _mathematical_likelihood_fast(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed) {
	error("Deprecated");
	int nspecies = 				Rf_nrows(_envcoefs);
	int nsamples = 				Rf_nrows(_env);
	int nenv = 					Rf_ncols(_envcoefs);
	int sizeX = 				nspecies * nsamples;

	double *env = 				NUMERIC_POINTER(_env);
	double *envcoefs = 			NUMERIC_POINTER(_envcoefs);
	double *spcoefs = 			NUMERIC_POINTER(_spcoefs);
	int *observed = 			INTEGER_POINTER(_observed);	

	double *envLP = 			malloc(sizeof(double) * sizeX);
	if(envLP == NULL) return R_NilValue;
	short *dependency = 		computeDependencyMatrix(_spcoefs);
	if(dependency == NULL) return R_NilValue;
	
	SEXP out;
	
	out = PROTECT(				NEW_NUMERIC(nsamples));
	double *pOut = 				NUMERIC_POINTER(out);
	

	// TODO this should be on load hook
/*	if(logInvLogitTab[0] == 0) {
		for(int i=0; i<LOOKUPTABLESIZE; i++) {
			logInvLogitTab[i] = log(1 / (1 + exp(-(double) (i - MIDDLELOOKUPTABLE) / DIVISOR)));
			logSymmInvLogitTab[i] = log(1- 1 / (1 + exp(-(double) (i - MIDDLELOOKUPTABLE) / DIVISOR)));
		}
	}
*/

//	float invLogitTab[LOOKUPTABLESIZE];
/*	createInverseLinkFunctionTableProb(invLogitTab);*/

	// Compute the environmental Linear Predictor in each sample, for all species
	matProd(env, envcoefs, envLP, nsamples, nspecies, nenv);

	for(int sa=0; sa<nsamples; sa++) {
		double prob = 0;
		for(int spi=0, ptr2=sa; spi<nspecies; spi++, ptr2 += nsamples) {
//			ptr = spi * nsamples + sa;
			double LP = envLP[ptr2];	// environmental contribution
			
			// For this species, iterate through all dependencies (species upon which this one depends)
			for(short depi=0; depi<nspecies; depi++) {
				short depsp = dependency[spi + depi * nspecies];	// this species on which it depends
				if(depsp == -1) break;		// no more dependencies
				// we want the hash of the observed presence-absence pattern, so we can fetch the conditional probability
				if(observed[depsp * nsamples + sa] == 1) {	// is this other species present? if yes, sum coefficient to LP
					LP += spcoefs[spi + depsp * nspecies];		// the coefficient of this species depending on each other
//					hash += 1 << depi;
				}
			}
//			Rprintf("LP: %f + %f = %f\n", LP, condP, LP + condP);
			
//			if(hashtable[hash] == -1) {
				// compute conditional probability
//				cond = 
//				spcoefs[spi + depsp * nspecies]
//			}

			if(observed[ptr2] == 1) {
				if(LP >= UPPERBOUND)
					prob += logInvLogitTab[LOOKUPTABLESIZE - 1];
				else if(LP <= LOWERBOUND)
					prob += logInvLogitTab[0];
				else
					prob += logInvLogitTab[(unsigned int) (LP * DIVISOR + MIDDLELOOKUPTABLE)];
			} else {
				if(LP >= UPPERBOUND)
					prob += logSymmInvLogitTab[LOOKUPTABLESIZE - 1];
				else if(LP <= LOWERBOUND)
					prob += logSymmInvLogitTab[0];
				else
					prob += logSymmInvLogitTab[(unsigned int) (LP * DIVISOR + MIDDLELOOKUPTABLE)];
			}

		}
		pOut[sa] = prob;
	}

	free(dependency);
	free(envLP);
	UNPROTECT(1);
	return out;
}

SEXP _likelihood_superfast(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed) {
	const int nspecies =		Rf_nrows(_envcoefs);
	const int nsamples = 		Rf_nrows(_env);
	const int nenv = 			Rf_ncols(_envcoefs);
	const int sizeX =			nspecies * nsamples;

	const double *env =			NUMERIC_POINTER(_env);
	const double *envcoefs = 	NUMERIC_POINTER(_envcoefs);
	const double *spcoefs = 	NUMERIC_POINTER(_spcoefs);
	const int *observed = 		INTEGER_POINTER(_observed);	

	double *envLP =		malloc(sizeof(double) * sizeX);	
	if(envLP == NULL) return R_NilValue;
	short *dependency = 		computeDependencyMatrix(_spcoefs);
	if(dependency == NULL) return R_NilValue;
	
	SEXP out = 					PROTECT(NEW_NUMERIC(1));
	double *pOut = 				NUMERIC_POINTER(out);
	double prob = 				0;
	

//	FILE *pFile2 = fopen("myfile2.txt", "a");

	// Compute the environmental Linear Predictor in each sample, for all species
	matProd(env, envcoefs, envLP, nsamples, nspecies, nenv);

	for(int spi=0; spi<nspecies; spi++) {
		const register int ptr2 = spi * nsamples;
		for(short depi=0; depi<nspecies; depi++) {
			const short depsp = dependency[spi + depi * nspecies];	// this species on which it depends
			if(depsp == -1) break;		// no more dependencies
			const register int ptr1 = depsp * nsamples;
			const register double spc = spcoefs[depsp * nspecies + spi];
	
			for(int sa=0; sa<nsamples; sa++) {
				if(observed[ptr1 + sa] == 1) {	// is this other species present? if yes, sum coefficient to LP
					envLP[sa + ptr2] += spc;		// the coefficient of this species depending on each other
				}
			}
		}
	}

/*	for(int i=0; i<sizeX; i++) {
		if(observed[i] == 1) 
			prob *= 0.01 * envLP[i] + 0.5;
		else
			prob *= ((-0.01) * envLP[i]) + 0.5;
	}
	*pOut = -prob;
*/

	for(int i=0; i<sizeX; i++) {
		if(observed[i] == 1) {
			if(envLP[i] >= UPPERBOUND)
				prob += logInvLogitTab[LOOKUPTABLESIZE - 1];
			else if(envLP[i] <= LOWERBOUND)
				prob += logInvLogitTab[0];
			else {
//				fprintf(pFile2, "%d\n", tmp);
				prob += logInvLogitTab[(unsigned int) (envLP[i] * DIVISOR + MIDDLELOOKUPTABLE)];
			}
		} else {
			if(envLP[i] >= UPPERBOUND)
				prob += logSymmInvLogitTab[LOOKUPTABLESIZE - 1];
			else if(envLP[i] <= LOWERBOUND)
				prob += logSymmInvLogitTab[0];
			else
				prob += logSymmInvLogitTab[(unsigned int) (envLP[i] * DIVISOR + MIDDLELOOKUPTABLE)];
		}
	}
	*pOut = prob;

//	fclose(pFile2);
	
	free(dependency);
	free(envLP);
	UNPROTECT(1);
	return out;
}

SEXP _likelihood_superfast_NAallowed(SEXP _env, SEXP _envcoefs, SEXP _spcoefs, SEXP _observed) {
	const int nspecies =		Rf_nrows(_envcoefs);
	const int nsamples = 		Rf_nrows(_env);
	const int nenv = 			Rf_ncols(_envcoefs);
	const int sizeX =			nspecies * nsamples;

	const double *env =			NUMERIC_POINTER(_env);
	const double *envcoefs = 	NUMERIC_POINTER(_envcoefs);
	const double *spcoefs = 	NUMERIC_POINTER(_spcoefs);
	const int *observed = 		INTEGER_POINTER(_observed);	

	double *envLP =		malloc(sizeof(double) * sizeX);
	if(envLP == NULL) return R_NilValue;
	short *dependency = 		computeDependencyMatrix(_spcoefs);
	if(dependency == NULL) return R_NilValue;
	
	SEXP out = 					PROTECT(NEW_NUMERIC(1));
	double *pOut = 				NUMERIC_POINTER(out);
	double prob = 				0;
	

//	FILE *pFile2 = fopen("myfile2.txt", "a");

	// Compute the environmental Linear Predictor in each sample, for all species
	matProd(env, envcoefs, envLP, nsamples, nspecies, nenv);

	for(int spi=0; spi<nspecies; spi++) {
		const register int ptr2 = spi * nsamples;
		if(observed[ptr2] == NA_INTEGER) continue;
		for(short depi=0; depi<nspecies; depi++) {
			const short depsp = dependency[spi + depi * nspecies];	// this species on which it depends
			if(depsp == -1) break;		// no more dependencies
			const register int ptr1 = depsp * nsamples;
			const register double spc = spcoefs[depsp * nspecies + spi];
	
			for(int sa=0; sa<nsamples; sa++) {
				if(observed[ptr1 + sa] == 1) {	// is this other species present? if yes, sum coefficient to LP
					envLP[sa + ptr2] += spc;		// the coefficient of this species depending on each other
				}
			}
		}
	}

/*	for(int i=0; i<sizeX; i++) {
		if(observed[i] == 1) 
			prob *= 0.01 * envLP[i] + 0.5;
		else
			prob *= ((-0.01) * envLP[i]) + 0.5;
	}
	*pOut = -prob;
*/

	for(int i=0; i<sizeX; i++) {
		if(observed[i] == NA_INTEGER)
			continue;
		else if(observed[i] == 1) {
			if(envLP[i] >= UPPERBOUND)
				prob += logInvLogitTab[LOOKUPTABLESIZE - 1];
			else if(envLP[i] <= LOWERBOUND)
				prob += logInvLogitTab[0];
			else {
//				fprintf(pFile2, "%d\n", tmp);
				prob += logInvLogitTab[(unsigned int) (envLP[i] * DIVISOR + MIDDLELOOKUPTABLE)];
			}
		} else {
			if(envLP[i] >= UPPERBOUND)
				prob += logSymmInvLogitTab[LOOKUPTABLESIZE - 1];
			else if(envLP[i] <= LOWERBOUND)
				prob += logSymmInvLogitTab[0];
			else
				prob += logSymmInvLogitTab[(unsigned int) (envLP[i] * DIVISOR + MIDDLELOOKUPTABLE)];
		}
	}
	*pOut = prob;

//	fclose(pFile2);
	
	free(dependency);
	free(envLP);
	UNPROTECT(1);
	return out;
}

short *computeDependencyMatrix(SEXP _spcoefs) {
	double *spcoefs = 			NUMERIC_POINTER(_spcoefs);
	int nspecies = 				Rf_nrows(_spcoefs);
	short *dependency = 		malloc(nspecies * nspecies * sizeof(short));
//	short *ndeps =		 		calloc(nspecies, sizeof(short));
	short ndeps[1000];
	int i, j, base;
	
	if(dependency == NULL) return NULL;
	
	memset(ndeps, 0, nspecies * sizeof(short));
	memset(dependency, 0xff, nspecies * nspecies * sizeof(short));
	
	for(i=0, base=0; i<nspecies; i++, base += nspecies) {
		for(j=0; j<nspecies; j++) {
			if(i == j || fabs(spcoefs[j + base]) < 0.001) continue;
			dependency[j + ndeps[j] * nspecies] = i;
			ndeps[j] ++;
		}
	}
/*	
	for(j=0; j<nspecies; j++) {
		Rprintf("Ndep %d: ", ndeps[j]);
//		for(i=0; i<ndeps[j]; i++)
		for(i=0; i<nspecies; i++)
			Rprintf("%2d ", dependency[j + i*nspecies]);
		Rprintf("\n");
	}
*/
	//free(ndeps);
	return dependency;
}

SEXP _isCyclic(SEXP _spcoefs) {
	double *spcoefs = 			NUMERIC_POINTER(_spcoefs);
	int nspecies = 				Rf_nrows(_spcoefs);
	bool *leafnodes = 			malloc(nspecies * sizeof(bool));
	bool *checkednodes =		calloc(nspecies, sizeof(bool));
	int i, j, base;
	int count, total = 			nspecies;

	for(;;) {	
		memset(leafnodes, true, nspecies * sizeof(bool));
		count = 0;
		for(i=0, base=0; i<nspecies; i++, base += nspecies) {
			if(checkednodes[i]) continue;
			for(j=0; j<nspecies; j++) {
				if(checkednodes[j]) continue;
				if(spcoefs[j + base] != 0) {
					leafnodes[i] = false;
					count ++;
					break;
				}
			}
		}
		if(count == total) {	// no leaf nodes!
			free(leafnodes);
			free(checkednodes);
			return ScalarLogical(true);
		} else if(count == 0) {	// all leaf nodes
			free(leafnodes);
			free(checkednodes);
			return ScalarLogical(false);
		}
	
		for(j=0, total=0; j<nspecies; j++) {
			if(leafnodes[j])
				checkednodes[j] = true;
			else
				total ++;
		}
	}
	free(leafnodes);
	free(checkednodes);
}

