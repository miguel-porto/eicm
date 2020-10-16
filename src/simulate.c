#include "eicm.h"
extern inline double invLogit(const double n);

SEXP _simulate_community_probability(SEXP _nrepetitions, SEXP _env, SEXP _sp, SEXP _envcoefs
	, SEXP _spcoefs, SEXP _cycles, SEXP _randomseed) {
	int nrepetitions = 					INTEGER_POINTER(_nrepetitions)[0];
	int randomseed = 					INTEGER_POINTER(_randomseed)[0];
	int nspecies = 						Rf_nrows(_envcoefs);
	int nsamples = 						Rf_nrows(_env);
	int nenv = 							Rf_ncols(_envcoefs);
	int ncycles = 						LENGTH(_cycles);
	int sizeout = 						nspecies * nsamples;
	
	double *env = 						NUMERIC_POINTER(_env);
	int *knownpresences = 				isNull(_sp) ? NULL : INTEGER_POINTER(_sp);
	double *envcoefs = 					NUMERIC_POINTER(_envcoefs);
	double *spcoefs = 					NUMERIC_POINTER(_spcoefs);
	
	double *spenv = 					malloc(sizeof(double) * sizeout);
	double *spenvori = 					malloc(sizeof(double) * sizeout);
	short *presences = 					malloc(sizeof(short) * sizeout);
	bool *excludeSpecies =				isNull(_sp) ? NULL : malloc(sizeof(bool) * nspecies);
	SEXP thiscycle						;
	double *pprobabilities				;
	int *pthiscycle;
//	unsigned int invLogitTab[LOOKUPTABLESIZE];
		
	GetRNGstate();
	fast_srand(randomseed);
	
//	createInverseLinkFunctionTable(invLogitTab);
	
	SEXP probabilities;
	PROTECT(probabilities = Rf_allocMatrix(REALSXP, nsamples, nspecies));
	pprobabilities = NUMERIC_POINTER(probabilities);
	
	memset(pprobabilities, 0, sizeof(double) * sizeout);

	int *filledSpp = calloc(nspecies, sizeof(int));
	int currentFilledSp;

	if(Rf_ncols(_env) != nenv) Rf_error("Number of predictors in environmental matrix is not the same as the number of columns in the coefficient matrix. Did you forget the intercept column?");

	if(knownpresences == NULL) {
		// Compute environmental contribution
		// spenvori (samples x species) is the LP before interactions
		matProd(env, envcoefs, spenvori, nsamples, nspecies, nenv);
	} else {
		// check the first element of the known species matrix. If it's not NA, then this species is to be
		// excluded from predictions, because it is given.
		// TODO we can allow NAs and givens in the same species? Sure we can.
		for(int i=0; i<nspecies; i++) {
			excludeSpecies[i] = (knownpresences[i * nsamples] != NA_INTEGER);
			if(excludeSpecies[i]) Rprintf("Excluded species %d\n", i + 1);
		}
		matProdExclSpecies(env, envcoefs, spenvori, nsamples, nspecies, nenv, excludeSpecies);
		// spenvori is undefined for excluded species
	}
	int r, i;

	for(r=0; r<nrepetitions; r++) {
		memset(filledSpp, 0, nspecies);
		currentFilledSp = 0;
		
		if(knownpresences != NULL) {
			// so we have already filled species, get the list of them
			for(i=0; i<nspecies; i++) {
				if(excludeSpecies[i]) {
					filledSpp[currentFilledSp] = i + 1;	// because filledSpp indices are 1-based
					currentFilledSp ++;
					// we fill in a priori the presences of the known species
					int ptr = i * nsamples;
					for(int sa = 0; sa < nsamples; sa++)
						presences[sa + ptr] = knownpresences[sa + ptr];
				}
			}
		}
		
		// restore original environmental matrix
		memcpy(spenv, spenvori, sizeof(double) * sizeout);

		// We fill in the presences by dependency cycles.
		// We start with the independent species.
		for(i = 0; i < ncycles; i++) {
			thiscycle = VECTOR_ELT(_cycles, i);
			int nsp = LENGTH(thiscycle);
			pthiscycle = INTEGER_POINTER(thiscycle);

/*Rprintf("\nCycle %d: ", i);
for(int j=0; j<LENGTH(thiscycle); j++) Rprintf("%d ", pthiscycle[j]);
Rprintf("\n");*/

			if(i > 0) {	// i==0 are those that depend only on the environment
				// sum the environmental contribution with the interaction contribution
				// TODO if some species presences are given
				matProdShort(presences, spcoefs, spenv, nsamples, nspecies, nspecies
					, filledSpp, currentFilledSp, pthiscycle, nsp);
			}
/*
	for(int ti=0; ti<nsamples; ti++) {
		for(int tj=0; tj<nspecies; tj++) {
			Rprintf("%.03f ", spenv[ti + tj*nsamples]);
		}
		Rprintf("\n");
	}
*/
			// for each species in this cycle
			for(int sp = 0, ptr2 = 0; sp < nsp; sp++) {
				// this points to the column of the species.
				// note that we must fetch the species index from the cycles
				// Species indexes in pthiscycle are 1-based
				ptr2 = (pthiscycle[sp] - 1) * nsamples;

				if(knownpresences != NULL && excludeSpecies[pthiscycle[sp] - 1])
					continue;

				// for each sample
				for(int sa = 0; sa < nsamples; sa++) {
					// apply the inverse logit function
					// this is the human-readable slow way
					const double prob = invLogit(spenv[sa + ptr2]);
					presences[sa + ptr2] = (fastrand() <= (prob * MAXRND)) ? 1 : 0;

					// this is the optimized way
/*					double LP = spenv[sa + ptr2];
					if(LP >= UPPERBOUND) {
						presences[sa + ptr2] = 1;
					} else if(LP <= LOWERBOUND) {
						presences[sa + ptr2] = 0;
					} else
						presences[sa + ptr2] = (fastrand() <= invLogitTab[(unsigned int) (LP * DIVISOR + MIDDLELOOKUPTABLE)]) ? 1 : 0;
*/
					//ppresences[sa + sp * nsamples + r * nsamples * nspecies] = rbinom(1, prob);
					// make one realization of the probability
					//presences[sa + ptr2] = (short) rbinom(1, prob);
				}
				if(knownpresences != NULL) {
					// incrementally add to filledspecies
					filledSpp[currentFilledSp] = pthiscycle[sp];
					currentFilledSp ++;
				}
			}
			if(knownpresences == NULL) {
				memcpy(&filledSpp[currentFilledSp], pthiscycle, LENGTH(thiscycle) * sizeof(int));
				currentFilledSp += LENGTH(thiscycle);
			}
	/*		Rprintf("\nFilled spp: ");
			for(int i=0; i<nspecies; i++)
				Rprintf("%d ", filledSpp[i]);*/
		}
		
		// make the presence count in all samples, all species
		for(int i=0; i < sizeout; i++)
			pprobabilities[i] += presences[i];
	}
	
	for(int i=0; i < sizeout; i++)
		pprobabilities[i] /= nrepetitions;

	free(spenv);
	free(spenvori);
	free(filledSpp);
	free(presences);
	if(excludeSpecies != NULL) free(excludeSpecies);
	PutRNGstate();
	UNPROTECT(1);
	return probabilities;
}

