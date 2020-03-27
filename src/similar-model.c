#include "eicm.h"

/**
* Returns the index of the model (within cachedModelList) that is most similar to each of the given model candidates
* popToEval: numeric matrix with 0/1, each row is a candidate model
* cachedModelList: a list whose names are the bit strings
*/
SEXP _getMostSimilarModel(SEXP _popToEval, SEXP _cachedModelList) {
	double *popToEval = 					REAL(_popToEval);
	int nModelsToEval = 					Rf_nrows(_popToEval);
	int nCachedModels = 					LENGTH(_cachedModelList);
	SEXP cachedModelBits = 					GET_NAMES(_cachedModelList);
	int nbits = 							Rf_ncols(_popToEval);
	int i, j;
	const char *name;
	SEXP out;
	int *pout;
	
	if(isNull(cachedModelBits))
		Rf_error("cachedModelList must have names.");

	int *dissimil = 						malloc(nCachedModels * sizeof(int));
	
	PROTECT(out = NEW_INTEGER(nModelsToEval));
	pout = INTEGER(out);
	
	for(j=0; j<nModelsToEval; j++) {
//		Rprintf("Model %d\n", j);
		int minV = 100000, minI = -1;
		for(i=0; i<nCachedModels; i++) {
			name = STRING_VALUE(STRING_ELT(cachedModelBits, i));
			int tmpdiss = 0;
			for(int k=0; k<nbits; k++)
				tmpdiss += (name[k] == '0') != (popToEval[k * nModelsToEval + j] == 0);

			dissimil[i] = tmpdiss;
			if(tmpdiss < minV) {
				minV = tmpdiss;
				minI = i;
			}
//			Rprintf("%s %d\n", name, tmpdiss);
		}
		
		if(minI == -1) {
			pout[j] = NA_INTEGER;
		} else {
			if(IS_LOGICAL(VECTOR_ELT(_cachedModelList, minI))) {
				// if the minimum is an impossible, must do further search.
				do {
					for(i=0, minV=100000, minI=-1; i<nCachedModels; i++) {
						if(dissimil[i] > -1 && dissimil[i] < minV) {
							minV = dissimil[i];
							minI = i;
						}
					}
					if(minI == -1) {
						pout[j] = NA_INTEGER;
						break;
					} else {
						if(!IS_LOGICAL(VECTOR_ELT(_cachedModelList, minI))) {
							pout[j] = minI + 1;
							break;
						}
						dissimil[minI] = -1;
					}
				} while(true);

			} else
				pout[j] = minI + 1;
		}
	}
	
	free(dissimil);
	UNPROTECT(1);
	return out;
}
