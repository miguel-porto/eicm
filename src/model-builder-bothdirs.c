#include "eicm.h"

/**
 * MODEL BUILDER
 * From a vector of numbers, construct the model according to mask and offset
 * mask: a list of two integer matrices. Zeroes will not be estimated, non-zeroes will.
 * The structure is (imagine 5 species, 2 env of which 1 is to be estimated; and 10 samples):
 * 00000 EEEEE FFFFF IIII III II I JJJJ JJJ JJ J (LLLLLLLLLL)
 * 0: intercepts
 * E: environmental coefficients 1st env
 * E: environmental coefficients 2nd env
 * I: species interaction coefficients (lower triangle)
 * J: species interaction coefficients (upper triangle)
 * (L): environmental predictor - NOTE it is only added afterwards, in the R code, not here.
 *
 * NOTE: in this approach we don't estimate direction. We estimate coefficients of both directions;
 * and give them as a result. It's up to the user to know how to handle cyclicity.
 *
 * NOTE: this is the only approach that works with gradient-based optimizers, because all parameters have
 * continuous responses and they don't interact with each other during estimation.
**/
SEXP _makeModelMatricesFromPars(SEXP _pars, SEXP _nspecies, SEXP _nenv, SEXP _mask, SEXP _offset) {
	double *pars = 			NUMERIC_POINTER(_pars);
	int nspecies = 			INTEGER_POINTER(_nspecies)[0];
	int nenv = 				INTEGER_POINTER(_nenv)[0];
	double *offsetEnv 		= NULL;
	double *offsetSp		= NULL;
	int *maskEnv			= NULL;
	int *maskSp				= NULL;
	int constMaskEnv		= -1;
	int constMaskSp			= -1;
	SEXP out, names;
	double *envmat, *spmat;
	int i, j, nSpCoefs, counter = 0, iptr;
	
	PROTECT(out = 			NEW_LIST(2));
	PROTECT(names =			NEW_CHARACTER(2));
	SET_STRING_ELT(names, 0, mkChar("env"));
	SET_STRING_ELT(names, 1, mkChar("sp"));
	SET_NAMES(out, names);
	UNPROTECT(1);
	SET_ELEMENT(out, 0, Rf_allocMatrix(REALSXP, nspecies, nenv));
	SET_ELEMENT(out, 1, Rf_allocMatrix(REALSXP, nspecies, nspecies));
	
	envmat = NUMERIC_POINTER(VECTOR_ELT(out, 0));
	spmat = NUMERIC_POINTER(VECTOR_ELT(out, 1));
	
	if(_offset == R_NilValue)
		Rf_error("Offset must be always supplied - if not used, set to all zero.");
		
	// TODO get by names
	offsetEnv = NUMERIC_POINTER(VECTOR_ELT(_offset, 0));
	offsetSp = NUMERIC_POINTER(VECTOR_ELT(_offset, 1));
	
	if(_mask != R_NilValue) {	// TODO get by names
		if(LENGTH(VECTOR_ELT(_mask, 0)) > 1)
			maskEnv = INTEGER_POINTER(VECTOR_ELT(_mask, 0));
		else
			constMaskEnv = INTEGER_POINTER(VECTOR_ELT(_mask, 0))[0];
			
		if(LENGTH(VECTOR_ELT(_mask, 1)) > 1)
			maskSp = INTEGER_POINTER(VECTOR_ELT(_mask, 1));
		else
			constMaskSp = INTEGER_POINTER(VECTOR_ELT(_mask, 1))[0];
	}
	
	memset(spmat, 0, sizeof(double) * nspecies * nspecies);
	memset(envmat, 0, sizeof(double) * nspecies * nenv);
	
	if((_mask == R_NilValue) != (_offset == R_NilValue))
		Rf_error("If offset is supplied, mask must be supplied and vice-versa.");
	
	if(_mask == R_NilValue || constMaskEnv == 1) {
		// all env coefs
		// set the intercept
		for(i=0; i<nspecies; i++) {
			envmat[i] = pars[counter];
			counter ++;
		}
		
		// set environmental coefs
		for(i=0; i<nspecies; i++) {
			for(j=1; j<nenv; j++) {
				envmat[i + j * nspecies] = pars[counter];
				counter ++;
			}
		}
	} else if(maskEnv != NULL) {
		// set the intercept
		for(i=0; i<nspecies; i++) {
			if(maskEnv[i]) {
				envmat[i] = pars[counter];
				counter ++;
			} else
				envmat[i] = offsetEnv[i];
		}

		int ptr;
		// set environmental coefs
		for(i=0; i<nspecies; i++) {
			for(j=1; j<nenv; j++) {
				ptr = i + j * nspecies;
				if(maskEnv[ptr]) {
					envmat[ptr] = pars[counter];
					counter ++;
				} else
					envmat[ptr] = offsetEnv[ptr];
			}
		}
	} // else do not set env coefs, cause env mask is all zero
	
	// now species interactions
	if(_mask == R_NilValue || constMaskSp == 1) {
		// full interactions
		// set species interaction coefs
		nSpCoefs = nspecies * (nspecies - 1) / 2;
		for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {
			for(j=i + 1; j<nspecies; j++) {
				spmat[j + iptr] = pars[counter];		// LT
				spmat[i + j * nspecies] = pars[counter + nSpCoefs];		// UT
				counter ++;
			}
		}
		counter += nSpCoefs;
	} else if(maskSp != NULL) {
		// set species interaction coefs
		int ptrLT, ptrUT;
		for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {	// columns
			for(j=i + 1; j<nspecies; j++) {
				ptrLT = j + iptr;
// Rprintf("%d x %d: LT %d UT %d maskLT %d maskUT %d", j, i, ptrLT, ptrUT, maskSp[ptrLT], maskSp[ptrUT]);
				if(maskSp[ptrLT]) {
					spmat[ptrLT] = pars[counter];
					counter ++;
				} else 
					spmat[ptrLT] = offsetSp[ptrLT];

			}
		}
		for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {	// columns
			for(j=i + 1; j<nspecies; j++) {
				ptrUT = i + j * nspecies;
				if(maskSp[ptrUT]) {
					spmat[ptrUT] = pars[counter];
					counter ++;
				} else 
					spmat[ptrUT] = offsetSp[ptrUT];
			}
		}
	} // else do nothing; don't estimate species interactions
			
/*	if(LENGTH(_pars) != counter) {
		UNPROTECT(1);
		Rf_error("The number of parameters %d does not conform to the provided model frame (which requires %d).", LENGTH(_pars), counter);
	}*/

	UNPROTECT(1);
	return out;
}

SEXP _makeParsFromModelMatrices(SEXP _matrices, SEXP _mask) {
	double *envmat =		NUMERIC_POINTER(VECTOR_ELT(_matrices, 0));	// TODO get by names
	double *spmat =			NUMERIC_POINTER(VECTOR_ELT(_matrices, 1));
	int nspecies = 			Rf_nrows(VECTOR_ELT(_matrices, 1));
	int nenv = 				Rf_ncols(VECTOR_ELT(_matrices, 0));
	int *maskEnv			= NULL;
	int *maskSp				= NULL;
	int constMaskEnv		= -1;
	int constMaskSp			= -1;
	SEXP out;
	int i, j, nSpCoefs, counter = 0, iptr;
	int maxNPars =			nspecies * nenv + nspecies * (nspecies - 1);
	double *pars =			calloc(maxNPars, sizeof(double));
	
	if(_mask != R_NilValue) {	// TODO get by names
		if(LENGTH(VECTOR_ELT(_mask, 0)) > 1)
			maskEnv = INTEGER_POINTER(VECTOR_ELT(_mask, 0));
		else
			constMaskEnv = INTEGER_POINTER(VECTOR_ELT(_mask, 0))[0];
			
		if(LENGTH(VECTOR_ELT(_mask, 1)) > 1)
			maskSp = INTEGER_POINTER(VECTOR_ELT(_mask, 1));
		else
			constMaskSp = INTEGER_POINTER(VECTOR_ELT(_mask, 1))[0];
	}

	if(_mask == R_NilValue || constMaskEnv == 1) {
		// set the intercept
		for(i=0; i<nspecies; i++) {
			pars[counter] = envmat[i];
			counter ++;
		}
		
		// set environmental coefs
		for(i=0; i<nspecies; i++) {
			for(j=1; j<nenv; j++) {
				pars[counter] = envmat[i + j * nspecies];
				counter ++;
			}
		}
	} else if(maskEnv != NULL) {
		// set the intercept
		for(i=0; i<nspecies; i++) {
			if(maskEnv[i]) {
				pars[counter] = envmat[i];
				counter ++;
			}
		}

		int ptr;
		// set environmental coefs
		for(i=0; i<nspecies; i++) {
			for(j=1; j<nenv; j++) {
				ptr = i + j * nspecies;
				if(maskEnv[ptr]) {
					pars[counter] = envmat[ptr];
					counter ++;
				}
			}
		}
	} // else do nothing
	
	if(_mask == R_NilValue || constMaskSp == 1) {
		// set species interaction coefs
		nSpCoefs = nspecies * (nspecies - 1) / 2;
		for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {
			for(j=i + 1; j<nspecies; j++) {
				pars[counter] = spmat[j + iptr];
				pars[counter + nSpCoefs] = spmat[i + j * nspecies];
				counter ++;
			}
		}
		counter += nSpCoefs;
	} else if(maskSp != NULL) {
		// set species interaction coefs
		// first count how many coefs will be estimated
		int ptrLT, ptrUT;
		for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {	// columns
			for(j=i + 1; j<nspecies; j++) {
				ptrLT = j + iptr;
				if(maskSp[ptrLT]) {
					pars[counter] = spmat[ptrLT];
					counter ++;
				}
			}
		}
		for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {	// columns
			for(j=i + 1; j<nspecies; j++) {
				ptrUT = i + j * nspecies;
				if(maskSp[ptrUT]) {
					pars[counter] = spmat[ptrUT];
					counter ++;
				}
			}
		}
	}
	
	PROTECT(out = NEW_NUMERIC(counter));
	double *pout = NUMERIC_POINTER(out);
	memcpy(pout, pars, counter * sizeof(double));
	UNPROTECT(1);

	return out;
}

SEXP _getNumberOfParameters(SEXP _nspecies, SEXP _nenv, SEXP _mask) {
	int nspecies = 			INTEGER_POINTER(_nspecies)[0];
	int nenv = 				INTEGER_POINTER(_nenv)[0];
	int *maskEnv			= NULL;
	int *maskSp				= NULL;
	int constMaskEnv		= -1;
	int constMaskSp			= -1;
	SEXP out;
	int i, j, counter = 0, iptr;
	PROTECT(out = 			NEW_INTEGER(1));
	int *pout =				INTEGER_POINTER(out);	

	if(_mask != R_NilValue) {	// TODO get by names
		if(LENGTH(VECTOR_ELT(_mask, 0)) > 1)
			maskEnv = INTEGER_POINTER(VECTOR_ELT(_mask, 0));
		else
			constMaskEnv = INTEGER_POINTER(VECTOR_ELT(_mask, 0))[0];
			
		if(LENGTH(VECTOR_ELT(_mask, 1)) > 1)
			maskSp = INTEGER_POINTER(VECTOR_ELT(_mask, 1));
		else
			constMaskSp = INTEGER_POINTER(VECTOR_ELT(_mask, 1))[0];
	}

		
	if(_mask == R_NilValue || (constMaskEnv == 1 && constMaskSp == 1)) {
		// full model
		pout[0] = nspecies * nenv + nspecies * (nspecies - 1);
	} else {	// we've got mask and offset
		if(maskEnv != NULL) {
			for(i=0; i<nspecies; i++) {
				for(j=0; j<nenv; j++) {
					if(maskEnv[i + j * nspecies]) counter ++;
				}
			}
		} else if(constMaskEnv == 1)
			counter += nspecies * nenv;
		
		if(maskSp != NULL) {
			int ptrLT, ptrUT;
			for(i=0, iptr=0; i<nspecies; i++, iptr += nspecies) {
				for(j=i + 1; j<nspecies; j++) {
					ptrLT = j + iptr;
					ptrUT = i + j * nspecies;

					if(maskSp[ptrLT] || maskSp[ptrUT])	// estimate the coefficient
						counter ++;
					
					if(maskSp[ptrLT] && maskSp[ptrUT])	// normal case, estimate the direction also
						counter ++;
				}	
			}
		} else if(constMaskSp == 1)
			counter += nspecies * (nspecies - 1);
			
		pout[0] = counter;
	}
	
	UNPROTECT(1);
	return out;
}

