#include "eicm.h"
//unsigned int invLogitTab[LOOKUPTABLESIZE];

void createInverseLinkFunctionTable(unsigned int *table) {
	for(int i=0; i<LOOKUPTABLESIZE; i++)
		table[i] = 1 / (1 + exp(-(double) (i - MIDDLELOOKUPTABLE) / DIVISOR)) * MAXRND;
}

void createInverseLinkFunctionTableProb(float *table) {
	for(int i=0; i<LOOKUPTABLESIZE; i++)
		table[i] = (float) (1 / (1 + exp(-(double) (i - MIDDLELOOKUPTABLE) / DIVISOR)));
}

/**
* Matrix product with transposed m2
* TODO: this can be optimized for speed, see https://github.com/deuxbot/fast-matrix-multiplication
* Usual parameters:
* m1: sample x env
* m2: species x env
* m1rows: nº samples
* m2rows: nº species
* m1cols: nº env vars
*
* Returns the env LP in a sample x species matrix
*/
void matProd(const double *m1, const double *m2, double *out, const int m1rows, const int m2rows, const int m1cols) {
	for(int i = 0, icol = 0; i < m2rows; i++, icol += m1rows) {
		for(int j = 0; j < m1rows; j++) {
			double value = 0;
			for(int k = 0; k < m1cols; k++)
				value += m2[i + k * m2rows] * m1[j + k * m1rows];
			out[j + icol] = value;
//			printf("%f ",value);
		}
	}
}

/**
* Same as matProd but skips species.
*/
void matProdExclSpecies(const double *m1, const double *m2, double *out, const int m1rows, const int m2rows, const int m1cols, const bool *exclSpecies) {
	for(int i = 0, icol = 0; i < m2rows; i++, icol += m1rows) {
		if(exclSpecies[i]) continue;
		for(int j = 0; j < m1rows; j++) {
			double value = 0;
			for(int k = 0; k < m1cols; k++)
				value += m2[i + k * m2rows] * m1[j + k * m1rows];
			out[j + icol] = value;
//			printf("%f ",value);
		}
	}
}

/**
* This is a special partial matrix prodcut, which adds the result to the output matrix.
* It is only used for one single purpose.
* NOTE: whichm1cols is 1-based!
*/
void matProdInt(int *m1, double *m2, double *out, int m1rows, int m2rows, int m1cols
	, int *whichm1cols, int lenwhichm1cols, int *whichm2rows, int lenwhichm2rows) {
	for(int i = 0; i < lenwhichm2rows; i++) {
		int icol = (whichm2rows[i] - 1) * m1rows;
		for(int j = 0; j < m1rows; j++) {
			double value = 0;
			for(int k = 0; k < lenwhichm1cols; k++) {
				int col = whichm1cols[k] - 1;
				value += m2[(whichm2rows[i] - 1) + col * m2rows] * (double) m1[j + col * m1rows];
			}
			out[j + icol] += value;
		}
	}
}

/**
* This is a special partial matrix prodcut, which adds the result to the output matrix.
* It is only used for one single purpose in prediction.
* Usual parameters:
* NOTE that environment has been replaced by specie here.
* m1: presences as sample x species
* m2: species coefs as species x species
* NOTE: whichm1cols is 1-based!
*/
void matProdShort(const short *m1, const double *m2, double *out, const int m1rows, const int m2rows, const int m1cols
	, const int *whichm1cols, const int lenwhichm1cols, const int *whichm2rows, const int lenwhichm2rows) {
	for(int i = 0; i < lenwhichm2rows; i++) {
		int icol = (whichm2rows[i] - 1) * m1rows;
		for(int j = 0; j < m1rows; j++) {
			double value = 0;
			for(int k = 0; k < lenwhichm1cols; k++) {
				int col = whichm1cols[k] - 1;
				value += m2[(whichm2rows[i] - 1) + col * m2rows] * (double) m1[j + col * m1rows];
			}
			out[j + icol] += value;
		}
	}
}


/**
* This is a special partial matrix prodcut, which adds the result to the output matrix, and allows excluding samples from the computation
* It is only used for one single purpose.
*/
void matProdShortExclSamp(short *prM, double *coM, double *out, int prMrows, int coMrows, int prMcols
	, int *whichprMcols, int lenwhichprMcols, int *whichcoMrows, int lenwhichcoMrows, short *excludeSamples, short maxErrors) {
	for(int i = 0; i < lenwhichcoMrows; i++) {
		int icol = (whichcoMrows[i] - 1) * prMrows;
		for(int j = 0; j < prMrows; j++) {
			if(excludeSamples[j] >= maxErrors) continue;
			double value = 0;
			for(int k = 0; k < lenwhichprMcols; k++) {
				int col = whichprMcols[k] - 1;
				value += coM[(whichcoMrows[i] - 1) + col * coMrows] * (double) prM[j + col * prMrows];
			}
			out[j + icol] += value;
		}
	}
}
/*
void matProdShortExclSamp2(short *prM, double *coM, double *out, int prMrows, int coMrows, int prMcols
	, bool *filledSpecies, int *whichcoMrows, int lenwhichcoMrows, bool *excludeSamples) {
//	for(int i=0;i<lenwhichcoMrows; i++) printf("%d ",whichcoMrows[i]);
	for(int i = 0; i < lenwhichcoMrows; i++) {
		int icol = (whichcoMrows[i] - 1) * prMrows;
		for(int j = 0; j < prMrows; j++) {
			if(excludeSamples[j]) continue;
			double value = 0;
			for(int k = 0; k < prMcols; k++) {
				if(!filledSpecies[k]) continue;
				value += coM[(whichcoMrows[i] - 1) + k * coMrows] * (double) prM[j + k * prMrows];
			}
			out[j + icol] += value;
		}
	}
}
*/

