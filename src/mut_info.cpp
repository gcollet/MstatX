/*
 *  mut_info.cpp
 *  MstatX++
 *
 *  Created by Guillaume Collet on 27/11/10.
 *
 */

#include "mut_info.h"

#include <iostream>
#include <string>

#define MIN(x,y)  (x < y ? x : y)

using namespace std;

/** Destructor */
MIStat :: ~MIStat()
{
	if (f_i != NULL){
		for (int i(0); i < ncol; ++i) {
			if (f_i[i] != NULL){
				free (f_i[i]);	
			}
		}
		free(f_i);
	}
	if (f_ij != NULL){
		for (int c(0); c < npair; ++c) {
			if (f_ij[c] != NULL){
				for(int i(0); i < alph_size; ++i){
					if (f_ij[c][i] != NULL){
						free(f_ij[c][i]);
					}
				}
				free (f_ij[c]);	
			}
		}		
		free(f_ij);
	}
}


/* Calculate the weight of sequence i in the msa 
 * by the formula from Henikoff & Henikoff (1994)
 * w_i = \frac{1}{L}\sum_{x=1}^{L}\frac{1}{k_x n_{x_i}}
 * We use these notations in the code below
 */
float 
MIStat :: calcSeqWeight(Msa & msa, int i)
{
	int x, seq;
	int k;								    /**< number of symbol types in a column */
	int n;								    /**< number of occurence of aa[i][x] in column x */
	int L = msa.getNcol();    /**< number of columns (i.e. length of the alignment) */
	int nseq = msa.getNseq(); /**< number of rows (i.e. number of sequences in the alignment) */
	float w;						    	/**< weight of sequence i */
	
	w = 0.0;
	for	(x = 0; x < L; ++x){
		k = msa.getNtype(x);
		/* Calculate the number of aa[i][x] in current column */
		n = 0;
		for(seq = 0; seq < nseq; ++seq){
		  if (msa.getSymbol(i, x) == msa.getSymbol(seq, x)){
				n++;
			}
		}
		w += (float) 1 / (float) (n * k); 
	}
	w /= (float) L;
	return w ;
}


void
MIStat :: calculateStatistic(Msa & msa)
{
	/* Init size */
	ncol = msa.getNcol();
	npair = ncol * (ncol - 1) / 2;
	nseq = msa.getNseq();
	string alphabet = msa.getAlphabet();
	alph_size = alphabet.size();
	
	f_i = NULL;
	f_ij = NULL;
	cout << alphabet << "\n";
	cout << "ncol = " << ncol << " -- nseq = " << nseq << " -- alph_size = " << alph_size << "\n";
	
	/* Allocate f_i (size = alph_size * ncol) */
	f_i = (float **) calloc(ncol, sizeof(float*));
	if (f_i == NULL){
		fprintf(stderr,"Cannot Allocate frequency matrix (f_i)\n");
		exit(0);
	}
	for (int i(0); i < ncol; i++){
		f_i[i] = (float *) calloc(alphabet.size(), sizeof(float));
		if (f_i[i] == NULL){
			fprintf(stderr,"Cannot Allocate frequency matrix (f_i)\n");
			exit(0);
		}
	}
	
	/* Allocate f_ij (size = alph_size * alph_size * npair) */
	f_ij = (float ***) calloc(npair, sizeof(float**));
	if (f_ij == NULL){
		fprintf(stderr,"Cannot Allocate pair frequency matrix (f_ij)\n");
		exit(0);
	}
	for (int c(0); c < npair; c++){
		f_ij[c] = (float **) calloc(alph_size, sizeof(float*));
		if (f_ij[c] == NULL){
			fprintf(stderr,"Cannot Allocate pair frequency matrix (f_ij)\n");
			exit(0);
		}
		for (int i(0); i < alph_size; i++){
			f_ij[c][i] = (float *) calloc(alph_size, sizeof(float));
			if (f_ij[c][i] == NULL){
				fprintf(stderr,"Cannot Allocate pair frequency matrix (f_ij)\n");
				exit(0);
			}
		}
	}
	
	
	/* Calculate Sequence Weights */
	cout << "Calculate Sequence Weights\n";
	seq_weight = vector<float>(nseq,0.0);
	int seq;
#pragma omp parallel for schedule (dynamic)
	for (seq = 0; seq < nseq; ++seq){
		seq_weight[seq] = calcSeqWeight(msa,seq);
	}
	
	
	/* Calculate frequencies for single positions */
	for (int x(0); x < ncol; x++){
		float total = 0.0;
		for (int a(0); a < alphabet.size(); a++){
			f_i[x][a] = 0.0;
			for (int j(0); j < nseq; j++){
				if(msa.getSymbol(j, x) == alphabet[a]){
					f_i[x][a] += seq_weight[j];
				}
			}
			if (f_i[x][a] == 0.0) {
				f_i[x][a] = (float) 1 / (float) nseq;
			}
		}
	}
	
	
	
}

